/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_MMA.cpp
 *
 * Native port of Svanberg's GCMMA (2002). See header for provenance. The numerical
 * arithmetic mirrors github.com/kkmaute/gcmma line-for-line; func()/grad() replace the
 * TPL's opt_alg_gcmma_func_wrap/grad_wrap with direct calls to the inherited MORIS
 * criteria evaluation, so MORIS owns the optimization loop.
 */

#include <cmath>
#include <vector>
#include <string>
#include <cstdio>

#include "cl_OPT_Algorithm_MMA.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris::opt
{
    //------------------------------------------------------------------------------

    Algorithm_MMA::Algorithm_MMA( const Parameter_List& aParameterList )
            : mMaxInnerIter( aParameterList.get< sint >( "max_inner_its" ) )
            , mNormDrop( aParameterList.get< real >( "norm_drop" ) )
            , mAsympAdapt0( aParameterList.get< real >( "asymp_adapt0" ) )
            , mAsympShrink( aParameterList.get< real >( "asymp_adaptb" ) )
            , mAsympExpand( aParameterList.get< real >( "asymp_adaptc" ) )
            , mStepSize( aParameterList.get< real >( "step_size" ) )
            , mPenalty( aParameterList.get< real >( "penalty" ) )
    {
        mRestartIndex         = aParameterList.get< sint >( "restart_index" );
        mMaxIterationsInitial = aParameterList.get< sint >( "max_its" );
        mMaxIterations        = mMaxIterationsInitial;
        mReinitHook           = aParameterList.get< bool >( "reinit_in_place" );
        mReinitOnAcceptOnly   = aParameterList.get< bool >( "reinit_on_accept_only" );
        mDebugInnerTrials     = aParameterList.get< bool >( "debug_inner_trials" );

        MORIS_LOG_INFO( "Native MMA initialized: max_its=%d, step_size=%.6f, penalty=%.2f%s%s%s",
                (int)mMaxIterations, mStepSize, mPenalty, mReinitHook ? ", in-place reinit hook ON" : "",
                mReinitHook && mReinitOnAcceptOnly ? ", redistance on accepted iterates only" : "",
                mDebugInnerTrials ? ", inner-trial diagnostics ON" : "" );
        MORIS_ERROR( !mReinitOnAcceptOnly || mReinitHook, "Native MMA: reinit_on_accept_only needs reinit_in_place = true." );
    }

    //------------------------------------------------------------------------------

    Algorithm_MMA::~Algorithm_MMA() {}

    //------------------------------------------------------------------------------

    uint
    Algorithm_MMA::solve(
            uint                                   aCurrentOptAlgInd,
            std::shared_ptr< moris::opt::Problem > aOptProb )
    {
        Tracer tTracer( "OptimizationAlgorithm", "MMA", "Solve" );

        mRunning          = opt::Task::wait;
        mCurrentOptAlgInd = aCurrentOptAlgInd;
        mProblem          = aOptProb;

        if ( mRestartIndex > 0 )
        {
            gLogger.set_opt_iteration( mRestartIndex );
        }

        if ( par_rank() == 0 )
        {
            this->mma_solve();

            mRunning = opt::Task::exit;
            this->communicate_running_status();
        }
        else
        {
            mPrint = false;
            this->dummy_solve();
        }

        uint tOptIter = gLogger.get_opt_iteration();
        gLogger.set_iteration( "OPT", "Manager", "Perform", tOptIter );

        return tOptIter;
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::func( double* aX, real& aF0, std::vector< double >& aF )
    {
        // Guard: the subproblem can overflow on a degenerate evaluation (huge objective/
        // gradient jump from a near-degenerate XFEM cut) and produce a non-finite trial
        // design. Evaluating it would compute garbage criteria and trip the gradient-
        // consistency assert. Request an optimizer restart instead (the solve loop checks
        // the flag after every func call and exits; the Manager then re-initializes from
        // the current design) and return the last evaluated objective and constraints.
        for ( int i = 0; i < mNumVar; ++i )
        {
            if ( !std::isfinite( aX[ i ] ) )
            {
                MORIS_LOG_WARNING( "Native MMA produced a non-finite design vector; requesting optimizer restart." );

                mProblem->request_restart();

                aF0 = mObjScale * this->get_objectives()( 0 );

                const Matrix< DDRMat >& tLastCon = this->get_constraints();
                for ( int j = 0; j < mNumCon; ++j ) aF[ j ] = tLastCon( j );

                return;
            }
        }

        Vector< real > tADVs( mNumVar );
        for ( int i = 0; i < mNumVar; ++i ) tADVs( i ) = aX[ i ];

        // Write the restart file for this design, as the TPL GCMMA wrapper does before every
        // criteria evaluation (MB-016: the native port wrote none, so no warm restart was possible).
        // Capture the index first: compute_design_criteria() advances the optimization iteration (MB-021)
        const uint tOptIter = gLogger.get_opt_iteration();
        this->write_advs_to_file( tADVs, tOptIter );

        this->compute_design_criteria( tADVs );

        aF0 = mObjScale * this->get_objectives()( 0 );

        const Matrix< DDRMat >& tCon = this->get_constraints();
        for ( int j = 0; j < mNumCon; ++j ) aF[ j ] = tCon( j );

        // In-place reinitialization hook: the workflow may have redistanced the design (overwritten
        // tADVs with the SDF coefficients) during criteria evaluation. Store those separately rather
        // than writing them back into aX -- this keeps xmma as the RAW subproblem solution so that
        // raaupdate() and the conservatism check operate on the point the subproblem actually
        // produced. The redistanced design is committed to xval after the step is accepted; the
        // criteria call already set mADVs = tADVs, so the subsequent gradient call stays consistent.
        if ( mReinitHook )
        {
            mLastRedistanced.assign( mNumVar, 0.0 );
            for ( int i = 0; i < mNumVar; ++i ) mLastRedistanced[ i ] = tADVs( i );

            // rewrite this iteration's restart file with the redistanced design the criteria were evaluated at, so
            // checkpoints (and restarts) hold phi = R(psi) rather than the raw trial point psi (written above).
            // Name it from the index captured before the evaluation: the logger has already advanced to the next
            // iteration, whose file the next trial point would overwrite (MB-021)
            this->write_advs_to_file( tADVs, tOptIter );
        }
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::grad( const double* aX )
    {
        Vector< real > tADVs( mNumVar );
        for ( int i = 0; i < mNumVar; ++i ) tADVs( i ) = aX[ i ];

        this->compute_design_criteria_gradients( tADVs );

        const Matrix< DDRMat >& tObjGrad = this->get_objective_gradients();      // 1 x numvar
        const Matrix< DDRMat >& tConGrad = this->get_constraint_gradients();     // numcon x numvar

        // Uniform objective scaling: engage only above the pathology threshold, quantized to
        // decades, and re-book the stored objective value (currently in the previous scale)
        // so the whole outer iteration is consistent in the new scale.
        real tRawMax = 0.0;
        for ( int i = 0; i < mNumVar; ++i ) tRawMax = mx( tRawMax, std::abs( tObjGrad( 0, i ) ) );

        real tNewScale = 1.0;
        if ( tRawMax > sObjGradEngage )
        {
            // smallest power of ten that brings tRawMax at or below the target
            tNewScale = std::pow( 10.0, -std::ceil( std::log10( tRawMax / sObjGradTarget ) ) );
        }

        if ( tNewScale != mObjScale )
        {
            MORIS_LOG_INFO( "Native MMA: objective scale -> %.1e (raw max|df0dx| = %.4e).",
                    tNewScale, tRawMax );
            f0val    *= tNewScale / mObjScale;
            // in-place redistance mode: keep the conservatism weight consistent with the objective's scale,
            // so a scale change does not jump the approximation's curvature (raa0 is in objective units)
            if ( mReinitHook )
            {
                raa0 *= tNewScale / mObjScale;
            }
            mObjScale = tNewScale;
        }

        for ( int i = 0; i < mNumVar; ++i ) df0dx[ i ] = mObjScale * tObjGrad( 0, i );
        for ( int j = 0; j < mNumCon; ++j )
            for ( int i = 0; i < mNumVar; ++i ) dfdx[ j ][ i ] = tConGrad( j, i );
    }

    //------------------------------------------------------------------------------

    bool
    Algorithm_MMA::accept_only_redistance() const
    {
        // without an inner loop there is no conservatism test, so evaluating the trial at R(psi) is already consistent
        // and the cheapest ordering is the old one (one evaluation per iteration)
        return mReinitHook && mReinitOnAcceptOnly && mMaxInnerIter > 0;
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::trial_func( int aIter, int aInner, int& aNumFunc )
    {
        const bool tRawTrial = this->accept_only_redistance();

        real                  tAltF0 = 0.0;
        std::vector< double > tAltF( mNumCon, 0.0 );
        real                  tDeltaR  = 0.0;
        bool                  tHaveAlt = false;

        // diagnostic evaluation of the OTHER design first (raw psi if the trial is redistanced, R(psi) if it is raw), so
        // that the Problem afterwards holds the state of the real trial evaluation. No restart file is written for it and
        // mLastRedistanced is not touched. It advances the evaluation counter (decks that ramp by evaluation must pin it).
        if ( mDebugInnerTrials && mReinitHook )
        {
            mProblem->set_reinit_suppressed( !tRawTrial );

            Vector< real > tADVs( mNumVar );
            for ( int i = 0; i < mNumVar; ++i ) tADVs( i ) = xmma[ i ];

            this->compute_design_criteria( tADVs );
            aNumFunc++;

            if ( !mProblem->restart_optimization() )
            {
                tAltF0 = mObjScale * this->get_objectives()( 0 );
                const Matrix< DDRMat >& tCon = this->get_constraints();
                for ( int j = 0; j < mNumCon; ++j ) tAltF[ j ] = tCon( j );
                for ( int i = 0; i < mNumVar; ++i ) tDeltaR = mx( tDeltaR, std::abs( tADVs( i ) - xmma[ i ] ) );
                tHaveAlt = true;
            }
        }

        // the real trial evaluation: raw psi in accept-only mode, otherwise the workflow redistances it as before
        if ( mReinitHook )
        {
            mProblem->set_reinit_suppressed( tRawTrial );
        }
        this->func( xmma.data(), f0valnew, fvalnew );
        aNumFunc++;

        if ( mReinitHook && !tRawTrial )
        {
            for ( int i = 0; i < mNumVar && (int)mLastRedistanced.size() == mNumVar; ++i )
            {
                tDeltaR = mx( tDeltaR, std::abs( mLastRedistanced[ i ] - xmma[ i ] ) );
            }
        }

        if ( mDebugInnerTrials && !mProblem->restart_optimization() )
        {
            this->log_trial( aIter, aInner, tAltF0, tAltF, tDeltaR, tHaveAlt );
        }
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::log_trial( int aIter, int aInner, real aAltF0, const std::vector< double >& aAltF, real aDeltaR, bool aHaveAlt ) const
    {
        // what the conservatism test compares: the approximation (f0app, fapp) at the subproblem solution psi = xmma
        // against the evaluated criteria (f0valnew, fvalnew); "at" names the design the evaluation used
        const bool  tRawTrial = this->accept_only_redistance() || !mReinitHook;
        std::string tLine     = "MMA-TRIAL iter=" + std::to_string( aIter ) + " inner=" + std::to_string( aInner )
                          + " at=" + ( tRawTrial ? "psi" : "R(psi)" ) + " alt=" + ( aHaveAlt ? ( tRawTrial ? "R(psi)" : "psi" ) : "none" );

        char tBuf[ 512 ];
        std::snprintf( tBuf, sizeof( tBuf ), " |R(psi)-psi|inf=%.6e conservative=%d | obj app=%.12e true=%.12e alt=%.12e raa0=%.6e",
                aDeltaR, this->is_conservative() ? 1 : 0, f0app, f0valnew, aHaveAlt ? aAltF0 : std::nan( "" ), raa0 );
        tLine += tBuf;
        for ( int j = 0; j < mNumCon; ++j )
        {
            std::snprintf( tBuf, sizeof( tBuf ), " | c%d app=%.12e true=%.12e alt=%.12e raa=%.6e",
                    j, fapp[ j ], fvalnew[ j ], aHaveAlt ? aAltF[ j ] : std::nan( "" ), raa[ j ] );
            tLine += tBuf;
        }
        MORIS_LOG_INFO( "%s", tLine.c_str() );
    }

    //------------------------------------------------------------------------------

    bool
    Algorithm_MMA::all_finite( real aF0, const std::vector< double >& aF ) const
    {
        if ( !std::isfinite( aF0 ) ) return false;
        for ( int j = 0; j < mNumCon; ++j )
            if ( !std::isfinite( aF[ j ] ) ) return false;
        return true;
    }

    //------------------------------------------------------------------------------

    bool
    Algorithm_MMA::criteria_finite() const
    {
        return this->all_finite( f0valnew, fvalnew );
    }

    //------------------------------------------------------------------------------

    bool
    Algorithm_MMA::is_conservative() const
    {
        if ( ( f0app + epsimin ) < f0valnew ) return false;
        for ( int j = 0; j < mNumCon; ++j )
            if ( ( fapp[ j ] + epsimin ) < fvalnew[ j ] ) return false;
        return true;
    }

    //------------------------------------------------------------------------------

    Algorithm_MMA::Trial_Status
    Algorithm_MMA::backtrack_non_finite( int aIter, int& aNumFunc )
    {
        // NaN guard: a non-finite trial criterion means the forward evaluation failed -- almost
        // always an XTK decomposition failure where the moving interface's cut spans two HMR
        // refinement levels ("not all child meshes on the same level"), which poisons every IQI
        // integral (including the geometric volume/perimeter). Reject the step and shrink the
        // move limit, re-solving until the criteria are finite. As the move limit -> 0 the trial
        // design -> xval, which was finite when committed.
        int       tBacktrack    = 0;
        const int tMaxBacktrack = 15;
        while ( !this->criteria_finite() && tBacktrack < tMaxBacktrack )
        {
            tBacktrack++;
            mMoveScale *= 0.5;
            MORIS_LOG_INFO( "Native MMA iter %d: non-finite criteria -- likely XTK decomposition failure "
                            "(interface cut spans different HMR refinement levels). Rejecting step; "
                            "move limit x%.4g (backtrack %d/%d).",
                    aIter, mMoveScale, tBacktrack, tMaxBacktrack );
            this->gcmmasub( aIter );
            this->trial_func( aIter, -tBacktrack, aNumFunc );

            if ( mProblem->restart_optimization() )
            {
                MORIS_LOG_INFO( "Native MMA iter %d: optimization restart requested during backtrack; exiting solve loop.", aIter );
                return Trial_Status::restart;
            }
        }

        if ( !this->criteria_finite() )
        {
            MORIS_LOG_INFO( "Native MMA iter %d: criteria still non-finite after %d backtracks.", aIter, tMaxBacktrack );
            return Trial_Status::failed;
        }

        // keep the (possibly shrunk) move limit for the rest of this outer iteration;
        // it is reset to 1.0 at the top of the next one.
        return Trial_Status::finite;
    }

    //------------------------------------------------------------------------------

    bool
    Algorithm_MMA::restore_accepted_state( int aIter, int& aNumFunc )
    {
        MORIS_LOG_INFO( "Native MMA iter %d: rejecting the iteration; re-evaluating the accepted design to restore "
                        "the problem state.",
                aIter );

        // evaluate into scratch criteria: the optimizer's own f0val/fval for xval stay as committed (in reinit mode
        // the restore evaluates R(xval), so the Problem holds that redistance with criteria consistent with it)
        real                  tF0 = 0.0;
        std::vector< double > tF( mNumCon, 0.0 );
        if ( mReinitHook )
        {
            // accept-only mode: restore exactly the committed design (it is already redistanced)
            mProblem->set_reinit_suppressed( this->accept_only_redistance() );
        }
        this->func( xval.data(), tF0, tF );
        aNumFunc++;

        if ( mProblem->restart_optimization() )
        {
            MORIS_LOG_INFO( "Native MMA iter %d: optimization restart requested while restoring; exiting solve loop.", aIter );
            return false;
        }

        // the accepted design evaluated finite when it was committed; failing now means the evaluation is not
        // reproducible. Fail loudly (nonzero exit for the run sentinels) rather than iterate on an unevaluable design
        MORIS_ERROR( this->all_finite( tF0, tF ),
                "Native MMA iter %d: the accepted design no longer evaluates to finite criteria; restart from the last "
                "checkpoint.",
                aIter );

        return true;
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::mma_initialize()
    {
        mNumVar = (int)mProblem->get_num_advs();
        mNumCon = (int)mProblem->get_num_constraints();

        epsimin = std::sqrt( (double)( mNumCon + mNumVar ) ) * 1.0e-9;
        raa0eps = 1.0e-7;
        a0      = 1.0;
        f0app   = 0.0;

        auto allocV = [ & ]( std::vector< double >& v, int n ) { v.assign( n, 0.0 ); };

        // size numvar
        for ( auto* v : { &xval, &xold1, &xold2, &xmax, &xmin, &low, &upp, &alfa, &beta,
                      &p0, &q0, &df0dx, &x, &xold, &delx, &diagx, &diagxinv, &xmma,
                      &xsi, &xsiold, &dxsi, &eta, &etaold, &deta } )
            allocV( *v, mNumVar );
        allocV( dx, mNumVar + 1 );

        // size numcon
        for ( auto* v : { &fval, &fvalnew, &a, &b, &c, &d, &r, &raa, &fapp, &gvec,
                      &y, &yold, &dy, &dely, &diagy, &diagyinv,
                      &lamold, &dellam, &dellamyi, &diaglam, &diaglamyi, &diaglamyiinv,
                      &mu, &muold, &dmu, &s, &sold, &ds } )
            allocV( *v, mNumCon );
        allocV( lam, mNumCon );
        allocV( dlam, mNumCon + 1 );

        // 2D numcon x numvar
        P.assign( mNumCon, std::vector< double >( mNumVar, 0.0 ) );
        Q.assign( mNumCon, std::vector< double >( mNumVar, 0.0 ) );
        GG.assign( mNumCon, std::vector< double >( mNumVar, 0.0 ) );
        dfdx.assign( mNumCon, std::vector< double >( mNumVar, 0.0 ) );

        // dense Newton system, k = min(numcon,numvar)+1
        int k = ( mNumCon < mNumVar ? mNumCon : mNumVar ) + 1;
        AA.assign( k * k, 0.0 );
        bb.assign( k, 0.0 );

        // initial values (get_advs/bounds return Vector< real >&)
        const Vector< real >& tAdv = mProblem->get_advs();
        const Vector< real >& tUpp = mProblem->get_upper_bounds();
        const Vector< real >& tLow = mProblem->get_lower_bounds();

        for ( int i = 0; i < mNumVar; ++i )
        {
            xval[ i ]  = tAdv( i );
            xmax[ i ]  = tUpp( i );
            xmin[ i ]  = tLow( i );
            xold1[ i ] = xval[ i ];
            xold2[ i ] = xval[ i ];
            low[ i ]   = xmin[ i ];
            upp[ i ]   = xmax[ i ];
            alfa[ i ]  = xmin[ i ];
            beta[ i ]  = xmax[ i ];
        }

        for ( int j = 0; j < mNumCon; ++j )
        {
            a[ j ] = 0.0;
            d[ j ] = 0.0;
        }
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::asymp( int iter )
    {
        if ( iter > 1.5 )
        {
            raa0 = mx( raa0eps, 0.1 * raa0 );
            for ( int j = 0; j < mNumCon; ++j ) raa[ j ] = mx( raa0eps, 0.1 * raa[ j ] );
        }

        if ( iter < 2.5 )
        {
            for ( int i = 0; i < mNumVar; ++i )
            {
                low[ i ] = xval[ i ] - mAsympAdapt0 * ( xmax[ i ] - xmin[ i ] );
                upp[ i ] = xval[ i ] + mAsympAdapt0 * ( xmax[ i ] - xmin[ i ] );
            }
        }
        else
        {
            for ( int i = 0; i < mNumVar; ++i )
            {
                double factor = 1.0;
                double xxx    = ( xval[ i ] - xold1[ i ] ) * ( xold1[ i ] - xold2[ i ] );

                if ( xxx > 0 ) factor = mAsympExpand;
                if ( xxx < 0 ) factor = mAsympShrink;

                low[ i ] = xval[ i ] - factor * ( xold1[ i ] - low[ i ] );
                upp[ i ] = xval[ i ] + factor * ( upp[ i ] - xold1[ i ] );

                double lowmin = xval[ i ] - 10.0 * ( xmax[ i ] - xmin[ i ] );
                double lowmax = xval[ i ] - 0.01 * ( xmax[ i ] - xmin[ i ] );
                double uppmin = xval[ i ] + 0.01 * ( xmax[ i ] - xmin[ i ] );
                double uppmax = xval[ i ] + 10.0 * ( xmax[ i ] - xmin[ i ] );

                low[ i ] = mx( low[ i ], lowmin );
                low[ i ] = mn( low[ i ], lowmax );
                upp[ i ] = mn( upp[ i ], uppmax );
                upp[ i ] = mx( upp[ i ], uppmin );
            }
        }

        // diagnostic: asymptote spread (upp-low)/range per variable. The floor is 0.02*range
        // (both asymptotes at their 0.01*range clamp). Variables collapsed to the floor have their
        // move limit choked to ~1% of range -- effectively locked. Reports how many are locked.
        {
            real tMinSp = 1.0e30, tMaxSp = 0.0, tSumSp = 0.0;
            int  tLocked = 0;
            for ( int i = 0; i < mNumVar; ++i )
            {
                real tSp = ( upp[ i ] - low[ i ] ) / ( xmax[ i ] - xmin[ i ] );
                tMinSp = mn( tMinSp, tSp );
                tMaxSp = mx( tMaxSp, tSp );
                tSumSp += tSp;
                if ( tSp < 0.05 ) tLocked++;
            }
            mNumLocked = tLocked;
            MORIS_LOG_INFO( "MMA asymptote spread iter %d: min=%.4f mean=%.4f max=%.4f | %d/%d vars at floor (locked)",
                    iter, tMinSp, tSumSp / (real)mNumVar, tMaxSp, tLocked, mNumVar );
        }
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::gcmmasub( int iter )
    {
        // move limit honors mMoveScale (shrunk by the NaN guard while backtracking off a bad step)
        const real tStep = mStepSize * mMoveScale;
        for ( int i = 0; i < mNumVar; ++i )
        {
            alfa[ i ] = mx( xmin[ i ], low[ i ] + 0.1 * ( xval[ i ] - low[ i ] ) );
            alfa[ i ] = mx( alfa[ i ], xval[ i ] - tStep * ( xmax[ i ] - xmin[ i ] ) );

            beta[ i ] = mn( xmax[ i ], upp[ i ] - 0.1 * ( upp[ i ] - xval[ i ] ) );
            beta[ i ] = mn( beta[ i ], xval[ i ] + tStep * ( xmax[ i ] - xmin[ i ] ) );
        }

        r0 = f0val;
        for ( int j = 0; j < mNumCon; ++j ) r[ j ] = fval[ j ];

        for ( int i = 0; i < mNumVar; ++i )
        {
            double ux1 = upp[ i ] - xval[ i ];
            double ux2 = ux1 * ux1;
            double xl1 = xval[ i ] - low[ i ];
            double xl2 = xl1 * xl1;
            double ul1 = upp[ i ] - low[ i ];

            double uxinv = 1.0 / ux1;
            double xlinv = 1.0 / xl1;
            double ulinv = 1.0 / ul1;

            p0[ i ] = 0.0;
            if ( df0dx[ i ] > 0 ) p0[ i ] = df0dx[ i ];
            p0[ i ] += 0.5 * raa0 * ulinv;
            p0[ i ] *= ux2;

            q0[ i ] = 0.0;
            if ( df0dx[ i ] < 0 ) q0[ i ] = -df0dx[ i ];
            q0[ i ] += 0.5 * raa0 * ulinv;
            q0[ i ] *= xl2;

            r0 -= p0[ i ] * uxinv + q0[ i ] * xlinv;

            for ( int j = 0; j < mNumCon; ++j )
            {
                P[ j ][ i ] = 0.0;
                if ( dfdx[ j ][ i ] > 0 ) P[ j ][ i ] = dfdx[ j ][ i ];
                P[ j ][ i ] += 0.5 * ulinv * raa[ j ];
                P[ j ][ i ] *= ux2;

                Q[ j ][ i ] = 0.0;
                if ( dfdx[ j ][ i ] < 0 ) Q[ j ][ i ] = -dfdx[ j ][ i ];
                Q[ j ][ i ] += 0.5 * ulinv * raa[ j ];
                Q[ j ][ i ] *= xl2;

                r[ j ] -= uxinv * P[ j ][ i ] + xlinv * Q[ j ][ i ];
            }
        }

        for ( int j = 0; j < mNumCon; ++j ) b[ j ] = -r[ j ];

        this->subsolve();

        f0app = r0;
        for ( int j = 0; j < mNumCon; ++j ) fapp[ j ] = r[ j ];

        for ( int i = 0; i < mNumVar; ++i )
        {
            double ux1   = upp[ i ] - xmma[ i ];
            double xl1   = xmma[ i ] - low[ i ];
            double uxinv = 1.0 / ux1;
            double xlinv = 1.0 / xl1;

            f0app += uxinv * p0[ i ] + xlinv * q0[ i ];
            for ( int j = 0; j < mNumCon; ++j ) fapp[ j ] += uxinv * P[ j ][ i ] + xlinv * Q[ j ][ i ];
        }
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::raaupdate()
    {
        double raacofmin = 1.0e-7;
        double raacof    = 0.0;

        for ( int i = 0; i < mNumVar; ++i )
        {
            double yyux = ( xmma[ i ] - xval[ i ] ) / ( upp[ i ] - xmma[ i ] );
            double yyxl = ( xmma[ i ] - xval[ i ] ) / ( xmma[ i ] - low[ i ] );
            raacof += 0.5 * yyux * yyxl;
        }
        raacof = mx( raacof, raacofmin );

        if ( f0valnew > ( f0app + 0.5 * epsimin ) )
        {
            double deltaraa0 = ( 1.0 / raacof ) * ( f0valnew - f0app );
            double zz0       = 1.1 * ( raa0 + deltaraa0 );
            zz0              = mn( zz0, 10.0 * raa0 );
            raa0             = zz0;
        }

        for ( int j = 0; j < mNumCon; ++j )
        {
            double fappe  = fapp[ j ] + 0.5 * epsimin;
            double fdelta = fvalnew[ j ] - fappe;
            if ( fdelta > 0.0 )
            {
                double deltaraa = ( 1.0 / raacof ) * ( fvalnew[ j ] - fapp[ j ] );
                double zzz      = 1.1 * ( raa[ j ] + deltaraa );
                raa[ j ]        = mn( zzz, 10.0 * raa[ j ] );
            }
        }
    }

    //------------------------------------------------------------------------------

    int
    Algorithm_MMA::kktcheck( real& aKKTrefres )
    {
        double residunorm = 0.0;
        double residumax  = -1.0;

        for ( int i = 0; i < mNumVar; ++i )
        {
            double rex   = df0dx[ i ] - xsi[ i ] + eta[ i ];
            // complementarity at the committed design xval (== xmma unless the in-place redistance moved it)
            double rexsi = xsi[ i ] * ( xval[ i ] - xmin[ i ] );
            double reeta = eta[ i ] * ( xmax[ i ] - xval[ i ] );

            for ( int j = 0; j < mNumCon; ++j ) rex += dfdx[ j ][ i ] * lam[ j ];

            residunorm += rex * rex;   residumax = mx( residumax, std::abs( rex ) );
            residunorm += rexsi * rexsi; residumax = mx( residumax, std::abs( rexsi ) );
            residunorm += reeta * reeta; residumax = mx( residumax, std::abs( reeta ) );
        }

        double rez   = a0 - zet;
        double rezet = zet * z;

        for ( int j = 0; j < mNumCon; ++j )
        {
            double rey   = c[ j ] + d[ j ] * y[ j ] - mu[ j ] - lam[ j ];
            double relam = fval[ j ] - z * a[ j ] - y[ j ] + s[ j ];
            double remu  = mu[ j ] * y[ j ];
            double res   = lam[ j ] * s[ j ];
            rez -= a[ j ] * lam[ j ];

            residunorm += rey * rey;     residumax = mx( residumax, std::abs( rey ) );
            residunorm += relam * relam; residumax = mx( residumax, std::abs( relam ) );
            residunorm += remu * remu;   residumax = mx( residumax, std::abs( remu ) );
            residunorm += res * res;     residumax = mx( residumax, std::abs( res ) );
        }

        residunorm += rez * rez;     residumax = mx( residumax, std::abs( rez ) );
        residunorm += rezet * rezet; residumax = mx( residumax, std::abs( rezet ) );

        residunorm = std::sqrt( residunorm );

        double KKTres = mx( residunorm, residumax );

        if ( aKKTrefres < 0 ) aKKTrefres = KKTres;

        return ( KKTres > mNormDrop * aKKTrefres ) ? 0 : 1;
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::solve_dense( std::vector< double >& A, std::vector< double >& X,
            const std::vector< double >& B, int n )
    {
        // Gaussian elimination with partial pivoting on row-major A (n x n); solves A x = B.
        // A is overwritten. Matches dgesv (LU partial pivot) numerics for these small systems.
        X = B;
        for ( int k = 0; k < n; ++k )
        {
            // pivot
            int    piv  = k;
            double pmax = std::abs( A[ k * n + k ] );
            for ( int i = k + 1; i < n; ++i )
            {
                double v = std::abs( A[ i * n + k ] );
                if ( v > pmax ) { pmax = v; piv = i; }
            }
            if ( piv != k )
            {
                for ( int j = 0; j < n; ++j ) std::swap( A[ k * n + j ], A[ piv * n + j ] );
                std::swap( X[ k ], X[ piv ] );
            }

            double akk = A[ k * n + k ];
            for ( int i = k + 1; i < n; ++i )
            {
                double f = A[ i * n + k ] / akk;
                for ( int j = k; j < n; ++j ) A[ i * n + j ] -= f * A[ k * n + j ];
                X[ i ] -= f * X[ k ];
            }
        }
        // back substitution
        for ( int i = n - 1; i >= 0; --i )
        {
            double sum = X[ i ];
            for ( int j = i + 1; j < n; ++j ) sum -= A[ i * n + j ] * X[ j ];
            X[ i ] = sum / A[ i * n + i ];
        }
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_MMA::subsolve()
    {
        const int numvar = mNumVar;
        const int numcon = mNumCon;

        for ( int i = 0; i < numvar; ++i )
        {
            x[ i ]   = 0.5 * ( alfa[ i ] + beta[ i ] );
            xsi[ i ] = mx( 1.0, 1.0 / ( x[ i ] - alfa[ i ] ) );
            eta[ i ] = mx( 1.0, 1.0 / ( beta[ i ] - x[ i ] ) );
        }
        for ( int j = 0; j < numcon; ++j )
        {
            y[ j ]   = 1.0;
            lam[ j ] = 1.0;
            mu[ j ]  = mx( 1.0, 0.5 * c[ j ] );
            s[ j ]   = 1.0;
        }
        z   = 1.0;
        zet = 1.0;

        double epsi  = 1.0;
        double dz = 0.0, dzet = 0.0, delz = 0.0;

        std::vector< double > dlocal;   // scratch for dense solve result (size k)

        while ( epsi > epsimin )
        {
            double residunorm = 0.0;
            double residumax  = -1.0;

            for ( int j = 0; j < numcon; ++j ) gvec[ j ] = 0.0;

            for ( int i = 0; i < numvar; ++i )
            {
                double ux1    = upp[ i ] - x[ i ];
                double xl1    = x[ i ] - low[ i ];
                double ux2    = ux1 * ux1;
                double xl2    = xl1 * xl1;
                double uxinv1 = 1.0 / ux1;
                double xlinv1 = 1.0 / xl1;

                double plam = p0[ i ];
                for ( int j = 0; j < numcon; ++j ) plam += P[ j ][ i ] * lam[ j ];
                double qlam = q0[ i ];
                for ( int j = 0; j < numcon; ++j ) qlam += Q[ j ][ i ] * lam[ j ];

                double dpsidx = plam / ux2 - qlam / xl2;
                double rex    = dpsidx - xsi[ i ] + eta[ i ];
                double rexsi  = xsi[ i ] * ( x[ i ] - alfa[ i ] ) - epsi;
                double reeta  = eta[ i ] * ( beta[ i ] - x[ i ] ) - epsi;

                for ( int j = 0; j < numcon; ++j ) gvec[ j ] += P[ j ][ i ] * uxinv1 + Q[ j ][ i ] * xlinv1;

                residunorm += rex * rex;     residumax = mx( residumax, std::abs( rex ) );
                residunorm += rexsi * rexsi; residumax = mx( residumax, std::abs( rexsi ) );
                residunorm += reeta * reeta; residumax = mx( residumax, std::abs( reeta ) );
            }

            double rez   = a0 - zet;
            double rezet = zet * z - epsi;
            for ( int j = 0; j < numcon; ++j )
            {
                double rey   = c[ j ] + d[ j ] * y[ j ] - mu[ j ] - lam[ j ];
                double relam = gvec[ j ] - z * a[ j ] - y[ j ] + s[ j ] - b[ j ];
                double remu  = mu[ j ] * y[ j ] - epsi;
                double res   = lam[ j ] * s[ j ] - epsi;
                rez -= a[ j ] * lam[ j ];

                residunorm += rey * rey;     residumax = mx( residumax, std::abs( rey ) );
                residunorm += relam * relam; residumax = mx( residumax, std::abs( relam ) );
                residunorm += remu * remu;   residumax = mx( residumax, std::abs( remu ) );
                residunorm += res * res;     residumax = mx( residumax, std::abs( res ) );
            }
            residunorm += rez * rez;     residumax = mx( residumax, std::abs( rez ) );
            residunorm += rezet * rezet; residumax = mx( residumax, std::abs( rezet ) );
            residunorm = std::sqrt( residunorm );

            int ittt = 0;
            while ( residumax > 0.9 * epsi && ittt < 200 )
            {
                ittt++;

                for ( int j = 0; j < numcon; ++j ) gvec[ j ] = 0.0;

                for ( int i = 0; i < numvar; ++i )
                {
                    double ux1 = upp[ i ] - x[ i ];
                    double xl1 = x[ i ] - low[ i ];
                    double ux2 = ux1 * ux1;
                    double xl2 = xl1 * xl1;
                    double ux3 = ux1 * ux2;
                    double xl3 = xl1 * xl2;

                    double uxinv1 = 1.0 / ux1;
                    double xlinv1 = 1.0 / xl1;
                    double uxinv2 = 1.0 / ux2;
                    double xlinv2 = 1.0 / xl2;

                    double plam = p0[ i ];
                    for ( int j = 0; j < numcon; ++j ) plam += P[ j ][ i ] * lam[ j ];
                    double qlam = q0[ i ];
                    for ( int j = 0; j < numcon; ++j ) qlam += Q[ j ][ i ] * lam[ j ];

                    double dpsidx = plam / ux2 - qlam / xl2;
                    delx[ i ]     = dpsidx - epsi / ( x[ i ] - alfa[ i ] ) + epsi / ( beta[ i ] - x[ i ] );

                    diagx[ i ]    = plam / ux3 + qlam / xl3;
                    diagx[ i ]    = 2.0 * diagx[ i ] + xsi[ i ] / ( x[ i ] - alfa[ i ] ) + eta[ i ] / ( beta[ i ] - x[ i ] );
                    diagxinv[ i ] = 1.0 / diagx[ i ];

                    for ( int j = 0; j < numcon; ++j )
                    {
                        gvec[ j ]   += P[ j ][ i ] * uxinv1 + Q[ j ][ i ] * xlinv1;
                        GG[ j ][ i ] = P[ j ][ i ] * uxinv2 - Q[ j ][ i ] * xlinv2;
                    }
                }

                delz = a0 - epsi / z;
                for ( int j = 0; j < numcon; ++j )
                {
                    dely[ j ]    = c[ j ] + d[ j ] * y[ j ] - lam[ j ] - epsi / y[ j ];
                    dellam[ j ]  = gvec[ j ] - z * a[ j ] - y[ j ] - b[ j ] + epsi / lam[ j ];
                    diagy[ j ]   = d[ j ] + mu[ j ] / y[ j ];
                    diagyinv[ j ] = 1.0 / diagy[ j ];
                    diaglam[ j ] = s[ j ] / lam[ j ];
                    diaglamyi[ j ] = diaglam[ j ] + diagyinv[ j ];
                    delz -= a[ j ] * lam[ j ];
                }

                if ( numcon < numvar )
                {
                    int    k  = numcon + 1;
                    auto   AT = [ & ]( int i2, int j2 ) -> double& { return AA[ i2 * k + j2 ]; };

                    AT( numcon, numcon ) = -zet / z;

                    for ( int j = 0; j < numcon; ++j )
                    {
                        bb[ j ] = dellam[ j ] + dely[ j ] / diagy[ j ];
                        for ( int i = 0; i < numvar; ++i ) bb[ j ] -= GG[ j ][ i ] * ( delx[ i ] / diagx[ i ] );
                    }
                    bb[ numcon ] = delz;

                    for ( int j = 0; j < numcon; ++j )
                    {
                        for ( int jj = 0; jj < numcon; ++jj )
                        {
                            double sum = 0.0;
                            for ( int i = 0; i < numvar; ++i ) sum += GG[ j ][ i ] * GG[ jj ][ i ] * diagxinv[ i ];
                            AT( j, jj ) = sum;
                        }
                        AT( j, j ) += diaglamyi[ j ];
                        AT( numcon, j ) = a[ j ];
                        AT( j, numcon ) = a[ j ];
                    }

                    solve_dense( AA, dlocal, bb, k );
                    for ( int j = 0; j < numcon; ++j ) dlam[ j ] = dlocal[ j ];
                    dz = dlocal[ numcon ];

                    for ( int i = 0; i < numvar; ++i )
                    {
                        dx[ i ] = -delx[ i ] / diagx[ i ];
                        for ( int j = 0; j < numcon; ++j ) dx[ i ] -= ( GG[ j ][ i ] * dlam[ j ] ) / diagx[ i ];
                    }
                }
                else
                {
                    int    k  = numvar + 1;
                    auto   AT = [ & ]( int i2, int j2 ) -> double& { return AA[ i2 * k + j2 ]; };

                    double azz = zet / z;
                    for ( int j = 0; j < numcon; ++j )
                    {
                        diaglamyiinv[ j ] = 1.0 / diaglamyi[ j ];
                        dellamyi[ j ]     = dellam[ j ] + dely[ j ] / diagy[ j ];
                        azz += a[ j ] * a[ j ] / diaglamyi[ j ];
                    }
                    AT( numvar, numvar ) = azz;

                    for ( int i = 0; i < numvar; ++i )
                    {
                        for ( int ii = 0; ii < numvar; ++ii )
                        {
                            double sum = 0.0;
                            for ( int j = 0; j < numcon; ++j ) sum += GG[ j ][ i ] * GG[ j ][ ii ] * diaglamyiinv[ j ];
                            AT( i, ii ) = sum;
                        }
                        AT( i, i ) += diagx[ i ];
                        double col = 0.0;
                        for ( int j = 0; j < numcon; ++j ) col -= GG[ j ][ i ] * a[ j ] / diaglamyi[ j ];
                        AT( numvar, i ) = col;
                        AT( i, numvar ) = col;
                        bb[ i ]         = -delx[ i ];
                        for ( int j = 0; j < numcon; ++j ) bb[ i ] -= GG[ j ][ i ] * dellamyi[ j ] / diaglamyi[ j ];
                    }
                    bb[ numvar ] = -delz;
                    for ( int j = 0; j < numcon; ++j ) bb[ numvar ] += a[ j ] * dellamyi[ j ] / diaglamyi[ j ];

                    solve_dense( AA, dlocal, bb, k );
                    for ( int i = 0; i < numvar; ++i ) dx[ i ] = dlocal[ i ];
                    dz = dlocal[ numvar ];

                    for ( int j = 0; j < numcon; ++j )
                    {
                        dlam[ j ] = -dz * a[ j ] / diaglamyi[ j ] + dellamyi[ j ] / diaglamyi[ j ];
                        for ( int i = 0; i < numvar; ++i ) dlam[ j ] += GG[ j ][ i ] * dx[ i ] / diaglamyi[ j ];
                    }
                }

                dzet = -zet + epsi / z - zet * dz / z;

                double stmxx   = -1.0e99;
                double stmalfa = -1.0e99;
                double stmbeta = -1.0e99;

                stmxx = mx( stmxx, -1.01 * dz / z );
                stmxx = mx( stmxx, -1.01 * dzet / zet );

                for ( int i = 0; i < numvar; ++i )
                {
                    dxsi[ i ] = -xsi[ i ] + epsi / ( x[ i ] - alfa[ i ] ) - ( xsi[ i ] * dx[ i ] ) / ( x[ i ] - alfa[ i ] );
                    deta[ i ] = -eta[ i ] + epsi / ( beta[ i ] - x[ i ] ) + ( eta[ i ] * dx[ i ] ) / ( beta[ i ] - x[ i ] );

                    stmxx   = mx( stmxx, -1.01 * dxsi[ i ] / xsi[ i ] );
                    stmxx   = mx( stmxx, -1.01 * deta[ i ] / eta[ i ] );
                    stmalfa = mx( stmalfa, -1.01 * dx[ i ] / ( x[ i ] - alfa[ i ] ) );
                    stmbeta = mx( stmbeta, 1.01 * dx[ i ] / ( beta[ i ] - x[ i ] ) );
                }

                for ( int j = 0; j < numcon; ++j )
                {
                    dy[ j ]  = -dely[ j ] / diagy[ j ] + dlam[ j ] / diagy[ j ];
                    dmu[ j ] = -mu[ j ] + epsi / y[ j ] - ( mu[ j ] * dy[ j ] ) / y[ j ];
                    ds[ j ]  = -s[ j ] + epsi / lam[ j ] - ( s[ j ] * dlam[ j ] ) / lam[ j ];

                    stmxx = mx( stmxx, -1.01 * dy[ j ] / y[ j ] );
                    stmxx = mx( stmxx, -1.01 * dlam[ j ] / lam[ j ] );
                    stmxx = mx( stmxx, -1.01 * dmu[ j ] / mu[ j ] );
                    stmxx = mx( stmxx, -1.01 * ds[ j ] / s[ j ] );
                }

                double stmalbe   = mx( stmalfa, stmbeta );
                double stmalbexx = mx( stmalbe, stmxx );
                double stminv    = mx( stmalbexx, 1.0 );
                double steg      = 1.0 / stminv;

                double zold   = z;
                double zetold = zet;
                for ( int i = 0; i < numvar; ++i ) { xold[ i ] = x[ i ]; xsiold[ i ] = xsi[ i ]; etaold[ i ] = eta[ i ]; }
                for ( int j = 0; j < numcon; ++j ) { yold[ j ] = y[ j ]; lamold[ j ] = lam[ j ]; muold[ j ] = mu[ j ]; sold[ j ] = s[ j ]; }

                int    itto    = 0;
                double resinew = 2.0 * residunorm;
                double resimax = residumax;

                while ( resinew > residunorm && itto < 50 )
                {
                    itto++;
                    resinew = 0.0;
                    resimax = -1.0;

                    z   = zold + steg * dz;
                    zet = zetold + steg * dzet;
                    for ( int i = 0; i < numvar; ++i )
                    {
                        x[ i ]   = xold[ i ] + steg * dx[ i ];
                        xsi[ i ] = xsiold[ i ] + steg * dxsi[ i ];
                        eta[ i ] = etaold[ i ] + steg * deta[ i ];
                    }
                    for ( int j = 0; j < numcon; ++j )
                    {
                        y[ j ]   = yold[ j ] + steg * dy[ j ];
                        lam[ j ] = lamold[ j ] + steg * dlam[ j ];
                        mu[ j ]  = muold[ j ] + steg * dmu[ j ];
                        s[ j ]   = sold[ j ] + steg * ds[ j ];
                        gvec[ j ] = 0.0;
                    }

                    for ( int i = 0; i < numvar; ++i )
                    {
                        double ux1    = upp[ i ] - x[ i ];
                        double xl1    = x[ i ] - low[ i ];
                        double ux2    = ux1 * ux1;
                        double xl2    = xl1 * xl1;
                        double uxinv1 = 1.0 / ux1;
                        double xlinv1 = 1.0 / xl1;

                        double plam = p0[ i ];
                        for ( int j = 0; j < numcon; ++j ) plam += P[ j ][ i ] * lam[ j ];
                        double qlam = q0[ i ];
                        for ( int j = 0; j < numcon; ++j ) qlam += Q[ j ][ i ] * lam[ j ];

                        double dpsidx = plam / ux2 - qlam / xl2;
                        double rex    = dpsidx - xsi[ i ] + eta[ i ];
                        double rexsi  = xsi[ i ] * ( x[ i ] - alfa[ i ] ) - epsi;
                        double reeta  = eta[ i ] * ( beta[ i ] - x[ i ] ) - epsi;

                        for ( int j = 0; j < numcon; ++j ) gvec[ j ] += P[ j ][ i ] * uxinv1 + Q[ j ][ i ] * xlinv1;

                        resinew += rex * rex;     resimax = mx( resimax, std::abs( rex ) );
                        resinew += rexsi * rexsi; resimax = mx( resimax, std::abs( rexsi ) );
                        resinew += reeta * reeta; resimax = mx( resimax, std::abs( reeta ) );
                    }

                    rez   = a0 - zet;
                    rezet = zet * z - epsi;
                    for ( int j = 0; j < numcon; ++j )
                    {
                        double rey   = c[ j ] + d[ j ] * y[ j ] - mu[ j ] - lam[ j ];
                        double relam = gvec[ j ] - z * a[ j ] - y[ j ] + s[ j ] - b[ j ];
                        double remu  = mu[ j ] * y[ j ] - epsi;
                        double res   = lam[ j ] * s[ j ] - epsi;
                        rez -= a[ j ] * lam[ j ];

                        resinew += rey * rey;     resimax = mx( resimax, std::abs( rey ) );
                        resinew += relam * relam; resimax = mx( resimax, std::abs( relam ) );
                        resinew += remu * remu;   resimax = mx( resimax, std::abs( remu ) );
                        resinew += res * res;     resimax = mx( resimax, std::abs( res ) );
                    }
                    resinew += rez * rez;     resimax = mx( resimax, std::abs( rez ) );
                    resinew += rezet * rezet; resimax = mx( resimax, std::abs( rezet ) );
                    resinew = std::sqrt( resinew );
                    steg    = steg / 2.0;
                }

                residunorm = resinew;
                residumax  = resimax;
                steg *= 2.0;
            }
            epsi *= 0.1;
        }

        zmma   = z;
        zetmma = zet;
        for ( int i = 0; i < numvar; ++i ) xmma[ i ] = x[ i ];
    }

    //------------------------------------------------------------------------------

    int
    Algorithm_MMA::mma_solve()
    {
        this->mma_initialize();

        // fresh solve (possibly after a restart): objective scale re-derives from the first grad()
        mObjScale = 1.0;

        int iter  = 0;
        int istop = 0;
        int nfunc = 0;
        int nsens = 0;

        // enforce bounds and validate
        int error1 = 0, error2 = 0;
        for ( int i = 0; i < mNumVar; ++i )
        {
            xval[ i ] = mx( mn( xval[ i ], xmax[ i ] ), xmin[ i ] );
            if ( std::abs( xmax[ i ] - xmin[ i ] ) < 1e-18 ) error1 = 1;
            if ( xmax[ i ] <= xmin[ i ] ) error2 = 1;
        }
        MORIS_ERROR( error1 == 0, "Native MMA: upper and lower bounds of some variables are numerically identical." );
        MORIS_ERROR( error2 == 0, "Native MMA: upper bound is smaller or equal to lower bound for some variables." );

        // the suppression flag lives in the shared interface: never leave it set when this solve returns (restart, end)
        struct Reinit_Suppression_Reset
        {
            std::shared_ptr< Problem > mProb;
            ~Reinit_Suppression_Reset() { mProb->set_reinit_suppressed( false ); }
        } tSuppressionReset{ mProblem };

        // initial evaluation (an accepted iterate: redistanced in every in-place mode)
        mProblem->set_reinit_suppressed( false );
        this->func( xval.data(), f0val, fval );

        // The workflow can request an optimization restart during a criteria evaluation
        // (adaptive remeshing or a scheduled reinitialization invalidates the current ADV
        // vector). Mirror the TPL GCMMA wrapper: exit the solve loop immediately so the
        // Manager re-initializes the problem (fresh mesh + ADVs) and re-enters solve.
        // Backtracking instead would iterate on a stale design AND on frozen sensitivities
        // (the pending-restart flag makes the criteria interface return cached gradients).
        if ( mProblem->restart_optimization() )
        {
            MORIS_LOG_INFO( "Native MMA: optimization restart requested during initial evaluation; exiting solve loop." );
            return istop;
        }

        // in-place redistance mode: the initial evaluation may have redistanced the design; start from the
        // redistanced point so the gradient call sees the ADVs the criteria were computed for (MB-013: the
        // consistency assert in Problem::compute_design_criteria_gradients fired here)
        if ( mReinitHook && (int)mLastRedistanced.size() == mNumVar )
        {
            for ( int i = 0; i < mNumVar; ++i ) xval[ i ] = mLastRedistanced[ i ];
        }

        this->grad( xval.data() );
        nfunc++; nsens++;

        // initialize raa0 and raa
        if ( mMaxInnerIter )
        {
            raa0 = 1.0;
            for ( int j = 0; j < mNumCon; ++j ) raa[ j ] = 1.0;
        }
        else
        {
            raa0 = 1.0e-5;
            for ( int j = 0; j < mNumCon; ++j ) raa[ j ] = 1.0e-5;
        }
        double sclraa = 0.1 / mNumVar;
        for ( int i = 0; i < mNumVar; ++i )
        {
            raa0 += sclraa * std::abs( df0dx[ i ] ) / ( xmax[ i ] - xmin[ i ] );
            for ( int j = 0; j < mNumCon; ++j ) raa[ j ] += sclraa * std::abs( dfdx[ j ][ i ] ) / ( xmax[ i ] - xmin[ i ] );
        }

        // constraint penalty c_j
        double penFac = ( mPenalty < 0 ) ? std::abs( f0val ) : mPenalty;
        for ( int j = 0; j < mNumCon; ++j ) c[ j ] = penFac;

        real tKKTrefres = -1.0;

        while ( iter < mMaxIterations && istop == 0 )
        {
            iter++;

            mMoveScale = 1.0;
            this->asymp( iter );
            this->gcmmasub( iter );

            const int tNumFuncAtStart = nfunc;
            this->trial_func( iter, 0, nfunc );

            // restart takes precedence over the NaN backtrack: the workflow has invalidated
            // the current design/mesh and the Manager must re-initialize the problem
            if ( mProblem->restart_optimization() )
            {
                MORIS_LOG_INFO( "Native MMA iter %d: optimization restart requested; exiting solve loop.", (int)iter );
                return istop;
            }

            // NaN guard on the first trial (see backtrack_non_finite)
            Trial_Status tTrialStatus = this->backtrack_non_finite( iter, nfunc );
            if ( tTrialStatus == Trial_Status::restart )
            {
                return istop;
            }

            bool tConservative = tTrialStatus == Trial_Status::finite && this->is_conservative();

            int innerit = 0;
            while ( tTrialStatus == Trial_Status::finite && !tConservative && innerit < mMaxInnerIter )
            {
                innerit++;
                this->raaupdate();
                this->gcmmasub( iter );
                this->trial_func( iter, innerit, nfunc );

                if ( mProblem->restart_optimization() )
                {
                    MORIS_LOG_INFO( "Native MMA iter %d: optimization restart requested during inner iteration; exiting solve loop.", (int)iter );
                    return istop;
                }

                // a failed inner trial must not reach the conservatism test: NaN comparisons are false, so it
                // would pass as conservative and be committed and differentiated (MB-023)
                tTrialStatus = this->backtrack_non_finite( iter, nfunc );
                if ( tTrialStatus == Trial_Status::restart )
                {
                    return istop;
                }
                if ( tTrialStatus == Trial_Status::failed )
                {
                    break;
                }

                tConservative = this->is_conservative();
            }

            if ( tTrialStatus == Trial_Status::failed )
            {
                // could not recover: hold the previous (finite) design. The failed evaluations overwrote the
                // shared Problem state, so re-evaluate xval before continuing; the Manager reads that state (MB-024)
                if ( !this->restore_accepted_state( iter, nfunc ) )
                {
                    return istop;    // restart requested
                }
                MORIS_LOG_INFO( "Native MMA iter %d: holding previous design and continuing.", (int)iter );
                continue;
            }

            // one summary line per outer iteration (inner-loop cap statistics)
            MORIS_LOG_INFO( "Native MMA iter %d: accepted after %d inner cycles (cap %d)%s, %d trial evaluations",
                    (int)iter, innerit, (int)mMaxInnerIter, tConservative ? ", conservative" : ", NOT conservative (cap reached)",
                    nfunc - tNumFuncAtStart );

            // SDF taper: tell the workflow about the accepted iterate (its step from the committed design, before any
            // redistance) so that it can wean the in-place redistance off
            if ( mReinitHook )
            {
                Vector< real > tPrevious( mNumVar );
                Vector< real > tAccepted( mNumVar );
                for ( int i = 0; i < mNumVar; ++i )
                {
                    tPrevious( i ) = xval[ i ];
                    tAccepted( i ) = xmma[ i ];
                }
                mProblem->notify_accepted_iterate( tPrevious, tAccepted );
            }

            // accept-only redistance: the trials were evaluated at the raw design; redistance the accepted xmma once and
            // re-evaluate there, so that the committed design, its criteria and (below) its gradient all belong to
            // phi = R(xmma). The trial criteria are no longer needed: f0valnew/fvalnew are overwritten with those at phi.
            // (skipped once the redistance is tapered off: the trial evaluation at xmma is then the committed one)
            if ( this->accept_only_redistance() && mProblem->is_reinit_active() )
            {
                const real            tF0Psi = f0valnew;
                std::vector< double > tFPsi  = fvalnew;

                mProblem->set_reinit_suppressed( false );
                this->func( xmma.data(), f0valnew, fvalnew );
                nfunc++;

                if ( mProblem->restart_optimization() )
                {
                    MORIS_LOG_INFO( "Native MMA iter %d: optimization restart requested during the accept redistance; exiting solve loop.", (int)iter );
                    return istop;
                }

                std::string tLine = "MMA-ACCEPT iter=" + std::to_string( iter );
                char        tBuf[ 256 ];
                std::snprintf( tBuf, sizeof( tBuf ), " obj psi=%.12e phi=%.12e", tF0Psi, f0valnew );
                tLine += tBuf;
                for ( int j = 0; j < mNumCon; ++j )
                {
                    std::snprintf( tBuf, sizeof( tBuf ), " | c%d psi=%.12e phi=%.12e", j, tFPsi[ j ], fvalnew[ j ] );
                    tLine += tBuf;
                }
                MORIS_LOG_INFO( "%s", tLine.c_str() );

                if ( !this->criteria_finite() )
                {
                    // the redistanced design does not evaluate: commit the accepted raw design instead (delta = 0)
                    MORIS_LOG_WARNING( "Native MMA iter %d: criteria at R(psi) are non-finite; committing the raw accepted design.", (int)iter );
                    mProblem->set_reinit_suppressed( true );
                    this->func( xmma.data(), f0valnew, fvalnew );
                    nfunc++;
                    if ( mProblem->restart_optimization() )
                    {
                        return istop;
                    }
                }
            }

            // advance design. In reinit mode commit the redistanced design (mLastRedistanced holds
            // the redistance of the accepted xmma); otherwise commit the raw subproblem solution.
            const bool tUseRedist = mReinitHook && (int)mLastRedistanced.size() == mNumVar;
            for ( int i = 0; i < mNumVar; ++i )
            {
                xold2[ i ] = xold1[ i ];
                xold1[ i ] = xval[ i ];
                xval[ i ]  = tUseRedist ? mLastRedistanced[ i ] : xmma[ i ];
            }

            // Warm restart across the redistance: express the optimizer's history in the redistanced
            // coordinates. The redistance R leaves the zero contour (hence the criteria) unchanged, so it
            // is not a step: shift xold1, xold2 and the asymptotes by delta = R(xmma) - xmma, AFTER the
            // rotation. asymp() then sees x_k - x_{k-1} = xmma - phi_{k-1} (the subproblem's own step) and
            // keeps xold1 - low / upp - xold1 exactly. Without the shift, far-field redistance jumps
            // (0.1-0.15 vs a 0.015 move limit) drive the oscillation test and collapse the asymptotes
            // onto the floor (measured: 641/4961 variables locked after 26 iterations).
            // Design note: runs/studies/lunch_and_learn_xiga (2026-10-10), MB-013.
            if ( tUseRedist && mShiftHistory )
            {
                for ( int i = 0; i < mNumVar; ++i )
                {
                    const double tDelta = xval[ i ] - xmma[ i ];
                    xold1[ i ] += tDelta;
                    xold2[ i ] += tDelta;
                    low[ i ] += tDelta;
                    upp[ i ] += tDelta;
                    MORIS_ASSERT( low[ i ] < xval[ i ] && xval[ i ] < upp[ i ],
                            "Native MMA: asymptotes on the wrong side of the redistanced design (variable %d).", i );
                }
            }
            f0val = f0valnew;
            for ( int j = 0; j < mNumCon; ++j ) fval[ j ] = fvalnew[ j ];

            {
                // committed design of this outer iteration (trajectory record; before grad(), which may re-book the scale)
                std::string tLine = "MMA-COMMIT iter=" + std::to_string( iter );
                char        tBuf[ 128 ];
                std::snprintf( tBuf, sizeof( tBuf ), " obj=%.12e", f0val / mObjScale );
                tLine += tBuf;
                for ( int j = 0; j < mNumCon; ++j )
                {
                    std::snprintf( tBuf, sizeof( tBuf ), " c%d=%.12e", j, fval[ j ] );
                    tLine += tBuf;
                }
                MORIS_LOG_INFO( "%s", tLine.c_str() );
            }

            this->grad( xval.data() );
            nsens++;

            istop = this->kktcheck( tKKTrefres );
        }

        return ( iter < mMaxIterations ) ? 0 : 1;
    }
}    // namespace moris::opt
