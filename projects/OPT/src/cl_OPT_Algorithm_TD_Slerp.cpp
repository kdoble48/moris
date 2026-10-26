/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_TD_Slerp.cpp
 *
 * Topological-derivative level-set evolution on the unit sphere (see header). MORIS owns the loop, as in Algorithm_MMA.
 */

#include <algorithm>
#include <cmath>

#include "cl_OPT_Algorithm_TD_Slerp.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

namespace moris::opt
{
    //------------------------------------------------------------------------------

    namespace
    {
        // degrees <-> radians, the factors of Python's math.degrees / math.radians (x / f and x * f), so the angle
        // round trip matches the offline driver this algorithm reproduces
        const real gDegToRad = M_PI / 180.0;
    }    // namespace

    //------------------------------------------------------------------------------

    std::string
    Algorithm_TD_Slerp::get_exit_name( Exit aExit )
    {
        switch ( aExit )
        {
            case Exit::NONE:
                return "NONE";
            case Exit::THETA:
                return "THETA";
            case Exit::THETA_PI:
                return "THETA_PI";
            case Exit::N_MAX:
                return "N_MAX";
            case Exit::LINE_SEARCH:
                return "LINE_SEARCH";
            case Exit::G_ZERO:
                return "G_ZERO";
        }
        return "UNKNOWN";
    }

    //------------------------------------------------------------------------------

    Algorithm_TD_Slerp::Algorithm_TD_Slerp( const Parameter_List& aParameterList )
            : mThetaTol( aParameterList.get< real >( "theta_tol" ) )
            , mThetaPiTol( aParameterList.get< real >( "theta_pi_tol" ) )
            , mKappaInit( aParameterList.get< real >( "kappa_init" ) )
            , mKappaMin( aParameterList.get< real >( "kappa_min" ) )
            , mAcceptRelTol( aParameterList.get< real >( "accept_rel_tol" ) )
            , mNormalizeStart( aParameterList.get< bool >( "normalize_start" ) )
    {
        mRestartIndex         = aParameterList.get< sint >( "restart_index" );
        mMaxIterationsInitial = aParameterList.get< sint >( "max_its" );
        mMaxIterations        = mMaxIterationsInitial;

        const std::string tPolicy = aParameterList.get< std::string >( "kappa_policy" );
        MORIS_ERROR( tPolicy == "reset" || tPolicy == "carry",
                "Algorithm_TD_Slerp: kappa_policy must be \"reset\" or \"carry\" (got \"%s\").",
                tPolicy.c_str() );
        mKappaCarry = tPolicy == "carry";

        MORIS_ERROR( mKappaInit > 0.0 && mKappaInit <= 1.0, "Algorithm_TD_Slerp: kappa_init must be in (0, 1]." );
        MORIS_ERROR( mKappaMin > 0.0 && mKappaMin <= mKappaInit, "Algorithm_TD_Slerp: kappa_min must be in (0, kappa_init]." );
        MORIS_ERROR( mThetaTol >= 0.0 && mThetaPiTol >= 0.0 && mThetaTol + mThetaPiTol < 180.0,
                "Algorithm_TD_Slerp: theta_tol and theta_pi_tol must be >= 0 with theta_tol + theta_pi_tol < 180 degrees." );
        MORIS_ERROR( mAcceptRelTol >= 0.0, "Algorithm_TD_Slerp: accept_rel_tol must be >= 0." );
        MORIS_ERROR( mMaxIterationsInitial >= 0 && mRestartIndex >= 0, "Algorithm_TD_Slerp: max_its and restart_index must be >= 0." );

        MORIS_LOG_INFO(
                "TD_SLERP initialized: max_its=%d, restart_index=%d, theta_tol=%.4g deg, theta_pi_tol=%.4g deg, "
                "kappa_init=%.6g, kappa_min=%.6g, accept_rel_tol=%.3e (J_t <= J_i (1 + tol)), kappa_policy=%s, "
                "normalize_start=%s",
                (int)mMaxIterations,
                (int)mRestartIndex,
                mThetaTol,
                mThetaPiTol,
                mKappaInit,
                mKappaMin,
                mAcceptRelTol,
                tPolicy.c_str(),
                mNormalizeStart ? "true" : "false" );
    }

    //------------------------------------------------------------------------------

    Algorithm_TD_Slerp::~Algorithm_TD_Slerp() {}

    //------------------------------------------------------------------------------

    uint
    Algorithm_TD_Slerp::solve( uint aCurrentOptAlgInd, std::shared_ptr< Problem > aOptProb )
    {
        Tracer tTracer( "OptimizationAlgorithm", "TDSlerp", "Solve" );

        mRunning          = opt::Task::wait;
        mCurrentOptAlgInd = aCurrentOptAlgInd;
        mProblem          = aOptProb;

        if ( par_rank() == 0 )
        {
            this->td_slerp_solve();

            mRunning = opt::Task::exit;
            this->communicate_running_status();
        }
        else
        {
            this->dummy_solve();
        }

        uint tOptIter = gLogger.get_opt_iteration();
        gLogger.set_iteration( "OPT", "Manager", "Perform", tOptIter );

        return tOptIter;
    }

    //------------------------------------------------------------------------------

    real
    Algorithm_TD_Slerp::mass_inner( const Vector< real >& aFirst, const Vector< real >& aSecond )
    {
        return mProblem->mass_inner_product( aFirst, aSecond );
    }

    //------------------------------------------------------------------------------

    real
    Algorithm_TD_Slerp::mass_norm( const Vector< real >& aVector )
    {
        real tSquare = this->mass_inner( aVector, aVector );
        MORIS_ERROR( std::isfinite( tSquare ) && tSquare >= 0.0,
                "Algorithm_TD_Slerp: non-finite or negative squared M norm (%e); is M symmetric positive definite?",
                tSquare );
        return std::sqrt( tSquare );
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_TD_Slerp::check_bounds( const Vector< real >& aADVs )
    {
        const Vector< real >& tLower = mProblem->get_lower_bounds();
        const Vector< real >& tUpper = mProblem->get_upper_bounds();
        for ( uint iADV = 0; iADV < aADVs.size(); iADV++ )
        {
            MORIS_ERROR( std::isfinite( aADVs( iADV ) ), "Algorithm_TD_Slerp: non-finite design coefficient %u.", iADV );
            MORIS_ERROR( aADVs( iADV ) > tLower( iADV ) && aADVs( iADV ) < tUpper( iADV ),
                    "Algorithm_TD_Slerp: design coefficient %u = %e at or beyond its bounds [%e, %e]; the update has no bound "
                    "projection (it would leave the sphere): raise discretization bounds.",
                    iADV,
                    aADVs( iADV ),
                    tLower( iADV ),
                    tUpper( iADV ) );
        }
    }

    //------------------------------------------------------------------------------

    real
    Algorithm_TD_Slerp::evaluate( Vector< real >& aADVs )
    {
        this->check_bounds( aADVs );

        // restart file of this design before the evaluation (MB-016); the evaluation advances the optimization
        // iteration, so name the file from the index captured before it (MB-021)
        const uint tOptIter = gLogger.get_opt_iteration();
        this->write_advs_to_file( aADVs, tOptIter );

        const Vector< real > tSubmitted = aADVs;
        this->compute_design_criteria( aADVs );
        mNumEvaluations++;

        // the workflow must not edit the design (in-place reinitialization): the iterate would leave the sphere
        real tChange = 0.0;
        for ( uint iADV = 0; iADV < aADVs.size(); iADV++ )
        {
            tChange = std::max( tChange, std::abs( aADVs( iADV ) - tSubmitted( iADV ) ) );
        }
        MORIS_ERROR( tChange == 0.0,
                "Algorithm_TD_Slerp: the criteria evaluation changed the design (max change %e); switch off level-set "
                "reinitialization (reinitialize_interface_iter, sdf_*) for TD_SLERP.",
                tChange );

        const real tObjective = this->get_objectives()( 0 );
        MORIS_ERROR( std::isfinite( tObjective ), "Algorithm_TD_Slerp: non-finite objective at evaluation %u.", mNumEvaluations );
        return tObjective;
    }

    //------------------------------------------------------------------------------

    void
    Algorithm_TD_Slerp::td_slerp_solve()
    {
        const uint tNumADVs = mProblem->get_num_advs();
        MORIS_ERROR( tNumADVs > 0, "Algorithm_TD_Slerp: no design variables." );

        mExit           = Exit::NONE;
        mNumEvaluations = 0;
        mObjectiveHistory.clear();
        mThetaHistory.clear();
        mThetaEuclideanHistory.clear();
        mKappaHistory.clear();

        // start: the problem's ADVs (the restart ADVs when a restart file is set), evaluated as given: M is the
        // interface's design mass matrix, available only once an evaluation has built it
        Vector< real > tPhi  = mProblem->get_advs();
        uint           tIter = (uint)mRestartIndex;    // accepted iterations
        real           tJ    = this->evaluate( tPhi );

        // on the unit sphere: normalize and evaluate again unless |phi_0|_M is 1 already (a restart)
        const real tNorm0 = this->mass_norm( tPhi );
        MORIS_ERROR( tNorm0 > 0.0, "Algorithm_TD_Slerp: the start design has zero M norm." );
        if ( mNormalizeStart && std::abs( tNorm0 - 1.0 ) > 1.0e-12 )
        {
            for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
            {
                tPhi( iADV ) = tPhi( iADV ) / tNorm0;
            }
            tJ = this->evaluate( tPhi );
        }
        MORIS_LOG_INFO( "TD_SLERP start: |phi_0|_M = %.15e as given, %.15e used (J = %.15e)", tNorm0, this->mass_norm( tPhi ), tJ );

        // objectives and constraints of the accepted design, restored in the problem after a failed line search
        Matrix< DDRMat > tAcceptedObjectives  = mProblem->get_objectives();
        Matrix< DDRMat > tAcceptedConstraints = mProblem->get_constraints();

        // g~ only right after an accepted evaluation (the interface holds the last evaluation)
        Vector< real > tG = mProblem->get_design_direction();
        MORIS_ERROR( tG.size() == tNumADVs, "Algorithm_TD_Slerp: design direction has %zu entries, %u ADVs.", tG.size(), tNumADVs );

        real tLastKappa = mKappaInit;

        while ( true )
        {
            const real tNormPhi = this->mass_norm( tPhi );
            const real tNormG   = this->mass_norm( tG );

            real tMaxAbs = 0.0;
            for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
            {
                tMaxAbs = std::max( tMaxAbs, std::abs( tPhi( iADV ) ) );
            }

            if ( tNormG == 0.0 )
            {
                MORIS_LOG_INFO( "TD_SLERP iter %u: J = %.15e, |g~|_M = 0", tIter, tJ );
                mObjectiveHistory.push_back( tJ );
                mExit = Exit::G_ZERO;
                break;
            }

            // angle (degrees), M inner product; Euclidean as a diagnostic
            const real tCos   = this->mass_inner( tPhi, tG ) / ( tNormPhi * tNormG );
            const real tTheta = std::acos( std::max( -1.0, std::min( 1.0, tCos ) ) ) / gDegToRad;

            real tDotE = 0.0, tPhiE = 0.0, tGE = 0.0;
            for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
            {
                tDotE += tPhi( iADV ) * tG( iADV );
                tPhiE += tPhi( iADV ) * tPhi( iADV );
                tGE += tG( iADV ) * tG( iADV );
            }
            const real tCosE   = tDotE / ( std::sqrt( tPhiE ) * std::sqrt( tGE ) );
            const real tThetaE = std::acos( std::max( -1.0, std::min( 1.0, tCosE ) ) ) / gDegToRad;

            mObjectiveHistory.push_back( tJ );
            mThetaHistory.push_back( tTheta );
            mThetaEuclideanHistory.push_back( tThetaE );

            MORIS_LOG_INFO( "TD_SLERP iter %u: J = %.15e theta = %.12f deg (euclidean %.12f) |phi|_M = %.15e |g~|_M = %.15e max|c| = %.6e",
                    tIter,
                    tJ,
                    tTheta,
                    tThetaE,
                    tNormPhi,
                    tNormG,
                    tMaxAbs );

            if ( tTheta < mThetaTol )
            {
                mExit = Exit::THETA;
                break;
            }
            if ( tTheta > 180.0 - mThetaPiTol )
            {
                if ( tIter == (uint)mRestartIndex )
                {
                    MORIS_LOG_WARNING(
                            "TD_SLERP: theta = %.4f deg at the first iteration (g~ ~ -s phi): check the sign "
                            "convention (material phase on the negative side, the TD's sign).",
                            tTheta );
                }
                mExit = Exit::THETA_PI;
                break;
            }
            if ( tIter >= (uint)mMaxIterations )
            {
                mExit = Exit::N_MAX;
                break;
            }

            // line search on the great circle through phi and g~: kappa halved until J does not rise
            const real tThetaRad = tTheta * gDegToRad;
            const real tSinTheta = std::sin( tThetaRad );

            real           tKappa    = mKappaCarry ? tLastKappa : mKappaInit;
            bool           tAccepted = false;
            real           tJTrial   = 0.0;
            Vector< real > tTrial( tNumADVs );

            while ( true )
            {
                const real tSinPhi = std::sin( ( 1.0 - tKappa ) * tThetaRad );
                const real tSinG   = std::sin( tKappa * tThetaRad );
                for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
                {
                    tTrial( iADV ) = ( tSinPhi * tPhi( iADV ) / tNormPhi + tSinG * tG( iADV ) / tNormG ) / tSinTheta;
                }

                // the one admissible rescaling: back onto the unit sphere
                const real tNormTrial = this->mass_norm( tTrial );
                for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
                {
                    tTrial( iADV ) = tTrial( iADV ) / tNormTrial;
                }

                tJTrial = this->evaluate( tTrial );

                tAccepted = tJTrial <= tJ * ( 1.0 + mAcceptRelTol );

                MORIS_LOG_INFO( "TD_SLERP iter %u trial: kappa = %.10g J = %.15e (J_i = %.15e) %s",
                        tIter + 1,
                        tKappa,
                        tJTrial,
                        tJ,
                        tAccepted ? "accepted" : "rejected" );

                if ( tAccepted )
                {
                    break;
                }

                tKappa *= 0.5;
                if ( tKappa < mKappaMin )
                {
                    break;
                }
            }

            if ( !tAccepted )
            {
                // restore phi_i: the last restart file holds the accepted design, the problem its objective and
                // constraints (its ADVs and the workflow's state stay those of the last trial; no re-evaluation)
                this->write_advs_to_file( tPhi );
                mProblem->set_objectives_and_constraints( tAcceptedObjectives, tAcceptedConstraints );
                MORIS_LOG_INFO( "TD_SLERP iter %u: line search failed below kappa_min = %.6g; restored phi_%u (J = %.15e)",
                        tIter + 1,
                        mKappaMin,
                        tIter,
                        tJ );
                mExit = Exit::LINE_SEARCH;
                break;
            }

            tPhi                 = tTrial;
            tJ                   = tJTrial;
            tAcceptedObjectives  = mProblem->get_objectives();
            tAcceptedConstraints = mProblem->get_constraints();
            tLastKappa           = tKappa;
            mKappaHistory.push_back( tKappa );
            tIter++;

            tG = mProblem->get_design_direction();
            MORIS_ERROR( tG.size() == tNumADVs, "Algorithm_TD_Slerp: design direction has %zu entries, %u ADVs.", tG.size(), tNumADVs );
        }

        mNumAccepted = tIter;

        MORIS_LOG_INFO( "TD_SLERP exit %s: accepted iterations %u (restart_index %d), evaluations %u, J = %.15e, theta = %.12f deg",
                get_exit_name( mExit ).c_str(),
                tIter,
                (int)mRestartIndex,
                mNumEvaluations,
                tJ,
                mThetaHistory.empty() ? 0.0 : mThetaHistory.back() );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::opt
