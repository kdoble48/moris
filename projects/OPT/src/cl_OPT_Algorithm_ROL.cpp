/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_ROL.cpp
 *
 */

#include "cl_OPT_Algorithm_ROL.hpp"
#include "cl_Communication_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include <iostream>

// Third party header files
#ifdef MORIS_HAVE_ROL
#include "ROL_Ptr.hpp"
#include "ROL_Types.hpp"
#include "ROL_ParameterList.hpp"
#include "ROL_StdVector.hpp"
#include "ROL_Bounds.hpp"
#include "ROL_Objective.hpp"
#include "ROL_Constraint.hpp"
#include "ROL_Problem.hpp"
#include "ROL_Solver.hpp"
#include "ROL_SecantFactory.hpp"
#include "ROL_TypeB_AlgorithmFactory.hpp"
#endif

#include <algorithm>
#include <cmath>

using namespace moris;

namespace moris::opt
{
    namespace
    {
        // Thrown by ensure_criteria when max_forward_solves is reached; caught in rol_solve. Not derived from
        // std::exception on purpose, so ROL's internal catch( std::exception& ) blocks do not swallow it.
        struct ROL_Budget_Exhausted
        {
        };
    }    // namespace

#ifdef MORIS_HAVE_ROL
    namespace
    {
        //--------------------------------------------------------------------------------------------------------------
        // Small helpers to read/write a MORIS ADV std::vector out of a ROL::Vector (always a StdVector here).

        inline const std::vector< real >&
        rol_std_vector( const ROL::Vector< real >& aVector )
        {
            return *dynamic_cast< const ROL::StdVector< real >& >( aVector ).getVector();
        }

        inline std::vector< real >&
        rol_std_vector( ROL::Vector< real >& aVector )
        {
            return *dynamic_cast< ROL::StdVector< real >& >( aVector ).getVector();
        }

        //--------------------------------------------------------------------------------------------------------------
        // ROL objective adapter: forwards value()/gradient() into the MORIS Problem via the Algorithm.

        class Moris_ROL_Objective : public ROL::Objective< real >
        {
          private:
            Algorithm_ROL* mAlgorithm;

          public:
            explicit Moris_ROL_Objective( Algorithm_ROL* aAlgorithm )
                    : mAlgorithm( aAlgorithm )
            {
            }

            // ROL signals the disposition of a design point here; write the restart file only when
            // this point becomes the new accepted iterate (not for the many trial evaluations).
            void
            update( const ROL::Vector< real >& aX, ROL::UpdateType aType, int aIter ) override
            {
                if ( aType == ROL::UpdateType::Accept || aType == ROL::UpdateType::Initial )
                {
                    mAlgorithm->notify_accepted_step( rol_std_vector( aX ) );
                }
            }

            real
            value( const ROL::Vector< real >& aX, real& aTol ) override
            {
                mAlgorithm->ensure_criteria( rol_std_vector( aX ) );

                return mAlgorithm->cached_objectives()( 0 );
            }

            void
            gradient( ROL::Vector< real >& aG, const ROL::Vector< real >& aX, real& aTol ) override
            {
                mAlgorithm->ensure_gradients( rol_std_vector( aX ) );

                const Matrix< DDRMat >& tObjectiveGradient = mAlgorithm->cached_objective_gradients();

                std::vector< real >& tG = rol_std_vector( aG );
                for ( uint iADV = 0; iADV < tG.size(); ++iADV )
                {
                    tG[ iADV ] = tObjectiveGradient( iADV );
                }
            }
        };

        //--------------------------------------------------------------------------------------------------------------
        // ROL constraint adapter: the full MORIS constraint vector g(x) and its dense Jacobian
        // d g_i / d adv_j. Sign convention matches MORIS directly (no flip): equality constraints
        // are g_i = 0 and inequality constraints are g_i <= 0; the equality-vs-inequality distinction
        // is carried by the bound placed on the constraint output in rol_solve (eq -> [0,0],
        // ineq -> [-inf,0]), so this adapter is sign-convention agnostic.

        class Moris_ROL_Constraint : public ROL::Constraint< real >
        {
          private:
            Algorithm_ROL* mAlgorithm;

          public:
            explicit Moris_ROL_Constraint( Algorithm_ROL* aAlgorithm )
                    : mAlgorithm( aAlgorithm )
            {
            }

            void
            value( ROL::Vector< real >& aC, const ROL::Vector< real >& aX, real& aTol ) override
            {
                mAlgorithm->ensure_criteria( rol_std_vector( aX ) );

                const Matrix< DDRMat >& tConstraints = mAlgorithm->cached_constraints();

                std::vector< real >& tC = rol_std_vector( aC );
                for ( uint iCon = 0; iCon < tC.size(); ++iCon )
                {
                    tC[ iCon ] = tConstraints( iCon );
                }
            }

            void
            applyJacobian(
                    ROL::Vector< real >&       aJv,
                    const ROL::Vector< real >& aV,
                    const ROL::Vector< real >& aX,
                    real&                      aTol ) override
            {
                mAlgorithm->ensure_gradients( rol_std_vector( aX ) );

                const Matrix< DDRMat >& tJacobian = mAlgorithm->cached_constraint_gradients();    // nCon x nADV

                const std::vector< real >& tV  = rol_std_vector( aV );     // nADV
                std::vector< real >&       tJv = rol_std_vector( aJv );    // nCon

                for ( uint iCon = 0; iCon < tJv.size(); ++iCon )
                {
                    real tSum = 0.0;
                    for ( uint jADV = 0; jADV < tV.size(); ++jADV )
                    {
                        tSum += tJacobian( iCon, jADV ) * tV[ jADV ];
                    }
                    tJv[ iCon ] = tSum;
                }
            }

            void
            applyAdjointJacobian(
                    ROL::Vector< real >&       aAjv,
                    const ROL::Vector< real >& aV,
                    const ROL::Vector< real >& aX,
                    real&                      aTol ) override
            {
                mAlgorithm->ensure_gradients( rol_std_vector( aX ) );

                const Matrix< DDRMat >& tJacobian = mAlgorithm->cached_constraint_gradients();    // nCon x nADV

                const std::vector< real >& tV   = rol_std_vector( aV );      // nCon
                std::vector< real >&       tAjv = rol_std_vector( aAjv );    // nADV

                for ( uint jADV = 0; jADV < tAjv.size(); ++jADV )
                {
                    real tSum = 0.0;
                    for ( uint iCon = 0; iCon < tV.size(); ++iCon )
                    {
                        tSum += tJacobian( iCon, jADV ) * tV[ iCon ];
                    }
                    tAjv[ jADV ] = tSum;
                }
            }
        };


        //--------------------------------------------------------------------------------------------------------------
        // Persistent-secant PHR augmented Lagrangian (al_mode = "persistent"), for MORIS constraints g_i (type 0:
        // g_i = 0, type 1: g_i <= 0) with multipliers lambda_i and penalty rho:
        //   equality   : psi = lambda g + rho/2 g^2,                          dpsi/dg = lambda + rho g
        //   inequality : psi = ( max(0, lambda + rho g)^2 - lambda^2 ) / (2 rho), dpsi/dg = max(0, lambda + rho g)
        // L(x) = f(x) + sum_i psi_i is minimized subject to the ADV bounds only (Lin-More), then the multipliers are
        // updated outside (lambda <- dpsi/dg at the new iterate).

        class Moris_ROL_AL_Objective : public ROL::Objective< real >
        {
          private:
            Algorithm_ROL*             mAlgorithm;
            const std::vector< real >& mLambda;
            const real&                mRho;
            const Matrix< DDSMat >&    mTypes;

          public:
            Moris_ROL_AL_Objective(
                    Algorithm_ROL*             aAlgorithm,
                    const std::vector< real >& aLambda,
                    const real&                aRho,
                    const Matrix< DDSMat >&    aTypes )
                    : mAlgorithm( aAlgorithm )
                    , mLambda( aLambda )
                    , mRho( aRho )
                    , mTypes( aTypes )
            {
            }

            void
            update( const ROL::Vector< real >& aX, ROL::UpdateType aType, int aIter ) override
            {
                if ( aType == ROL::UpdateType::Accept || aType == ROL::UpdateType::Initial )
                {
                    mAlgorithm->notify_accepted_step( rol_std_vector( aX ) );
                }
            }

            real
            multiplier_slope( uint aIndex, real aConstraint ) const
            {
                real tSlope = mLambda[ aIndex ] + mRho * aConstraint;
                return ( mTypes( aIndex ) == 0 ) ? tSlope : std::max( 0.0, tSlope );
            }

            real
            value( const ROL::Vector< real >& aX, real& aTol ) override
            {
                mAlgorithm->ensure_criteria( rol_std_vector( aX ) );

                real                    tValue       = mAlgorithm->cached_objectives()( 0 );
                const Matrix< DDRMat >& tConstraints = mAlgorithm->cached_constraints();
                for ( uint iCon = 0; iCon < mLambda.size(); ++iCon )
                {
                    real tG = tConstraints( iCon );
                    if ( mTypes( iCon ) == 0 )
                    {
                        tValue += mLambda[ iCon ] * tG + 0.5 * mRho * tG * tG;
                    }
                    else
                    {
                        real tSlope = std::max( 0.0, mLambda[ iCon ] + mRho * tG );
                        tValue += ( tSlope * tSlope - mLambda[ iCon ] * mLambda[ iCon ] ) / ( 2.0 * mRho );
                    }
                }
                return tValue;
            }

            void
            gradient( ROL::Vector< real >& aG, const ROL::Vector< real >& aX, real& aTol ) override
            {
                mAlgorithm->ensure_gradients( rol_std_vector( aX ) );

                const Matrix< DDRMat >& tDObjective   = mAlgorithm->cached_objective_gradients();     // 1 x nADV
                const Matrix< DDRMat >& tDConstraints = mAlgorithm->cached_constraint_gradients();    // nCon x nADV
                const Matrix< DDRMat >& tConstraints  = mAlgorithm->cached_constraints();

                std::vector< real >& tG = rol_std_vector( aG );
                for ( uint iADV = 0; iADV < tG.size(); ++iADV )
                {
                    tG[ iADV ] = tDObjective( iADV );
                }
                for ( uint iCon = 0; iCon < mLambda.size(); ++iCon )
                {
                    real tSlope = this->multiplier_slope( iCon, tConstraints( iCon ) );
                    for ( uint iADV = 0; iADV < tG.size(); ++iADV )
                    {
                        tG[ iADV ] += tSlope * tDConstraints( iCon, iADV );
                    }
                }
            }
        };

        //--------------------------------------------------------------------------------------------------------------

    }    // namespace
#endif

    //----------------------------------------------------------------------------------------------------------------------

    Algorithm_ROL::Algorithm_ROL( const Parameter_List& aParameterList )
            : mStepType( aParameterList.get< std::string >( "step_type" ) )
            , mSubproblemModel( aParameterList.get< std::string >( "subproblem_model" ) )
            , mSubproblemSolver( aParameterList.get< std::string >( "subproblem_solver" ) )
            , mSecantType( aParameterList.get< std::string >( "secant_type" ) )
            , mSecantStorage( aParameterList.get< sint >( "secant_storage" ) )
            , mInitialTRRadius( aParameterList.get< real >( "initial_tr_radius" ) )
            , mMaxTRRadius( aParameterList.get< real >( "max_tr_radius" ) )
            , mUseAugmentedLagrangian( aParameterList.get< bool >( "use_augmented_lagrangian" ) )
            , mSubproblemIterationLimit( aParameterList.get< sint >( "subproblem_iteration_limit" ) )
            , mALInitialPenalty( aParameterList.get< real >( "al_initial_penalty" ) )
            , mALPenaltyGrowth( aParameterList.get< real >( "al_penalty_growth" ) )
            , mGradientTol( aParameterList.get< real >( "gradient_tol" ) )
            , mConstraintTol( aParameterList.get< real >( "constraint_tol" ) )
            , mStepTol( aParameterList.get< real >( "step_tol" ) )
            , mXmlFile( aParameterList.get< std::string >( "xml_file" ) )
            , mALMode( aParameterList.get< std::string >( "al_mode" ) )
            , mALInitialMultiplier( aParameterList.get< real >( "al_initial_multiplier" ) )
            , mALFeasibilityReduction( aParameterList.get< real >( "al_feasibility_reduction" ) )
            , mALMaxPenalty( aParameterList.get< real >( "al_max_penalty" ) )
            , mMaxForwardSolves( aParameterList.get< sint >( "max_forward_solves" ) )
    {
        MORIS_ERROR( mALMode == "rol" || mALMode == "persistent",
                "Algorithm_ROL: al_mode must be \"rol\" or \"persistent\", got \"%s\"", mALMode.c_str() );
#ifndef MORIS_HAVE_ROL
        MORIS_ERROR( false, "MORIS was compiled without ROL support" );
#endif

        mRestartIndex         = aParameterList.get< sint >( "restart_index" );
        mMaxIterationsInitial = aParameterList.get< sint >( "max_its" );
        mMaxIterations        = mMaxIterationsInitial;

        MORIS_LOG_INFO( "ROL initialized: step='%s', subproblem_model='%s', max_its=%d, al_mode='%s', max_forward_solves=%d",
                mStepType.c_str(), mSubproblemModel.c_str(), (int)mMaxIterations, mALMode.c_str(), (int)mMaxForwardSolves );
    }

    //----------------------------------------------------------------------------------------------------------------------

    Algorithm_ROL::~Algorithm_ROL()
    {
    }

    //----------------------------------------------------------------------------------------------------------------------

    uint
    Algorithm_ROL::solve(
            uint                       aCurrentOptAlgInd,
            std::shared_ptr< Problem > aOptProb )
    {
        // Trace optimization
        Tracer tTracer( "OptimizationAlgorithm", "ROL", "Solve" );

        // running status has to be wait when starting a solve
        mRunning = opt::Task::wait;

        mCurrentOptAlgInd = aCurrentOptAlgInd;    // set index of current optimization algorithm
        mProblem          = aOptProb;             // set the member variable mProblem to aOptProb

        // fresh solve: invalidate the evaluation cache and the forward-solve budget counter
        mLatest           = Eval_Point();
        mAccepted         = Eval_Point();
        mActivePoint      = &mLatest;
        mNumForwardSolves = 0;
        mLastAcceptedADVs.clear();

        // Set optimization iteration index for restart
        if ( mRestartIndex > 0 )
        {
            gLogger.set_opt_iteration( mRestartIndex );
        }

        // Solve optimization problem
        if ( par_rank() == 0 )
        {
            // Run ROL algorithm
            this->rol_solve();

            // Communicate that optimization has finished
            mRunning = opt::Task::exit;

            this->communicate_running_status();
        }
        else
        {
            // Run dummy solve (follows rank 0's forward/gradient broadcasts until exit)
            this->dummy_solve();
        }

        uint tOptIter = gLogger.get_opt_iteration();

        gLogger.set_iteration( "OPT", "Manager", "Perform", tOptIter );

        return tOptIter;
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ROL::notify_accepted_step( const std::vector< real >& aADVs )
    {
        // keep the accepted iterate's criteria / gradients (normally the latest evaluation is at aADVs)
        if ( mLatest.mHaveCriteria && mLatest.mADVs == aADVs )
        {
            mAccepted = mLatest;
        }
        mLastAcceptedADVs = aADVs;

        // grep-able accepted-iterate marker (fsolve = forward solve index of this design, 1-based)
        MORIS_LOG_INFO( "ROL_ACCEPT fsolve=%d", (int)( mAccepted.mADVs == aADVs ? mAccepted.mFSolve : -1 ) );

        Vector< real > tADVs( aADVs.size() );
        for ( uint iADV = 0; iADV < aADVs.size(); ++iADV )
        {
            tADVs( iADV ) = aADVs[ iADV ];
        }

        this->write_advs_to_file( tADVs );
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ROL::ensure_criteria( const std::vector< real >& aADVs )
    {
        // skip the forward solve if criteria are already current at these ADVs
        if ( mLatest.mHaveCriteria && aADVs == mLatest.mADVs )
        {
            mActivePoint = &mLatest;
            return;
        }
        if ( mAccepted.mHaveCriteria && aADVs == mAccepted.mADVs )
        {
            mActivePoint = &mAccepted;
            return;
        }

        // forward-solve budget (max_forward_solves): stop ROL; rol_solve returns the last accepted iterate
        if ( mMaxForwardSolves > 0 && mNumForwardSolves >= mMaxForwardSolves )
        {
            throw ROL_Budget_Exhausted();
        }

        // copy into a MORIS ADV vector
        Vector< real > tADVs( aADVs.size() );
        for ( uint iADV = 0; iADV < aADVs.size(); ++iADV )
        {
            tADVs( iADV ) = aADVs[ iADV ];
        }

        // Compute the design criteria. NOTE: unlike GCMMA/SQP we do NOT write a restart file here.
        // ROL evaluates the objective at many trial (rejected) points per iteration, so a per-eval
        // restart write would emit hundreds of large HDF5 files; the restart is written only on
        // ACCEPTED steps via notify_accepted_step (called from the objective adapter's update()).
        this->compute_design_criteria( tADVs );
        ++mNumForwardSolves;

        mLatest.mADVs          = aADVs;
        mLatest.mObjectives    = this->get_objectives();
        mLatest.mConstraints   = this->get_constraints();
        mLatest.mHaveCriteria  = true;
        mLatest.mHaveGradients = false;    // criteria changed -> cached gradients are stale
        mLatest.mFSolve        = mNumForwardSolves;
        mActivePoint           = &mLatest;
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ROL::ensure_gradients( const std::vector< real >& aADVs )
    {
        this->ensure_criteria( aADVs );

        if ( mActivePoint->mHaveGradients )
        {
            return;
        }

        // the adjoint reuses the forward state: if the point is only cached (the Problem was moved to a later
        // trial), re-solve the forward problem at aADVs first
        if ( mActivePoint == &mAccepted )
        {
            if ( mMaxForwardSolves > 0 && mNumForwardSolves >= mMaxForwardSolves )
            {
                throw ROL_Budget_Exhausted();
            }
            Vector< real > tRe( aADVs.size() );
            for ( uint iADV = 0; iADV < aADVs.size(); ++iADV )
            {
                tRe( iADV ) = aADVs[ iADV ];
            }
            this->compute_design_criteria( tRe );
            ++mNumForwardSolves;
            mLatest                = mAccepted;
            mLatest.mObjectives    = this->get_objectives();
            mLatest.mConstraints   = this->get_constraints();
            mLatest.mHaveGradients = false;
            mLatest.mFSolve        = mNumForwardSolves;
            mActivePoint           = &mLatest;
        }

        Vector< real > tADVs( aADVs.size() );
        for ( uint iADV = 0; iADV < aADVs.size(); ++iADV )
        {
            tADVs( iADV ) = aADVs[ iADV ];
        }

        this->compute_design_criteria_gradients( tADVs );

        mLatest.mObjectiveGradients  = this->get_objective_gradients();
        mLatest.mConstraintGradients = this->get_constraint_gradients();
        mLatest.mHaveGradients       = true;
        if ( mAccepted.mHaveCriteria && mAccepted.mADVs == aADVs )
        {
            mAccepted = mLatest;
        }
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ROL::rol_solve()
    {
#ifdef MORIS_HAVE_ROL
        if ( mALMode == "persistent" && mProblem->get_num_constraints() > 0 && mUseAugmentedLagrangian )
        {
            this->rol_solve_persistent_al();
            return;
        }

        const uint tNumADVs = mProblem->get_num_advs();

        // --- design vector and bounds as ROL StdVectors backed by std::vectors -------------------
        auto tAdvVec = ROL::makePtr< std::vector< real > >( tNumADVs );
        auto tLoVec  = ROL::makePtr< std::vector< real > >( tNumADVs );
        auto tHiVec  = ROL::makePtr< std::vector< real > >( tNumADVs );

        for ( uint iADV = 0; iADV < tNumADVs; ++iADV )
        {
            ( *tAdvVec )[ iADV ] = mProblem->get_advs()( iADV );
            ( *tLoVec )[ iADV ]  = mProblem->get_lower_bounds()( iADV );
            ( *tHiVec )[ iADV ]  = mProblem->get_upper_bounds()( iADV );
        }

        ROL::Ptr< ROL::Vector< real > > tX     = ROL::makePtr< ROL::StdVector< real > >( tAdvVec );
        ROL::Ptr< ROL::Vector< real > > tLower = ROL::makePtr< ROL::StdVector< real > >( tLoVec );
        ROL::Ptr< ROL::Vector< real > > tUpper = ROL::makePtr< ROL::StdVector< real > >( tHiVec );

        ROL::Ptr< ROL::BoundConstraint< real > > tBounds =
                ROL::makePtr< ROL::Bounds< real > >( tLower, tUpper );

        // --- objective adapter -------------------------------------------------------------------
        ROL::Ptr< ROL::Objective< real > > tObjective =
                ROL::makePtr< Moris_ROL_Objective >( this );

        // --- assemble the ROL problem (bound-constrained; general constraints added in M2) -------
        ROL::Ptr< ROL::Problem< real > > tRolProblem =
                ROL::makePtr< ROL::Problem< real > >( tObjective, tX );
        tRolProblem->addBoundConstraint( tBounds );

        // --- general constraints (Augmented-Lagrangian outer loop) -------------------------------
        const uint tNumConstraints = mProblem->get_num_constraints();
        const bool tUseConstraints = ( tNumConstraints > 0 ) && mUseAugmentedLagrangian;

        if ( tNumConstraints > 0 && !mUseAugmentedLagrangian )
        {
            MORIS_LOG_WARNING( "ROL: %d constraint(s) present but use_augmented_lagrangian=false; "
                               "solving bound-constrained and IGNORING the general constraints.",
                    (int)tNumConstraints );
        }

        if ( tUseConstraints )
        {
            // Bounds on the constraint output encode the MORIS constraint types:
            //   equality   (type 0): g_i = 0    -> [ 0, 0 ]
            //   inequality (type 1): g_i <= 0   -> [ -inf, 0 ]
            Matrix< DDSMat > tConstraintTypes = mProblem->get_constraint_types();

            auto tConLo = ROL::makePtr< std::vector< real > >( tNumConstraints );
            auto tConHi = ROL::makePtr< std::vector< real > >( tNumConstraints );
            for ( uint iCon = 0; iCon < tNumConstraints; ++iCon )
            {
                ( *tConHi )[ iCon ] = 0.0;
                ( *tConLo )[ iCon ] = ( tConstraintTypes( iCon ) == 0 ) ? 0.0 : -ROL::ROL_INF< real >();
            }

            ROL::Ptr< ROL::Vector< real > > tConLoV = ROL::makePtr< ROL::StdVector< real > >( tConLo );
            ROL::Ptr< ROL::Vector< real > > tConHiV = ROL::makePtr< ROL::StdVector< real > >( tConHi );
            ROL::Ptr< ROL::BoundConstraint< real > > tConBound =
                    ROL::makePtr< ROL::Bounds< real > >( tConLoV, tConHiV );

            ROL::Ptr< ROL::Constraint< real > > tConstraint =
                    ROL::makePtr< Moris_ROL_Constraint >( this );

            auto tMulVec = ROL::makePtr< std::vector< real > >( tNumConstraints, 0.0 );
            ROL::Ptr< ROL::Vector< real > > tMultiplier = ROL::makePtr< ROL::StdVector< real > >( tMulVec );

            tRolProblem->addConstraint( "constraint", tConstraint, tMultiplier, tConBound );
        }

        // --- parameter list: XML passthrough, else curated trust-region / Lin-More keys ----------
        ROL::ParameterList tParameterList;

        if ( !mXmlFile.empty() )
        {
            ROL::Ptr< ROL::ParameterList > tXmlParameters = ROL::getParametersFromXmlFile( mXmlFile );
            tParameterList                                = *tXmlParameters;
        }
        else
        {
            // Print ROL's iteration history (the Lin-More trust-region table: value, gnorm, snorm,
            // trust-region radius delta, #fval/#grad, tr_flag/iterCG) so the trust-region behavior
            // on the immersed problem can be harvested. Level 1 = per-iteration history.
            tParameterList.sublist( "General" ).set( "Output Level", 1 );

            ROL::ParameterList& tStep = tParameterList.sublist( "Step" );
            // With general constraints, drive an Augmented-Lagrangian outer loop whose inner
            // subproblem is the (bound-constrained) trust region below; otherwise a plain TR.
            tStep.set( "Type", tUseConstraints ? std::string( "Augmented Lagrangian" ) : mStepType );

            ROL::ParameterList& tTrustRegion = tStep.sublist( "Trust Region" );
            tTrustRegion.set( "Subproblem Model", mSubproblemModel );
            tTrustRegion.set( "Subproblem Solver", mSubproblemSolver );
            tTrustRegion.set( "Initial Radius", mInitialTRRadius );
            tTrustRegion.set( "Maximum Radius", mMaxTRRadius );

            if ( tUseConstraints )
            {
                ROL::ParameterList& tAugLag = tStep.sublist( "Augmented Lagrangian" );
                // ROL ignores "Initial Penalty Parameter" unless the default-penalty heuristic is switched off
                // (ROL_TypeG_AugmentedLagrangianAlgorithm_Def.hpp: "Use Default Initial Penalty Parameter"
                // defaults to true), which made al_initial_penalty a silent no-op.
                tAugLag.set( "Use Default Initial Penalty Parameter", false );
                tAugLag.set( "Initial Penalty Parameter", mALInitialPenalty );
                tAugLag.set( "Penalty Parameter Growth Factor", mALPenaltyGrowth );
                tAugLag.set( "Subproblem Step Type", "Trust Region" );
                // Cap the inner (bound-constrained TR) subproblem so each AL outer iteration does a
                // bounded number of forward/adjoint solves; without this the inner solver dominates
                // the run cost (the immersed forward+adjoint solve is expensive).
                tAugLag.set( "Subproblem Iteration Limit", (int)mSubproblemIterationLimit );
                // Print the inner (Lin-More TR) subproblem history — this is where the trust-region
                // radius trajectory and step accept/reject signals live (the applicability signals).
                tAugLag.set( "Print Intermediate Optimization History", true );
            }

            ROL::ParameterList& tSecant = tParameterList.sublist( "General" ).sublist( "Secant" );
            tSecant.set( "Type", mSecantType );
            tSecant.set( "Maximum Storage", (int)mSecantStorage );
            // Use the secant approximation AS the Hessian. Without this ROL finite-differences the
            // Hessian-vector product from the gradient, firing many extra (expensive) adjoint solves
            // per iteration; the immersed MORIS gradient is far too costly for that. Matches Plato.
            tSecant.set( "Use as Hessian", true );

            ROL::ParameterList& tStatus = tParameterList.sublist( "Status Test" );
            tStatus.set( "Iteration Limit", (int)mMaxIterations );
            tStatus.set( "Gradient Tolerance", mGradientTol );
            tStatus.set( "Constraint Tolerance", mConstraintTol );
            tStatus.set( "Step Tolerance", mStepTol );
        }

        // finalize and solve (ROL log to stdout so the trust-region trajectory can be harvested)
        tRolProblem->finalize( false, true, std::cout );

        ROL::Solver< real > tSolver( tRolProblem, tParameterList );
        try
        {
            tSolver.solve( std::cout );

            ROL::Ptr< const ROL::AlgorithmState< real > > tState = tSolver.getAlgorithmState();
            MORIS_LOG_INFO( "ROL finished: iterations=%d, objective=%.6e, gnorm=%.3e",
                    tState->iter, tState->value, tState->gnorm );
        }
        catch ( const ROL_Budget_Exhausted& )
        {
            MORIS_LOG_INFO( "ROL stopped: forward-solve budget max_forward_solves=%d reached; returning the last accepted iterate",
                    (int)mMaxForwardSolves );
            if ( !mLastAcceptedADVs.empty() )
            {
                *tAdvVec = mLastAcceptedADVs;
            }
        }

        // --- leave the MORIS Problem at the final design ----------------------------------------
        Vector< real > tFinalADVs( tNumADVs );
        for ( uint iADV = 0; iADV < tNumADVs; ++iADV )
        {
            tFinalADVs( iADV ) = ( *tAdvVec )[ iADV ];
        }
        this->compute_design_criteria( tFinalADVs );
#else
        MORIS_ERROR( false, "MORIS was compiled without ROL support" );
#endif
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ROL::rol_solve_persistent_al()
    {
#ifdef MORIS_HAVE_ROL
        const uint tNumADVs        = mProblem->get_num_advs();
        const uint tNumConstraints = mProblem->get_num_constraints();

        auto tAdvVec = ROL::makePtr< std::vector< real > >( tNumADVs );
        auto tLoVec  = ROL::makePtr< std::vector< real > >( tNumADVs );
        auto tHiVec  = ROL::makePtr< std::vector< real > >( tNumADVs );
        for ( uint iADV = 0; iADV < tNumADVs; ++iADV )
        {
            ( *tAdvVec )[ iADV ] = mProblem->get_advs()( iADV );
            ( *tLoVec )[ iADV ]  = mProblem->get_lower_bounds()( iADV );
            ( *tHiVec )[ iADV ]  = mProblem->get_upper_bounds()( iADV );
        }
        // start inside the bounds (Lin-More assumes a feasible start for the bound constraints)
        for ( uint iADV = 0; iADV < tNumADVs; ++iADV )
        {
            ( *tAdvVec )[ iADV ] = std::min( std::max( ( *tAdvVec )[ iADV ], ( *tLoVec )[ iADV ] ), ( *tHiVec )[ iADV ] );
        }

        ROL::StdVector< real >                   tX( tAdvVec );
        ROL::Ptr< ROL::BoundConstraint< real > > tBounds = ROL::makePtr< ROL::Bounds< real > >(
                ROL::makePtr< ROL::StdVector< real > >( tLoVec ), ROL::makePtr< ROL::StdVector< real > >( tHiVec ) );

        Matrix< DDSMat > tTypes = mProblem->get_constraint_types();

        // --- multipliers and penalty --------------------------------------------------------------
        std::vector< real > tLambda( tNumConstraints, 0.0 );
        real                tRho = mALInitialPenalty;

        this->ensure_gradients( *tAdvVec );
        const Matrix< DDRMat > tG0 = this->cached_constraints();
        if ( mALInitialMultiplier >= 0.0 )
        {
            std::fill( tLambda.begin(), tLambda.end(), mALInitialMultiplier );
        }
        else
        {
            // least-squares multiplier estimate at the start: min || df + sum lambda_i dg_i ||, then projected
            // (lambda_i >= 0 for inequalities); for one constraint lambda = -<df, dg> / <dg, dg>
            const Matrix< DDRMat >& tDf = this->cached_objective_gradients();
            const Matrix< DDRMat >& tDg = this->cached_constraint_gradients();
            for ( uint iCon = 0; iCon < tNumConstraints; ++iCon )
            {
                real tFG = 0.0, tGG = 0.0;
                for ( uint iADV = 0; iADV < tNumADVs; ++iADV )
                {
                    tFG += tDf( iADV ) * tDg( iCon, iADV );
                    tGG += tDg( iCon, iADV ) * tDg( iCon, iADV );
                }
                real tL         = ( tGG > 0.0 ) ? -tFG / tGG : 0.0;
                tLambda[ iCon ] = ( tTypes( iCon ) == 0 ) ? tL : std::max( 0.0, tL );
            }
        }

        // violation measure: equality |g|, inequality max(0, g)
        auto tViolation = [ & ]( const Matrix< DDRMat >& aG ) -> real {
            real tV = 0.0;
            for ( uint iCon = 0; iCon < tNumConstraints; ++iCon )
            {
                tV = std::max( tV, ( tTypes( iCon ) == 0 ) ? std::abs( aG( iCon ) ) : std::max( 0.0, aG( iCon ) ) );
            }
            return tV;
        };
        real tPrevViolation = tViolation( tG0 );

        // --- ROL parameter list for the bound-constrained inner solve (one secant for the whole run) ---------
        ROL::ParameterList tParameterList;
        if ( !mXmlFile.empty() )
        {
            tParameterList = *ROL::getParametersFromXmlFile( mXmlFile );
        }
        else
        {
            tParameterList.sublist( "General" ).set( "Output Level", 1 );
            tParameterList.sublist( "Step" ).set( "Type", mStepType == "Augmented Lagrangian" ? std::string( "Trust Region" ) : mStepType );
            ROL::ParameterList& tTrustRegion = tParameterList.sublist( "Step" ).sublist( "Trust Region" );
            tTrustRegion.set( "Subproblem Model", mSubproblemModel );
            tTrustRegion.set( "Subproblem Solver", mSubproblemSolver );
            tTrustRegion.set( "Maximum Radius", mMaxTRRadius );
            ROL::ParameterList& tSecantList = tParameterList.sublist( "General" ).sublist( "Secant" );
            tSecantList.set( "Type", mSecantType );
            tSecantList.set( "Maximum Storage", (int)mSecantStorage );
            tSecantList.set( "Use as Hessian", true );
        }
        tParameterList.sublist( "Status Test" ).set( "Iteration Limit", (int)mSubproblemIterationLimit );
        tParameterList.sublist( "Status Test" ).set( "Gradient Tolerance", mGradientTol );
        tParameterList.sublist( "Status Test" ).set( "Step Tolerance", mStepTol );

        ROL::Ptr< ROL::Secant< real > > tSecant = ROL::SecantFactory< real >( tParameterList );

        Moris_ROL_AL_Objective tALObjective( this, tLambda, tRho, tTypes );

        real tRadius = mInitialTRRadius;
        bool tBudget = false;
        uint tOuter  = 0;

        for ( tOuter = 0; tOuter < (uint)mMaxIterations; ++tOuter )
        {
            tParameterList.sublist( "Step" ).sublist( "Trust Region" ).set( "Initial Radius", tRadius );

            ROL::Ptr< ROL::TypeB::Algorithm< real > > tAlgorithm =
                    ROL::TypeB::AlgorithmFactory< real >( tParameterList, tSecant );
            int  tInnerIterations = 0;
            real tGNorm           = 0.0;
            try
            {
                tAlgorithm->run( tX, tALObjective, *tBounds, std::cout );
                tInnerIterations = tAlgorithm->getState()->iter;
                tGNorm           = tAlgorithm->getState()->gnorm;
                tRadius          = std::max( tAlgorithm->getState()->searchSize, 1.0e-8 * mInitialTRRadius );
            }
            catch ( const ROL_Budget_Exhausted& )
            {
                tBudget = true;
            }

            // continue from the last accepted iterate (the inner solve may have ended on a rejected trial)
            if ( !mLastAcceptedADVs.empty() )
            {
                *tAdvVec = mLastAcceptedADVs;
            }
            if ( tBudget )
            {
                break;
            }

            // multiplier / penalty update at the accepted iterate (criteria cached -> no extra solve)
            this->ensure_criteria( *tAdvVec );
            const Matrix< DDRMat > tG         = this->cached_constraints();
            const real             tObjective = this->cached_objectives()( 0 );
            real                   tViol      = tViolation( tG );
            real                   tDLambda   = 0.0;    // multiplier change = complementarity / stationarity residual
            for ( uint iCon = 0; iCon < tNumConstraints; ++iCon )
            {
                real tNew       = tALObjective.multiplier_slope( iCon, tG( iCon ) );
                tDLambda        = std::max( tDLambda, std::abs( tNew - tLambda[ iCon ] ) / std::max( 1.0, std::abs( tNew ) ) );
                tLambda[ iCon ] = tNew;
            }
            if ( tViol > mALFeasibilityReduction * tPrevViolation && tViol > mConstraintTol )
            {
                tRho = std::min( tRho * mALPenaltyGrowth, mALMaxPenalty );
            }
            tPrevViolation = tViol;

            MORIS_LOG_INFO( "ROLAL outer=%d fsolves=%d inner_its=%d obj=%.8e violation=%.3e lambda0=%.6e dlambda=%.3e rho=%.3e radius=%.3e gnorm=%.3e",
                    (int)tOuter + 1, (int)mNumForwardSolves, tInnerIterations, tObjective, tViol,
                    tNumConstraints > 0 ? tLambda[ 0 ] : 0.0, tDLambda, tRho, tRadius, tGNorm );

            if ( tViol <= mConstraintTol && tDLambda <= mConstraintTol && tGNorm <= mGradientTol )
            {
                MORIS_LOG_INFO( "ROLAL converged: violation %.3e, relative multiplier change %.3e (<= %.3e) and projected AL gradient %.3e <= %.3e",
                        tViol, tDLambda, mConstraintTol, tGNorm, mGradientTol );
                break;
            }
        }

        if ( tBudget )
        {
            MORIS_LOG_INFO( "ROLAL stopped: forward-solve budget max_forward_solves=%d reached at outer %d; returning the last accepted iterate",
                    (int)mMaxForwardSolves, (int)tOuter + 1 );
        }

        // leave the MORIS Problem at the final (last accepted) design
        Vector< real > tFinalADVs( tNumADVs );
        for ( uint iADV = 0; iADV < tNumADVs; ++iADV )
        {
            tFinalADVs( iADV ) = ( *tAdvVec )[ iADV ];
        }
        this->compute_design_criteria( tFinalADVs );
#else
        MORIS_ERROR( false, "MORIS was compiled without ROL support" );
#endif
    }

    //----------------------------------------------------------------------------------------------------------------------

}    // namespace moris::opt
