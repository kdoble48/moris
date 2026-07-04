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
#endif

using namespace moris;

namespace moris::opt
{
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

                return mAlgorithm->get_objectives()( 0 );
            }

            void
            gradient( ROL::Vector< real >& aG, const ROL::Vector< real >& aX, real& aTol ) override
            {
                mAlgorithm->ensure_gradients( rol_std_vector( aX ) );

                const Matrix< DDRMat >& tObjectiveGradient = mAlgorithm->get_objective_gradients();

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

                const Matrix< DDRMat >& tConstraints = mAlgorithm->get_constraints();

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

                const Matrix< DDRMat >& tJacobian = mAlgorithm->get_constraint_gradients();    // nCon x nADV

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

                const Matrix< DDRMat >& tJacobian = mAlgorithm->get_constraint_gradients();    // nCon x nADV

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
    {
#ifndef MORIS_HAVE_ROL
        MORIS_ERROR( false, "MORIS was compiled without ROL support" );
#endif

        mRestartIndex         = aParameterList.get< sint >( "restart_index" );
        mMaxIterationsInitial = aParameterList.get< sint >( "max_its" );
        mMaxIterations        = mMaxIterationsInitial;

        MORIS_LOG_INFO( "ROL initialized: step='%s', subproblem_model='%s', max_its=%d",
                mStepType.c_str(), mSubproblemModel.c_str(), (int)mMaxIterations );
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

        // fresh solve: invalidate the evaluation cache
        mHaveCriteria  = false;
        mHaveGradients = false;

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
        if ( mHaveCriteria && aADVs == mCriteriaADVs )
        {
            return;
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

        mCriteriaADVs  = aADVs;
        mHaveCriteria  = true;
        mHaveGradients = false;    // criteria changed -> cached gradients are stale
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ROL::ensure_gradients( const std::vector< real >& aADVs )
    {
        // the adjoint reuses the forward state, so make sure criteria are current first
        this->ensure_criteria( aADVs );

        // skip the adjoint solve if gradients are already current at these ADVs
        if ( mHaveGradients && aADVs == mGradientADVs )
        {
            return;
        }

        Vector< real > tADVs( aADVs.size() );
        for ( uint iADV = 0; iADV < aADVs.size(); ++iADV )
        {
            tADVs( iADV ) = aADVs[ iADV ];
        }

        this->compute_design_criteria_gradients( tADVs );

        mGradientADVs  = aADVs;
        mHaveGradients = true;
    }

    //----------------------------------------------------------------------------------------------------------------------

    void
    Algorithm_ROL::rol_solve()
    {
#ifdef MORIS_HAVE_ROL
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
        tSolver.solve( std::cout );

        ROL::Ptr< const ROL::AlgorithmState< real > > tState = tSolver.getAlgorithmState();
        MORIS_LOG_INFO( "ROL finished: iterations=%d, objective=%.6e, gnorm=%.3e",
                tState->iter, tState->value, tState->gnorm );

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

}    // namespace moris::opt
