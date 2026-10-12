/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * opt_test.cpp
 *
 */

#include <catch.hpp>

#include <cmath>
#include <fstream>
#include <string>
#include <cstdio>

#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_OPT_Manager.hpp"
#include "fn_OPT_create_interface.hpp"
#include "cl_OPT_Problem_User_Defined.hpp"
#include "cl_OPT_Interface_User_Defined.hpp"
#include "cl_OPT_Algorithm_MMA.hpp"
#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "fn_OPT_Rosenbrock.hpp"
#include "fn_OPT_Test_Interface.hpp"

namespace moris::opt
{
    // ---------------------------------------------------------------------------------------------------------
    // Bimodal test problem for the gradient-explosion clip: one exploding sensitivity (1e8),
    // three healthy entries (~0.01), one exactly-zero (inactive) entry -- the structure of a
    // level-set design gradient with a small-cut-cell explosion.

    namespace bimodal
    {
        void initialize(
                Vector< real >& aADVs,
                Vector< real >& aLowerBounds,
                Vector< real >& aUpperBounds )
        {
            if ( par_rank() == 0 )
            {
                aADVs        = { 0.5, 0.5, 0.5, 0.5, 0.5 };
                aLowerBounds = { 0.0, 0.0, 0.0, 0.0, 0.0 };
                aUpperBounds = { 1.0, 1.0, 1.0, 1.0, 1.0 };
            }
        }

        Vector< real > get_criteria( const Vector< real >& aADVs )
        {
            if ( par_rank() == 0 )
            {
                return { 1.0 };
            }
            return { {} };
        }

        Matrix< DDRMat > get_dcriteria_dadv( const Vector< real >& aADVs )
        {
            if ( par_rank() == 0 )
            {
                return { { 1.0e8, 0.01, 0.012, 0.008, 0.0 } };
            }
            return { {} };
        }

        Matrix< DDSMat > get_constraint_types()
        {
            return { { 1 } };
        }

        Matrix< DDRMat > compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { aCriteria( 0 ) } };
        }

        Matrix< DDRMat > compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { 0.0 } };
        }

        Matrix< DDRMat > compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return Matrix< DDRMat >( 1, 5, 0.0 );
        }

        Matrix< DDRMat > compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { 1.0 } };
        }

        Matrix< DDRMat > compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return Matrix< DDRMat >( 1, 5, 0.0 );
        }

        Matrix< DDRMat > compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { 0.0 } };
        }
    }    // namespace bimodal

    // ---------------------------------------------------------------------------------------------------------
    // Binding-constraint quadratic for the persistent-AL ROL mode (task 34):
    //   min (x0 - 2)^2 + (x1 - 2)^2  s.t.  x0 + x1 - 2 <= 0,  -3 <= x <= 3
    // The unconstrained minimum (2,2) is infeasible; the constrained minimizer is (1,1), f* = 2, lambda* = 2.
    // Criteria = the ADVs themselves.

    namespace binding
    {
        void initialize(
                Vector< real >& aADVs,
                Vector< real >& aLowerBounds,
                Vector< real >& aUpperBounds )
        {
            if ( par_rank() == 0 )
            {
                aADVs        = { 0.0, -1.0 };
                aLowerBounds = { -3.0, -3.0 };
                aUpperBounds = { 3.0, 3.0 };
            }
        }

        Vector< real > get_criteria( const Vector< real >& aADVs )
        {
            if ( par_rank() == 0 )
            {
                return { aADVs( 0 ), aADVs( 1 ) };
            }
            return { {} };
        }

        Matrix< DDRMat > get_dcriteria_dadv( const Vector< real >& aADVs )
        {
            if ( par_rank() == 0 )
            {
                return { { 1.0, 0.0 }, { 0.0, 1.0 } };
            }
            return { {} };
        }

        Matrix< DDSMat > get_constraint_types()
        {
            return Matrix< DDSMat >( 1, 1, 1 );
        }

        Matrix< DDRMat > compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { std::pow( aCriteria( 0 ) - 2.0, 2 ) + std::pow( aCriteria( 1 ) - 2.0, 2 ) } };
        }

        Matrix< DDRMat > compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { aCriteria( 0 ) + aCriteria( 1 ) - 2.0 } };
        }

        Matrix< DDRMat > compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return Matrix< DDRMat >( 1, 2, 0.0 );
        }

        Matrix< DDRMat > compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { 2.0 * ( aCriteria( 0 ) - 2.0 ), 2.0 * ( aCriteria( 1 ) - 2.0 ) } };
        }

        Matrix< DDRMat > compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return Matrix< DDRMat >( 1, 2, 0.0 );
        }

        Matrix< DDRMat > compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { 1.0, 1.0 } };
        }

        std::shared_ptr< Problem >
        make_problem()
        {
            Parameter_List tProblemParameterList = moris::prm::create_opt_problem_parameter_list();

            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize, &get_criteria, &get_dcriteria_dadv );

            return std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types,
                    &compute_objectives,
                    &compute_constraints,
                    &compute_dobjective_dadv,
                    &compute_dobjective_dcriteria,
                    &compute_dconstraint_dadv,
                    &compute_dconstraint_dcriteria );
        }
    }    // namespace binding

    TEST_CASE( "[optimization]" )
    {

        // ---------------------------------------------------------------------------------------------------------

        // Native in-tree GCMMA port (cl_OPT_Algorithm_MMA). No TPL dependency, so not guarded by
        // MORIS_HAVE_GCMMA. Validates the native solver reaches the Rosenbrock minimum like the TPL.
        SECTION( "MMA native" )
        {
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();

            tAlgorithmParameterList.set( "norm_drop", 2.5e-5 );
            tAlgorithmParameterList.set( "asymp_adaptc", 1.1 );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 2E-7 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-4 );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        // Oscillation-adaptive per-ADV move limit (move_adapt): the update rule itself, then the native solver
        // with the rule on (Rosenbrock minimum reached), and shrink = grow = 1 reproducing the fixed box exactly.
        SECTION( "MMA adaptive move limit rule" )
        {
            // range 1. ADV 0 flips (+0.1, -0.05) on the bound: shrink relative to the step taken, 0.7 * min(0.08, 0.05);
            // ADV 1 agrees on the bound: grow; ADV 2 does not move; ADV 3 flips at the floor; ADV 4 grows into the cap;
            // ADV 5 flips but its steps were interior (not on the move bound): unchanged
            std::vector< double > tMove  = { 0.08, 0.01, 0.01, 1.0e-4, 0.02, 0.01 };
            std::vector< double > tXOld2 = { 0.0, 0.0, 0.5, 0.0, 0.0, 0.0 };
            std::vector< double > tXOld1 = { 0.1, 0.1, 0.5, 0.1, 0.1, 0.005 };
            std::vector< double > tX     = { 0.05, 0.2, 0.5, 0.0, 0.2, 0.0 };
            std::vector< double > tRange( 6, 1.0 );
            std::vector< int >    tHit   = { 1, 1, 1, 1, 1, 0 };

            int tShrunk = Algorithm_MMA::adapt_move_limits( tMove, tX, tXOld1, tXOld2, tRange, tHit, 0.7, 1.2, 1.0e-4, 0.02 );

            REQUIRE( tShrunk == 2 );
            REQUIRE( std::abs( tMove[ 0 ] - 0.02 ) < 1e-15 );     // 0.7 * 0.05 = 0.035, then the cap 0.02
            REQUIRE( std::abs( tMove[ 1 ] - 0.012 ) < 1e-15 );    // grow
            REQUIRE( std::abs( tMove[ 2 ] - 0.01 ) < 1e-15 );     // zero increment: unchanged
            REQUIRE( std::abs( tMove[ 3 ] - 1.0e-4 ) < 1e-15 );   // clamped at the floor
            REQUIRE( std::abs( tMove[ 4 ] - 0.02 ) < 1e-15 );     // clamped at the cap
            REQUIRE( std::abs( tMove[ 5 ] - 0.01 ) < 1e-15 );     // interior steps: gated out

            // shrink is relative to the step actually taken: range 2, step 0.02 (= 0.01 of range) on a 0.05 limit
            std::vector< double > tM2 = { 0.05 }, tA = { 0.0 }, tB = { 0.02 }, tC = { 0.0 }, tR2 = { 2.0 };
            std::vector< int >    tH2 = { 1 };
            Algorithm_MMA::adapt_move_limits( tM2, tC, tB, tA, tR2, tH2, 0.5, 1.2, 1.0e-4, 0.05 );
            REQUIRE( std::abs( tM2[ 0 ] - 0.005 ) < 1e-15 );
        }

        SECTION( "MMA rise guard test" )
        {
            std::vector< double > tFeas = { -0.01, 0.001 }, tInfeas = { 0.01 };
            REQUIRE( Algorithm_MMA::guard_rejects( 1.0, 1.2, tFeas, 0.1, 2e-3 ) );       // feasible, +20% > 10%
            REQUIRE( !Algorithm_MMA::guard_rejects( 1.0, 1.05, tFeas, 0.1, 2e-3 ) );     // +5%: accepted
            REQUIRE( !Algorithm_MMA::guard_rejects( 1.0, 1.5, tInfeas, 0.1, 2e-3 ) );    // infeasible march: never guarded
            REQUIRE( !Algorithm_MMA::guard_rejects( 1.0, 1.5, tFeas, -1.0, 2e-3 ) );     // off
            REQUIRE( Algorithm_MMA::guard_rejects( -1.0, -0.8, tFeas, 0.1, 2e-3 ) );     // relative to |f0|
        }

        SECTION( "MMA native adaptive move limit" )
        {
            auto tRun = []( bool aAdapt, real aShrink, real aGrow, sint aMaxIts, bool aGate = true, real aGuard = -1.0 ) -> std::pair< real, Vector< real > >
            {
                Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
                Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();
                tAlgorithmParameterList.set( "norm_drop", 2.5e-5 );
                tAlgorithmParameterList.set( "asymp_adaptc", 1.1 );
                tAlgorithmParameterList.set( "max_its", aMaxIts );
                tAlgorithmParameterList.set( "move_adapt", aAdapt );
                tAlgorithmParameterList.set( "move_adapt_shrink", aShrink );
                tAlgorithmParameterList.set( "move_adapt_grow", aGrow );
                tAlgorithmParameterList.set( "move_adapt_min", 1.0e-6 );
                tAlgorithmParameterList.set( "move_adapt_gate", aGate );
                tAlgorithmParameterList.set( "move_guard_rise", aGuard );

                std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                        &initialize_rosenbrock,
                        &get_criteria_rosenbrock,
                        &get_dcriteria_dadv_rosenbrock );
                std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &get_constraint_types_rosenbrock,
                        &compute_objectives_rosenbrock,
                        &compute_constraints_rosenbrock,
                        &compute_dobjective_dadv_rosenbrock,
                        &compute_dobjective_dcriteria_rosenbrock,
                        &compute_dconstraint_dadv_rosenbrock,
                        &compute_dconstraint_dcriteria_rosenbrock );

                Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
                tAlgorithms.add_parameter_list( tAlgorithmParameterList );
                Manager tManager( tAlgorithms, tProblem );
                tManager.perform();
                return { tManager.get_objectives()( 0 ), tManager.get_advs() };
            };

            // the rule on (gated and ungated): still reaches the Rosenbrock minimum
            for ( bool tGate : { true, false } )
            {
                auto [ tObj, tADVs ] = tRun( true, 0.7, 1.2, 100, tGate, 0.5 );    // with the rise guard
                if ( par_rank() == 0 )
                {
                    REQUIRE( std::abs( tObj ) < 2E-7 );
                    for ( auto iADV : tADVs )
                    {
                        REQUIRE( std::abs( iADV - 1.0 ) < 1E-4 );
                    }
                }
            }

            // shrink = grow = 1: every move limit stays step_size, so the trajectory equals the fixed box exactly
            auto [ tObjFix, tADVsFix ] = tRun( false, 0.7, 1.2, 12 );
            auto [ tObjNeu, tADVsNeu ] = tRun( true, 1.0, 1.0, 12 );
            if ( par_rank() == 0 )
            {
                REQUIRE( tObjFix == tObjNeu );
                for ( uint i = 0; i < tADVsFix.size(); ++i )
                {
                    REQUIRE( tADVsFix( i ) == tADVsNeu( i ) );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        // Gradient clip gate: off by default (raw gradient passes through untouched), and when
        // enabled via "grad_clip_factor" the exploding entry is capped at factor x the median of
        // the nonzero |entries| while healthy and zero entries are left alone.
        SECTION( "Gradient clip gate" )
        {
            for ( bool tClipOn : { false, true } )
            {
                Parameter_List tProblemParameterList = moris::prm::create_opt_problem_parameter_list();

                if ( tClipOn )
                {
                    tProblemParameterList.set( "grad_clip_factor", 20.0 );
                }

                std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                        &bimodal::initialize,
                        &bimodal::get_criteria,
                        &bimodal::get_dcriteria_dadv );

                std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                        tProblemParameterList,
                        tInterface,
                        &bimodal::get_constraint_types,
                        &bimodal::compute_objectives,
                        &bimodal::compute_constraints,
                        &bimodal::compute_dobjective_dadv,
                        &bimodal::compute_dobjective_dcriteria,
                        &bimodal::compute_dconstraint_dadv,
                        &bimodal::compute_dconstraint_dcriteria );

                tProblem->initialize();

                Vector< real > tADVs = tProblem->get_advs();
                tProblem->compute_design_criteria( tADVs );
                tProblem->compute_design_criteria_gradients( tADVs );

                if ( par_rank() == 0 )
                {
                    const Matrix< DDRMat >& tObjGrad = tProblem->get_objective_gradients();

                    if ( tClipOn )
                    {
                        // median of nonzero |entries| {1e8, 0.01, 0.012, 0.008} is 0.012 -> clip at 20 x 0.012
                        REQUIRE( tObjGrad( 0, 0 ) == Approx( 20.0 * 0.012 ) );
                    }
                    else
                    {
                        // clip disabled: exploding entry passes through untouched
                        REQUIRE( tObjGrad( 0, 0 ) == Approx( 1.0e8 ) );
                    }

                    // healthy and inactive entries are never modified
                    REQUIRE( tObjGrad( 0, 1 ) == Approx( 0.01 ) );
                    REQUIRE( tObjGrad( 0, 2 ) == Approx( 0.012 ) );
                    REQUIRE( tObjGrad( 0, 3 ) == Approx( 0.008 ) );
                    REQUIRE( tObjGrad( 0, 4 ) == 0.0 );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

#ifdef MORIS_HAVE_GCMMA
        SECTION( "GCMMA" )
        {
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_gcmma_parameter_list();

            tAlgorithmParameterList.set( "version", 1 );
            tAlgorithmParameterList.set( "norm_drop", 2.5e-5 );
            tAlgorithmParameterList.set( "asymp_adaptc", 1.1 );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 2E-7 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-4 );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        SECTION( "GCMMA 07" )
        {
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_gcmma_parameter_list();

            tAlgorithmParameterList.set( "version", 2 );
            tAlgorithmParameterList.set( "norm_drop", 1e-6 );
            tAlgorithmParameterList.set( "max_inner_its", 3 );
            tAlgorithmParameterList.set( "max_its", 11 );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 2E-7 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-4 );
                }
            }
        }
#endif

        // ---------------------------------------------------------------------------------------------------------

#ifdef MORIS_HAVE_SNOPT
        SECTION( "SQP" )
        {
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_sqp_parameter_list();

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 5E-4 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-3 );
                }
            }
        }
#endif

        // ---------------------------------------------------------------------------------------------------------

#ifdef MORIS_HAVE_LBFGS
        SECTION( "LBFGS" )
        {
            // This optimization problem does not use the constraints!
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_lbfgs_parameter_list();

            tAlgorithmParameterList.set( "max_its", 10 );
            tAlgorithmParameterList.set( "num_corr", 5 );
            tAlgorithmParameterList.set( "norm_drop", 1.0e-12 );

            tAlgorithmParameterList.set( "step_size", "0.2,0.1" );
            tAlgorithmParameterList.set( "outer_iteration_index", "0,1" );
            tAlgorithmParameterList.set( "number_inner_iterations", "10,20" );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 2E-7 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-4 );
                }
            }
        }
#endif

        // ---------------------------------------------------------------------------------------------------------

#ifdef MORIS_HAVE_ROL
        // Trilinos ROL trust-region (Lin-More) bridge, BOUND-CONSTRAINED. With
        // use_augmented_lagrangian = false the general constraints are ignored (like the LBFGS
        // case above); the bound-constrained Rosenbrock minimum is (1,1).
        SECTION( "ROL bound-constrained" )
        {
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_rol_parameter_list();

            tAlgorithmParameterList.set( "use_augmented_lagrangian", false );
            tAlgorithmParameterList.set( "max_its", 200 );
            tAlgorithmParameterList.set( "initial_tr_radius", 1.0 );
            tAlgorithmParameterList.set( "gradient_tol", 1.0e-12 );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 1E-6 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-3 );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        // ROL with the AUGMENTED-LAGRANGIAN outer loop over the two (inequality) Rosenbrock
        // constraints. Both constraints are active (= 0) at the unconstrained minimum (1,1), which
        // is therefore also the constrained minimizer -- so the AL path must land at (1,1) too.
        SECTION( "ROL augmented Lagrangian" )
        {
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_rol_parameter_list();

            tAlgorithmParameterList.set( "max_its", 100 );
            tAlgorithmParameterList.set( "initial_tr_radius", 1.0 );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 1E-4 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-2 );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        // Persistent-secant PHR augmented Lagrangian (al_mode = "persistent") on a BINDING constraint: the
        // unconstrained minimum is infeasible, so the multiplier update must actually drive feasibility.
        SECTION( "ROL persistent augmented Lagrangian (binding constraint)" )
        {
            Parameter_List tAlgorithmParameterList = moris::prm::create_rol_parameter_list();
            tAlgorithmParameterList.set( "al_mode", "persistent" );
            tAlgorithmParameterList.set( "max_its", 30 );
            tAlgorithmParameterList.set( "subproblem_iteration_limit", 10 );
            tAlgorithmParameterList.set( "initial_tr_radius", 1.0 );
            tAlgorithmParameterList.set( "al_initial_penalty", 1.0 );
            tAlgorithmParameterList.set( "al_penalty_growth", 10.0 );
            tAlgorithmParameterList.set( "gradient_tol", 1.0e-8 );
            tAlgorithmParameterList.set( "constraint_tol", 1.0e-8 );

            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, binding::make_problem() );
            tManager.perform();

            if ( par_rank() == 0 )
            {
                Vector< real > tADVs = tManager.get_advs();
                REQUIRE( std::abs( tADVs( 0 ) - 1.0 ) < 1E-4 );
                REQUIRE( std::abs( tADVs( 1 ) - 1.0 ) < 1E-4 );
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) - 2.0 ) < 1E-4 );
            }
        }

        // Adam with exact constraint projection on the binding quadratic: the iterate stays on x0 + x1 = 2 (linear
        // constraint -> the linearized projection is exact) and the decaying rate converges to (1,1).
        SECTION( "Adam projected (binding constraint)" )
        {
            Parameter_List tAlgorithmParameterList = moris::prm::create_adam_proj_parameter_list();
            tAlgorithmParameterList.set( "max_its", 400 );
            tAlgorithmParameterList.set( "learning_rate", 0.2 );
            tAlgorithmParameterList.set( "lr_decay_final", 0.001 );
            tAlgorithmParameterList.set( "projection_tol", 1.0e-10 );

            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, binding::make_problem() );
            tManager.perform();

            if ( par_rank() == 0 )
            {
                Vector< real > tADVs = tManager.get_advs();
                REQUIRE( std::abs( tADVs( 0 ) + tADVs( 1 ) - 2.0 ) < 1E-8 );    // exactly on the constraint
                REQUIRE( std::abs( tADVs( 0 ) - 1.0 ) < 2E-3 );
                REQUIRE( std::abs( tADVs( 1 ) - 1.0 ) < 2E-3 );
            }
        }

        // max_forward_solves stops ROL cleanly and leaves the Problem at the last accepted iterate (both modes).
        SECTION( "ROL forward-solve budget" )
        {
            for ( std::string tMode : { "rol", "persistent" } )
            {
                Parameter_List tAlgorithmParameterList = moris::prm::create_rol_parameter_list();
                tAlgorithmParameterList.set( "al_mode", tMode );
                tAlgorithmParameterList.set( "max_its", 100 );
                tAlgorithmParameterList.set( "subproblem_iteration_limit", 10 );
                tAlgorithmParameterList.set( "initial_tr_radius", 1.0 );
                tAlgorithmParameterList.set( "max_forward_solves", 3 );

                Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
                tAlgorithms.add_parameter_list( tAlgorithmParameterList );
                Manager tManager( tAlgorithms, binding::make_problem() );
                tManager.perform();

                if ( par_rank() == 0 )
                {
                    // moved off the start (0, -1) toward the minimizer, and stopped well before converging
                    Vector< real > tADVs = tManager.get_advs();
                    real           tF    = tManager.get_objectives()( 0 );
                    REQUIRE( tF < 13.0 );      // f(start) = 4 + 9 = 13
                    REQUIRE( std::abs( tADVs( 0 ) - 1.0 ) + std::abs( tADVs( 1 ) - 1.0 ) > 1E-6 );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------

        // ROL driven by a full ROL/Teuchos ParameterList XML (the "xml_file" escape hatch, e.g.
        // Plato's rol_inputs.xml). Confirms the passthrough path parses the XML and drives ROL;
        // bound-constrained so the XML "Trust Region" step matches (constraints off).
        SECTION( "ROL xml passthrough" )
        {
            // Write a minimal ROL parameter XML (rank 0 writes; the solver runs on rank 0)
            const std::string tXmlPath = "./rol_passthrough_test.xml";
            if ( par_rank() == 0 )
            {
                std::ofstream tXml( tXmlPath );
                tXml << "<ParameterList>\n"
                     << "  <ParameterList name=\"Step\">\n"
                     << "    <Parameter name=\"Type\" type=\"string\" value=\"Trust Region\"/>\n"
                     << "    <ParameterList name=\"Trust Region\">\n"
                     << "      <Parameter name=\"Subproblem Model\" type=\"string\" value=\"Lin-More\"/>\n"
                     << "      <Parameter name=\"Subproblem Solver\" type=\"string\" value=\"Truncated CG\"/>\n"
                     << "      <Parameter name=\"Initial Radius\" type=\"double\" value=\"1.0\"/>\n"
                     << "    </ParameterList>\n"
                     << "  </ParameterList>\n"
                     << "  <ParameterList name=\"General\">\n"
                     << "    <ParameterList name=\"Secant\">\n"
                     << "      <Parameter name=\"Type\" type=\"string\" value=\"Limited-Memory BFGS\"/>\n"
                     << "      <Parameter name=\"Use as Hessian\" type=\"bool\" value=\"true\"/>\n"
                     << "    </ParameterList>\n"
                     << "  </ParameterList>\n"
                     << "  <ParameterList name=\"Status Test\">\n"
                     << "    <Parameter name=\"Iteration Limit\" type=\"int\" value=\"200\"/>\n"
                     << "    <Parameter name=\"Gradient Tolerance\" type=\"double\" value=\"1.0e-12\"/>\n"
                     << "  </ParameterList>\n"
                     << "</ParameterList>\n";
            }
            barrier();

            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_rol_parameter_list();

            tAlgorithmParameterList.set( "use_augmented_lagrangian", false );
            tAlgorithmParameterList.set( "xml_file", tXmlPath );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::abs( tManager.get_objectives()( 0 ) ) < 1E-6 );    // check value of objective
                Vector< real > tADVs = tManager.get_advs();
                for ( auto iADV : tADVs )
                {
                    REQUIRE( std::abs( iADV - 1.0 ) < 1E-3 );
                }

                std::remove( tXmlPath.c_str() );
            }
        }
#endif

        // ---------------------------------------------------------------------------------------------------------

        SECTION( "Sweep" )
        {
            // Set up default parameter lists
            Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
            Parameter_List tAlgorithmParameterList = moris::prm::create_sweep_parameter_list();

            // Change parameters
            tAlgorithmParameterList.set( "hdf5_path", "./sweep.hdf5" );
            tAlgorithmParameterList.set( "num_evaluations_per_adv", "3, 2" );
            tAlgorithmParameterList.set( "finite_difference_type", "all" );
            tAlgorithmParameterList.set( "finite_difference_epsilons", "0.001, 0.01; 0.00001, 0.00001" );
            tAlgorithmParameterList.set( "finite_difference_adv_indices", "0,1" );

            // Create interface
            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    &get_dcriteria_dadv_rosenbrock );

            // Create Problem
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    &compute_dobjective_dadv_rosenbrock,
                    &compute_dobjective_dcriteria_rosenbrock,
                    &compute_dconstraint_dadv_rosenbrock,
                    &compute_dconstraint_dcriteria_rosenbrock );

            // Create manager
            Submodule_Parameter_Lists tAlgorithms( "Algorithms" );
            tAlgorithms.add_parameter_list( tAlgorithmParameterList );
            Manager tManager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();

            // Sweep without sensitivities
            tAlgorithms.set( "evaluate_objective_gradients", false );
            tAlgorithms.set( "evaluate_constraint_gradients", false );

            // Create interface
            tInterface = std::make_shared< Interface_User_Defined >(
                    &initialize_rosenbrock,
                    &get_criteria_rosenbrock,
                    nullptr );

            // Create Problem
            tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &get_constraint_types_rosenbrock,
                    &compute_objectives_rosenbrock,
                    &compute_constraints_rosenbrock,
                    nullptr,
                    nullptr,
                    nullptr,
                    nullptr );

            // Create manager
            tManager = Manager( tAlgorithms, tProblem );

            // Solve optimization problem
            tManager.perform();
        }

        // ---------------------------------------------------------------------------------------------------------

        SECTION( "Interface" )
        {
            if ( par_size() == 4 or par_size() == 8 )
            {
                // Create interfaces
                Vector< std::shared_ptr< Criteria_Interface > > tInterfaces( 4 );
                tInterfaces( 0 ) = std::make_shared< Interface_User_Defined >(
                        &initialize_test_1,
                        &get_criteria_test,
                        &get_dcriteria_dadv_test );
                tInterfaces( 1 ) = std::make_shared< Interface_User_Defined >(
                        &initialize_test_2,
                        &get_criteria_test,
                        &get_dcriteria_dadv_test );
                tInterfaces( 2 ) = std::make_shared< Interface_User_Defined >(
                        &initialize_test_3,
                        &get_criteria_test,
                        &get_dcriteria_dadv_test );
                tInterfaces( 3 ) = std::make_shared< Interface_User_Defined >(
                        &initialize_test_4,
                        &get_criteria_test,
                        &get_dcriteria_dadv_test );

                // Interface manager parameter list
                Parameter_List tInterfaceManagerParameterList = moris::prm::create_opt_interface_manager_parameter_list();

                // Set manager parameters
                tInterfaceManagerParameterList.set( "parallel", true );
                tInterfaceManagerParameterList.set( "shared_advs", false );
                if ( par_size() == 4 )
                {
                    tInterfaceManagerParameterList.set( "num_processors_per_interface", "1, 1, 1, 1" );
                }
                else
                {
                    tInterfaceManagerParameterList.set( "num_processors_per_interface", "1, 2, 3, 2" );
                }

                // Create manager
                Submodule_Parameter_Lists tInterfaceManager( "Interface" );
                tInterfaceManager.add_parameter_list( tInterfaceManagerParameterList );
                std::shared_ptr< Criteria_Interface > tInterface = create_interface(
                        tInterfaceManager,
                        tInterfaces );

                // Test manager in parallel
                Vector< real >  tADVs;
                Vector< real >  tLowerBounds;
                Vector< real >  tUpperBounds;
                Matrix< IdMat > tDummy;
                tInterface->initialize( tADVs, tLowerBounds, tUpperBounds, tDummy );
                for ( uint tADVIndex = 0; tADVIndex < 8; tADVIndex++ )
                {
                    REQUIRE( tADVs( tADVIndex ) == tADVIndex + 1 );
                }
                tInterface->get_criteria( tADVs );
                Matrix< DDRMat > tCriteriaGradients = tInterface->get_dcriteria_dadv();
                for ( uint tADVIndex = 0; tADVIndex < 8; tADVIndex++ )
                {
                    REQUIRE( tCriteriaGradients( tADVIndex, tADVIndex ) == tADVIndex + 1 );
                }
            }
        }

        // ---------------------------------------------------------------------------------------------------------
    }
}    // namespace moris::opt
