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
#include <filesystem>
#include <limits>
#include <cmath>

#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_OPT_Manager.hpp"
#include "fn_OPT_create_interface.hpp"
#include "cl_OPT_Problem_User_Defined.hpp"
#include "cl_OPT_Interface_User_Defined.hpp"
#include "cl_OPT_Algorithm_MMA.hpp"
#include "cl_Communication_Tools.hpp"
#include "paths.hpp"

#include "fn_OPT_Rosenbrock.hpp"
#include "cl_OPT_Algorithm_MMA.hpp"
#include "cl_OPT_Algorithm_Sweep.hpp"
#include "cl_OPT_Algorithm_TD_Slerp.hpp"
#include "fn_OPT_create_algorithm.hpp"
#include "cl_Tracer.hpp"
#include "fn_OPT_Test_Interface.hpp"
#include "HDF5_Tools.hpp"

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

            // restart files are written per iteration (MB-016); remove a stale one first
            std::remove( "ADV_Alg_0_Iter_1.hdf5" );

            // Solve optimization problem
            tManager.perform();

            // Check Solution
            if ( par_rank() == 0 )
            {
                REQUIRE( std::ifstream( "ADV_Alg_0_Iter_1.hdf5" ).good() );    // MB-016: native GCMMA writes restart files
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

        // Optional user-defined callbacks (Deck API v2, Stage 1): when a deck omits them,
        // get_constraint_types is taken from the 'constraint_types' parameter and the
        // explicit ADV gradients default to zeros — reproducing the hand-written
        // zero-returning bimodal callbacks exactly.
        SECTION( "Optional OPT callbacks default correctly" )
        {
            Parameter_List tProblemParameterList = moris::prm::create_opt_problem_parameter_list();
            tProblemParameterList.set( "constraint_types", "1" );

            std::shared_ptr< Criteria_Interface > tInterface = std::make_shared< Interface_User_Defined >(
                    &bimodal::initialize,
                    &bimodal::get_criteria,
                    &bimodal::get_dcriteria_dadv );

            // nullptr for get_constraint_types and both explicit ADV gradient callbacks
            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    nullptr,
                    &bimodal::compute_objectives,
                    &bimodal::compute_constraints,
                    nullptr,
                    &bimodal::compute_dobjective_dcriteria,
                    nullptr,
                    &bimodal::compute_dconstraint_dcriteria );

            tProblem->initialize();

            // constraint types come from the parameter
            Matrix< DDSMat > tConstraintTypes = tProblem->get_constraint_types();
            REQUIRE( tConstraintTypes.numel() == 1 );
            CHECK( tConstraintTypes( 0 ) == 1 );

            Vector< real > tADVs = tProblem->get_advs();
            tProblem->compute_design_criteria( tADVs );
            tProblem->compute_design_criteria_gradients( tADVs );

            if ( par_rank() == 0 )
            {
                // objective gradient = 0 (defaulted dObj/dADV) + dObj/dCriteria * dCriteria/dADV,
                // identical to the explicit bimodal callbacks
                const Matrix< DDRMat >& tObjGrad = tProblem->get_objective_gradients();
                REQUIRE( tObjGrad.numel() == 5 );
                CHECK( tObjGrad( 0, 0 ) == Approx( 1.0e8 ) );
                CHECK( tObjGrad( 0, 1 ) == Approx( 0.01 ) );
                CHECK( tObjGrad( 0, 4 ) == 0.0 );
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

    //------------------------------------------------------------------------------------------------------------------
    // Warm restart across an in-place redistance (MB-013): Rosenbrock in (x0, x1) plus two far-field variables
    // (x2, x3) that no criterion depends on. A test interface overwrites x2, x3 with bound-scale jumps on every
    // evaluation, the way the workflow's in-place SDF redistance rewrites far-field level-set coefficients.

    void
    initialize_rosenbrock_farfield( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds )
    {
        if ( par_rank() == 0 )
        {
            aADVs        = { 0.8, 1.2, 0.5, -0.5 };
            aLowerBounds = { -2.0, -2.0, -2.0, -2.0 };
            aUpperBounds = { 2.0, 2.0, 2.0, 2.0 };
        }
    }

    Matrix< DDRMat >
    get_dcriteria_dadv_rosenbrock_farfield( const Vector< real >& aADVs )
    {
        if ( par_rank() == 0 )
        {
            Matrix< DDRMat > tDCriteria( 2, 4, 0.0 );
            tDCriteria( 0, 0 ) = -1;
            tDCriteria( 1, 0 ) = -2 * aADVs( 0 );
            tDCriteria( 1, 1 ) = 1;
            return tDCriteria;
        }
        return { {} };
    }

    Matrix< DDRMat >
    compute_dobjective_dadv_rosenbrock_farfield( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tD2 = compute_dobjective_dadv_rosenbrock( aADVs, aCriteria );
        Matrix< DDRMat > tD4( 1, 4, 0.0 );
        tD4( 0 ) = tD2( 0 );
        tD4( 1 ) = tD2( 1 );
        return tD4;
    }

    Matrix< DDRMat >
    compute_dconstraint_dadv_rosenbrock_farfield( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tD2 = compute_dconstraint_dadv_rosenbrock( aADVs, aCriteria );
        Matrix< DDRMat > tD4( 2, 4, 0.0 );
        for ( uint iRow = 0; iRow < 2; iRow++ )
        {
            tD4( iRow, 0 ) = tD2( iRow, 0 );
            tD4( iRow, 1 ) = tD2( iRow, 1 );
        }
        return tD4;
    }

    /** User-defined interface that rewrites the far-field variables after every criteria evaluation. */
    class Interface_Farfield_Jump : public Interface_User_Defined
    {
      public:
        bool mJump  = false;
        uint mCalls = 0;

        Interface_Farfield_Jump()
                : Interface_User_Defined( &initialize_rosenbrock_farfield, &get_criteria_rosenbrock, &get_dcriteria_dadv_rosenbrock_farfield )
        {
        }

        Vector< real >
        perform( Vector< real >& aNewADVs ) override
        {
            Vector< real > tCriteria = Interface_User_Defined::perform( aNewADVs );
            if ( mJump && aNewADVs.size() == 4 )
            {
                mCalls++;
                aNewADVs( 2 ) = ( mCalls % 2 == 0 ) ? 1.8 : -1.8;    // far-field "redistance": large, alternating
                aNewADVs( 3 ) = ( mCalls % 3 == 0 ) ? 1.5 : -1.5;
            }
            return tCriteria;
        }
    };

    struct Farfield_Result
    {
        real mX0       = 0.0;
        real mX1       = 0.0;
        real mObjective = 0.0;
        int  mLocked   = 0;
    };

    Farfield_Result
    solve_rosenbrock_farfield( bool aJump, bool aShift )
    {
        Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
        Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();
        tAlgorithmParameterList.set( "norm_drop", 2.5e-5 );
        tAlgorithmParameterList.set( "asymp_adaptc", 1.1 );
        tAlgorithmParameterList.set( "reinit_in_place", true );

        std::shared_ptr< Interface_Farfield_Jump > tInterface = std::make_shared< Interface_Farfield_Jump >();
        tInterface->mJump = aJump;

        std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                tProblemParameterList,
                tInterface,
                &get_constraint_types_rosenbrock,
                &compute_objectives_rosenbrock,
                &compute_constraints_rosenbrock,
                &compute_dobjective_dadv_rosenbrock_farfield,
                &compute_dobjective_dcriteria_rosenbrock,
                &compute_dconstraint_dadv_rosenbrock_farfield,
                &compute_dconstraint_dcriteria_rosenbrock );
        tProblem->initialize();

        Algorithm_MMA tAlgorithm( tAlgorithmParameterList );
        tAlgorithm.set_shift_history( aShift );
        {
            Tracer tTracer( "OPT", "Manager", "Perform" );    // solve() books its iteration on the Manager's trace, as in Manager::perform()
            tAlgorithm.solve( 0, tProblem );
        }

        Farfield_Result tResult;
        if ( par_rank() == 0 )
        {
            Vector< real > tADVs = tProblem->get_advs();
            tResult.mX0        = tADVs( 0 );
            tResult.mX1        = tADVs( 1 );
            tResult.mObjective = tProblem->get_objectives()( 0 );
            tResult.mLocked    = tAlgorithm.get_num_locked_asymptotes();
        }
        return tResult;
    }

    TEST_CASE( "Native GCMMA warm restart across an in-place redistance", "[optimization][MB-013]" )
    {
        Farfield_Result tReference = solve_rosenbrock_farfield( false, true );    // no redistance
        Farfield_Result tShifted   = solve_rosenbrock_farfield( true, true );     // redistance + history shift
        Farfield_Result tUnshifted = solve_rosenbrock_farfield( true, false );    // redistance, old behaviour

        if ( par_rank() == 0 )
        {
            // the optimum is reached either way (the far-field variables carry no sensitivity)
            REQUIRE( std::abs( tReference.mX0 - 1.0 ) < 1E-4 );
            REQUIRE( std::abs( tShifted.mX0 - 1.0 ) < 1E-4 );

            // with the shift the redistance is invisible: same trajectory, no extra locked asymptotes
            REQUIRE( std::abs( tShifted.mX0 - tReference.mX0 ) < 1E-8 );
            REQUIRE( std::abs( tShifted.mX1 - tReference.mX1 ) < 1E-8 );
            REQUIRE( std::abs( tShifted.mObjective - tReference.mObjective ) < 1E-10 );
            REQUIRE( tShifted.mLocked == tReference.mLocked );

            // negative control: without the shift the jumps collapse the far-field asymptotes onto the floor
            REQUIRE( tUnshifted.mLocked > tReference.mLocked );
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    // In-place redistance checkpoints (MB-021): with reinit_in_place, the restart file of every criteria evaluation
    // must hold the design the criteria were evaluated at (the interface's redistanced output), not the raw trial
    // point. The evaluation advances the optimization iteration, so the rewrite has to target the index captured
    // before it.

    /** Far-field-jump interface that records the redistanced design returned by every criteria evaluation. */
    class Interface_Recording_Redistance : public Interface_Farfield_Jump
    {
      public:
        std::vector< std::vector< real > > mRedistanced;

        Vector< real >
        perform( Vector< real >& aNewADVs ) override
        {
            Vector< real > tCriteria = Interface_Farfield_Jump::perform( aNewADVs );

            std::vector< real > tCopy( aNewADVs.size() );
            for ( uint iADV = 0; iADV < aNewADVs.size(); iADV++ )
            {
                tCopy[ iADV ] = aNewADVs( iADV );
            }
            mRedistanced.push_back( tCopy );

            return tCriteria;
        }
    };

    std::string
    mb021_restart_file_name( uint aIter )
    {
        return "ADV_Alg_0_Iter_" + std::to_string( aIter ) + ".hdf5";
    }

    std::vector< real >
    mb021_read_restart_advs( uint aIter )
    {
        hid_t  tFileID = open_hdf5_file( mb021_restart_file_name( aIter ), false, true );
        herr_t tStatus = 0;

        std::vector< real > tADVs;
        load_vector_from_hdf5_file( tFileID, "ADVs", tADVs, tStatus );
        close_hdf5_file( tFileID );

        return tADVs;
    }

    TEST_CASE( "Native GCMMA in-place redistance checkpoints hold the redistanced design", "[optimization][MB-021]" )
    {
        Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
        Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();
        tAlgorithmParameterList.set( "max_its", 3 );
        tAlgorithmParameterList.set( "reinit_in_place", true );

        std::shared_ptr< Interface_Recording_Redistance > tInterface = std::make_shared< Interface_Recording_Redistance >();

        // redistance on every evaluation
        tInterface->mJump = true;

        std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                tProblemParameterList,
                tInterface,
                &get_constraint_types_rosenbrock,
                &compute_objectives_rosenbrock,
                &compute_constraints_rosenbrock,
                &compute_dobjective_dadv_rosenbrock_farfield,
                &compute_dobjective_dcriteria_rosenbrock,
                &compute_dconstraint_dadv_rosenbrock_farfield,
                &compute_dconstraint_dcriteria_rosenbrock );
        tProblem->initialize();

        // remove stale restart files of earlier tests in this directory
        for ( uint iIter = 0; iIter < 50; iIter++ )
        {
            std::remove( mb021_restart_file_name( iIter ).c_str() );
        }

        Algorithm_MMA tAlgorithm( tAlgorithmParameterList );
        {
            Tracer tTracer( "OPT", "Manager", "Perform" );
            tAlgorithm.solve( 0, tProblem );
        }

        if ( par_rank() == 0 )
        {
            // the initial design plus at least one evaluation per outer iteration
            const uint tNumEvaluations = tInterface->mRedistanced.size();
            REQUIRE( tNumEvaluations >= 4 );

            // evaluation k is checkpointed in Iter_k, and the file holds its redistanced design, the last one included
            for ( uint iEval = 0; iEval < tNumEvaluations; iEval++ )
            {
                REQUIRE( std::ifstream( mb021_restart_file_name( iEval ) ).good() );

                std::vector< real > tFileADVs = mb021_read_restart_advs( iEval );
                REQUIRE( tFileADVs.size() == 4 );
                for ( uint iADV = 0; iADV < 4; iADV++ )
                {
                    CHECK( tFileADVs[ iADV ] == tInterface->mRedistanced[ iEval ][ iADV ] );
                }
            }

            // no file past the last evaluation (the rewrite used to land in Iter_{k+1})
            REQUIRE_FALSE( std::ifstream( mb021_restart_file_name( tNumEvaluations ) ).good() );
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    // Redistance on accepted iterates only (reinit_on_accept_only): with the GCMMA inner loop on, every conservatism trial
    // must be evaluated at the raw subproblem solution psi (the design the approximation was built for), and only the
    // accepted iterate is redistanced, by one extra evaluation at that same psi. The interface below redistances (jumps
    // the far-field variables) only when the optimizer has not suppressed it, and records every call.

    /** Far-field-jump interface that honours the optimizer's redistance suppression and records every call. */
    class Interface_Suppressible_Redistance : public Interface_Farfield_Jump
    {
      public:
        std::vector< bool >                mSuppressed;    ///< per call: was the redistance suppressed
        std::vector< std::vector< real > > mInput;         ///< per call: the design handed in (before any redistance)

        Vector< real >
        perform( Vector< real >& aNewADVs ) override
        {
            std::vector< real > tIn( aNewADVs.size() );
            for ( uint iADV = 0; iADV < aNewADVs.size(); iADV++ )
            {
                tIn[ iADV ] = aNewADVs( iADV );
            }
            mInput.push_back( tIn );
            mSuppressed.push_back( this->is_reinit_suppressed() );

            mJump = !this->is_reinit_suppressed();
            return Interface_Farfield_Jump::perform( aNewADVs );
        }
    };

    std::shared_ptr< Interface_Suppressible_Redistance >
    solve_suppressible_redistance( bool aAcceptOnly )
    {
        Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
        Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();
        tAlgorithmParameterList.set( "max_its", 12 );
        tAlgorithmParameterList.set( "max_inner_its", 5 );
        tAlgorithmParameterList.set( "step_size", 0.3 );
        tAlgorithmParameterList.set( "reinit_in_place", true );
        tAlgorithmParameterList.set( "reinit_on_accept_only", aAcceptOnly );

        std::shared_ptr< Interface_Suppressible_Redistance > tInterface = std::make_shared< Interface_Suppressible_Redistance >();

        std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                tProblemParameterList,
                tInterface,
                &get_constraint_types_rosenbrock,
                &compute_objectives_rosenbrock,
                &compute_constraints_rosenbrock,
                &compute_dobjective_dadv_rosenbrock_farfield,
                &compute_dobjective_dcriteria_rosenbrock,
                &compute_dconstraint_dadv_rosenbrock_farfield,
                &compute_dconstraint_dcriteria_rosenbrock );
        tProblem->initialize();

        Algorithm_MMA tAlgorithm( tAlgorithmParameterList );
        {
            Tracer tTracer( "OPT", "Manager", "Perform" );
            tAlgorithm.solve( 0, tProblem );
        }

        // the suppression flag is never left set after the solve
        REQUIRE_FALSE( tInterface->is_reinit_suppressed() );

        return tInterface;
    }

    TEST_CASE( "Native GCMMA redistances only accepted iterates with reinit_on_accept_only", "[optimization][reinit_on_accept_only]" )
    {
        std::shared_ptr< Interface_Suppressible_Redistance > tAcceptOnly = solve_suppressible_redistance( true );
        std::shared_ptr< Interface_Suppressible_Redistance > tEveryEval  = solve_suppressible_redistance( false );

        if ( par_rank() == 0 )
        {
            const std::vector< bool >& tS = tAcceptOnly->mSuppressed;
            REQUIRE( tS.size() >= 3 );

            // the initial evaluation is an accepted iterate: redistanced
            REQUIRE_FALSE( tS[ 0 ] );

            uint tNumTrials   = 0;
            uint tNumAccepted = 0;
            for ( uint iCall = 1; iCall < tS.size(); iCall++ )
            {
                if ( tS[ iCall ] )
                {
                    tNumTrials++;
                    continue;
                }

                // a redistanced evaluation follows a raw trial and re-evaluates exactly that (accepted) design
                tNumAccepted++;
                REQUIRE( tS[ iCall - 1 ] );
                for ( uint iADV = 0; iADV < tAcceptOnly->mInput[ iCall ].size(); iADV++ )
                {
                    CHECK( tAcceptOnly->mInput[ iCall ][ iADV ] == tAcceptOnly->mInput[ iCall - 1 ][ iADV ] );
                }
            }

            // every outer iteration ended with its accept evaluation, and the inner loop ran: more trials than iterations
            REQUIRE( !tS.back() );
            REQUIRE( tNumAccepted >= 1 );
            REQUIRE( tNumTrials > tNumAccepted );

            // negative control: without the key, no evaluation is suppressed (every trial is redistanced)
            for ( uint iCall = 0; iCall < tEveryEval->mSuppressed.size(); iCall++ )
            {
                REQUIRE_FALSE( tEveryEval->mSuppressed[ iCall ] );
            }
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    // SDF taper hooks: the optimizer reports every accepted iterate (previous committed design, accepted subproblem
    // solution) to the interface, and in accept-only mode stops spending the accept evaluation once the interface
    // reports the redistance inactive.

    /** Suppressible far-field interface whose redistance tapers off after mActiveFor accepted iterates. */
    class Interface_Tapering_Redistance : public Interface_Suppressible_Redistance
    {
      public:
        uint                               mActiveFor = 3;
        std::vector< std::vector< real > > mNotifiedPrevious;
        std::vector< std::vector< real > > mNotifiedAccepted;

        void
        notify_accepted_iterate( const Vector< real >& aPrevious, const Vector< real >& aAccepted ) override
        {
            // an empty aAccepted marks an iterate accepted at the inner-loop cap
            std::vector< real > tPrev( aPrevious.size() ), tAcc( aAccepted.size() );
            for ( uint iADV = 0; iADV < aPrevious.size(); iADV++ )
            {
                tPrev[ iADV ] = aPrevious( iADV );
            }
            for ( uint iADV = 0; iADV < aAccepted.size(); iADV++ )
            {
                tAcc[ iADV ] = aAccepted( iADV );
            }
            mNotifiedPrevious.push_back( tPrev );
            mNotifiedAccepted.push_back( tAcc );
        }

        bool
        is_reinit_active() const override
        {
            return mNotifiedAccepted.size() <= mActiveFor;    // active for the first mActiveFor accepted iterates
        }
    };

    TEST_CASE( "Native GCMMA reports accepted iterates and skips the accept redistance once tapered off", "[optimization][sdf_taper]" )
    {
        Parameter_List tProblemParameterList   = moris::prm::create_opt_problem_parameter_list();
        Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();
        tAlgorithmParameterList.set( "max_its", 8 );
        tAlgorithmParameterList.set( "max_inner_its", 5 );
        tAlgorithmParameterList.set( "step_size", 0.3 );
        tAlgorithmParameterList.set( "reinit_in_place", true );
        tAlgorithmParameterList.set( "reinit_on_accept_only", true );

        std::shared_ptr< Interface_Tapering_Redistance > tInterface = std::make_shared< Interface_Tapering_Redistance >();

        std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                tProblemParameterList,
                tInterface,
                &get_constraint_types_rosenbrock,
                &compute_objectives_rosenbrock,
                &compute_constraints_rosenbrock,
                &compute_dobjective_dadv_rosenbrock_farfield,
                &compute_dobjective_dcriteria_rosenbrock,
                &compute_dconstraint_dadv_rosenbrock_farfield,
                &compute_dconstraint_dcriteria_rosenbrock );
        tProblem->initialize();

        Algorithm_MMA tAlgorithm( tAlgorithmParameterList );
        {
            Tracer tTracer( "OPT", "Manager", "Perform" );
            tAlgorithm.solve( 0, tProblem );
        }

        if ( par_rank() == 0 )
        {
            // one notification per accepted outer iteration
            REQUIRE( tInterface->mNotifiedAccepted.size() == 8 );

            // the accepted iterate is the raw trial that preceded each accept evaluation, and the previous design is the
            // committed (redistanced) one, i.e. differs from the raw accepted point in the far field while active
            uint tRedistanced = 0;
            for ( uint iCall = 1; iCall < tInterface->mSuppressed.size(); iCall++ )
            {
                tRedistanced += tInterface->mSuppressed[ iCall ] ? 0 : 1;
            }

            // the accept evaluation runs only while the interface reports the redistance active: 3 accepted iterates
            CHECK( tRedistanced == tInterface->mActiveFor );

            // after the taper every evaluation is a raw trial
            CHECK( tInterface->mSuppressed.back() );
        }
    }

    //------------------------------------------------------------------------------------------------------------------
    // Failed evaluations (MB-023, MB-024) and SWEEP finite-difference modes (MB-025), on the one-variable quadratic
    // f(x) = 100 (x - 0.45)^2 started at x = 0.5. A failed evaluation is modelled as a NaN criterion, the signature
    // of an XTK decomposition failure.

    namespace failed_evaluation
    {
        enum class Fault
        {
            none,                  ///< every evaluation is finite
            single_evaluation,     ///< only evaluation number mFaultyEvaluation fails
            off_initial_design     ///< every design other than the initial x = 0.5 fails
        };

        Fault mFault             = Fault::none;
        int   mFaultyEvaluation  = 0;
        int   mEvaluations       = 0;
        int   mGradientsOnFailed = 0;
        bool  mLastFailed        = false;

        void
        reset( Fault aFault, int aFaultyEvaluation = 0 )
        {
            mFault             = aFault;
            mFaultyEvaluation  = aFaultyEvaluation;
            mEvaluations       = 0;
            mGradientsOnFailed = 0;
            mLastFailed        = false;
        }

        void
        initialize(
                Vector< real >& aADVs,
                Vector< real >& aLowerBounds,
                Vector< real >& aUpperBounds )
        {
            aADVs        = { 0.5 };
            aLowerBounds = { 0.0 };
            aUpperBounds = { 1.0 };
        }

        Vector< real >
        get_criteria( const Vector< real >& aADVs )
        {
            mEvaluations++;
            mLastFailed = ( mFault == Fault::single_evaluation and mEvaluations == mFaultyEvaluation )
                       or ( mFault == Fault::off_initial_design and aADVs( 0 ) != 0.5 );

            if ( mLastFailed )
            {
                return { std::numeric_limits< real >::quiet_NaN() };
            }
            return { 100.0 * std::pow( aADVs( 0 ) - 0.45, 2 ) };
        }

        Matrix< DDRMat >
        get_dcriteria_dadv( const Vector< real >& aADVs )
        {
            if ( mLastFailed )
            {
                mGradientsOnFailed++;
            }
            return { { 200.0 * ( aADVs( 0 ) - 0.45 ) } };
        }

        Matrix< DDRMat >
        zero_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return { { 0.0 } };
        }

        /** remove the per-evaluation restart files the native MMA writes into the working directory */
        void
        remove_restart_files()
        {
            // the iteration counter is global to the process, so the indices depend on the tests run before
            for ( const std::filesystem::directory_entry& tEntry : std::filesystem::directory_iterator( "." ) )
            {
                const std::string tName = tEntry.path().filename().string();
                if ( tName.rfind( "ADV_Alg_0_Iter_", 0 ) == 0 )
                {
                    std::filesystem::remove( tEntry.path() );
                }
            }
        }

        std::shared_ptr< Problem >
        create_problem()
        {
            Parameter_List tProblemParameterList = moris::prm::create_opt_problem_parameter_list();

            std::shared_ptr< Criteria_Interface > tInterface =
                    std::make_shared< Interface_User_Defined >( &initialize, &get_criteria, &get_dcriteria_dadv );

            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &bimodal::get_constraint_types,
                    &bimodal::compute_objectives,
                    &bimodal::compute_constraints,
                    &zero_dadv,
                    &bimodal::compute_dobjective_dcriteria,
                    &zero_dadv,
                    &bimodal::compute_dconstraint_dcriteria );
            tProblem->initialize();

            return tProblem;
        }
    }    // namespace failed_evaluation

    TEST_CASE( "Native GCMMA rejects a failed conservatism inner trial", "[optimization][MB-023]" )
    {
        if ( par_rank() == 0 )
        {
            // evaluation 1 is the initial design, 2 the (non-conservative) first trial, 3 the first inner trial
            failed_evaluation::reset( failed_evaluation::Fault::single_evaluation, 3 );
            failed_evaluation::remove_restart_files();
            std::shared_ptr< Problem > tProblem = failed_evaluation::create_problem();

            Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();
            tAlgorithmParameterList.set( "max_its", 1 );
            tAlgorithmParameterList.set( "max_inner_its", 3 );
            tAlgorithmParameterList.set( "step_size", 0.3 );

            Algorithm_MMA tAlgorithm( tAlgorithmParameterList );
            {
                Tracer tTracer( "OPT", "Manager", "Perform" );
                tAlgorithm.solve( 0, tProblem );
            }

            // the inner loop was entered and its failed trial was evaluated
            REQUIRE( failed_evaluation::mEvaluations > 3 );

            // ... but never differentiated or committed: the backtrack re-solved to a finite design
            CHECK( failed_evaluation::mGradientsOnFailed == 0 );
            const real tX = tProblem->get_advs()( 0 );
            CHECK( std::isfinite( tProblem->get_objectives()( 0 ) ) );
            CHECK( tProblem->get_objectives()( 0 ) == Approx( 100.0 * std::pow( tX - 0.45, 2 ) ) );

            failed_evaluation::remove_restart_files();
        }
    }

    TEST_CASE( "Native GCMMA restores the accepted state after exhausted backtracking", "[optimization][MB-024]" )
    {
        if ( par_rank() == 0 )
        {
            // every trial fails, however short the step; only the initial design evaluates
            failed_evaluation::reset( failed_evaluation::Fault::off_initial_design );
            failed_evaluation::remove_restart_files();
            std::shared_ptr< Problem > tProblem = failed_evaluation::create_problem();

            Parameter_List tAlgorithmParameterList = moris::prm::create_mma_parameter_list();
            tAlgorithmParameterList.set( "max_its", 1 );
            tAlgorithmParameterList.set( "max_inner_its", 0 );
            tAlgorithmParameterList.set( "step_size", 0.3 );

            Algorithm_MMA tAlgorithm( tAlgorithmParameterList );
            {
                Tracer tTracer( "OPT", "Manager", "Perform" );
                tAlgorithm.solve( 0, tProblem );
            }

            // initial design + first trial + 15 backtracks + the restore evaluation of the held design
            CHECK( failed_evaluation::mEvaluations == 18 );
            CHECK( failed_evaluation::mGradientsOnFailed == 0 );

            // the Problem the Manager reads holds the accepted design and its criteria, not the last rejected trial
            CHECK( tProblem->get_advs()( 0 ) == 0.5 );
            CHECK( tProblem->get_objectives()( 0 ) == Approx( 0.25 ) );

            failed_evaluation::remove_restart_files();
        }
    }

    TEST_CASE( "SWEEP evaluates every requested finite difference scheme", "[optimization][MB-025][MB-026]" )
    {
        // f(x) = 100 (x - 0.45)^2 is quadratic, so its finite differences are exact up to round-off:
        // forward f' + 100 eps, backward f' - 100 eps, central f'. Two evaluation points and two epsilons make every
        // scheme run repeatedly, which exposes a perturbation sign that drifts between calls (MB-026).
        const std::vector< real > tPoints   = { 0.5, 0.6 };
        const std::vector< real > tEpsilons = { 1e-3, 1e-4 };

        // criteria evaluations per point with one ADV: the base point plus the perturbed points of every epsilon
        std::vector< std::pair< std::string, int > > tModes = { { "none", 1 }, { "forward", 3 }, { "backward", 3 }, { "central", 5 }, { "all", 5 } };

        for ( const std::pair< std::string, int >& tMode : tModes )
        {
            failed_evaluation::reset( failed_evaluation::Fault::none );
            std::shared_ptr< Problem > tProblem = failed_evaluation::create_problem();

            const std::string tFileName = "mb025_sweep_" + tMode.first + ".hdf5";

            Parameter_List tAlgorithmParameterList = moris::prm::create_sweep_parameter_list();
            tAlgorithmParameterList.set( "custom_adv_evaluations", "0.5, 0.6" );
            tAlgorithmParameterList.set( "finite_difference_type", tMode.first );
            tAlgorithmParameterList.set( "finite_difference_epsilons", "1e-3, 1e-4" );
            tAlgorithmParameterList.set( "hdf5_path", tFileName );
            tAlgorithmParameterList.set( "evaluate_constraint_gradients", false );

            Algorithm_Sweep tAlgorithm( tAlgorithmParameterList );
            {
                Tracer tTracer( "OPT", "Manager", "Perform" );
                tAlgorithm.solve( 0, tProblem );
            }

            if ( par_rank() == 0 )
            {
                INFO( "mode " << tMode.first );
                CHECK( failed_evaluation::mEvaluations == 2 * tMode.second );

                // datasets expected: the analytic reference for none/all, and every requested scheme
                std::vector< std::string > tSchemes;
                if ( tMode.first == "all" )
                {
                    tSchemes = { "fd_forward", "fd_backward", "fd_central" };
                }
                else if ( tMode.first != "none" )
                {
                    tSchemes = { "fd_" + tMode.first };
                }
                const bool tAnalytical = tMode.first == "none" or tMode.first == "all";

                hid_t  tFileID = open_hdf5_file( tFileName, true, true );
                herr_t tStatus = 0;
                for ( uint iPoint = 0; iPoint < tPoints.size(); iPoint++ )
                {
                    const std::string tEval    = "objective_gradients eval_" + std::to_string( iPoint + 1 ) + "-2";
                    const real        tExact   = 200.0 * ( tPoints[ iPoint ] - 0.45 );
                    Matrix< DDRMat >  tGradient;

                    CHECK( ( H5Lexists( tFileID, ( tEval + " analytical" ).c_str(), H5P_DEFAULT ) > 0 ) == tAnalytical );

                    for ( uint iEps = 0; iEps < tEpsilons.size(); iEps++ )
                    {
                        for ( const std::string& tScheme : tSchemes )
                        {
                            const std::string tLabel = tEval + " epsilon_" + std::to_string( iEps + 1 ) + "-2 " + tScheme;
                            INFO( "dataset '" << tLabel << "'" );
                            REQUIRE( H5Lexists( tFileID, tLabel.c_str(), H5P_DEFAULT ) > 0 );

                            real tExpected = tExact;
                            if ( tScheme == "fd_forward" )
                            {
                                tExpected += 100.0 * tEpsilons[ iEps ];
                            }
                            else if ( tScheme == "fd_backward" )
                            {
                                tExpected -= 100.0 * tEpsilons[ iEps ];
                            }

                            load_matrix_from_hdf5_file( tFileID, tLabel, tGradient, tStatus );
                            CHECK( tGradient( 0 ) == Approx( tExpected ).epsilon( 0.0 ).margin( 1e-6 ) );
                        }
                    }
                }
                close_hdf5_file( tFileID );
                std::remove( tFileName.c_str() );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    // A restart file written on another rank count lists the ADVs in that layout's order (MB-077): the values follow their
    // ijkl IDs, an ADV without an ID keeps its position, and an ID missing from the current layout is an error
    TEST_CASE( "Restart values are mapped to the current layout by ijkl ID", "[optimization][MB-077]" )
    {
        Matrix< IdMat > tFileIds    = { { 30 }, { 10 }, { gNoID }, { 20 } };
        Matrix< IdMat > tCurrentIds = { { 10 }, { 20 }, { gNoID }, { 30 } };
        Vector< real >  tFileValues = { 3.0, 1.0, 9.0, 2.0 };

        Vector< real > tValues = Problem::map_restart_values_to_layout( tFileIds, tCurrentIds, tFileValues );
        REQUIRE( tValues.size() == 4 );
        CHECK( tValues( 0 ) == 1.0 );
        CHECK( tValues( 1 ) == 2.0 );
        CHECK( tValues( 2 ) == 9.0 );
        CHECK( tValues( 3 ) == 3.0 );

        // the same layout maps to itself
        Vector< real > tSame = Problem::map_restart_values_to_layout( tCurrentIds, tCurrentIds, tValues );
        CHECK( tSame( 0 ) == 1.0 );
        CHECK( tSame( 3 ) == 3.0 );

        Matrix< IdMat > tOtherIds = { { 10 }, { 20 }, { gNoID }, { 40 } };
        CHECK_THROWS( Problem::map_restart_values_to_layout( tFileIds, tOtherIds, tFileValues ) );
    }
    //------------------------------------------------------------------------------------------------------------------
    // TD_SLERP (sandbox-topopt IMPL_PLAN class C, U12-U17): topological-derivative level-set evolution on the unit sphere
    // of the M inner product, driven through a mock criteria interface that supplies J, the design direction g~ and M.

    namespace td_slerp
    {
        enum class Mode
        {
            OVERSHOOT,      // J = -<phi, a>_M / |phi|_M; g~ at angle -1.5 alpha from a (alpha = angle(phi, a)): k = 1 overshoots
            FIXED_POINT,    // g~ = 2.5 phi (stationary point g~ = s phi, s > 0)
            ANTIPODE,       // g~ = -2.5 phi
            RISING          // J = |phi / |phi|_M - phi_0 / |phi_0|_M|_M: every trial raises J; g~ = a
        };

        class Mock_Interface : public Criteria_Interface
        {
          public:
            Matrix< DDRMat > mM;        // SPD mass matrix
            Vector< real >   mStart;    // initial ADVs
            Vector< real >   mA;        // target direction
            Mode             mMode;
            real             mSign = 1.0;    // sign flip: J_s( phi ) = J( s phi ), g_s( phi ) = s g( s phi )

            std::vector< std::vector< real > > mEvaluated;     // every evaluated design (as submitted)
            std::vector< real >                mObjectives;    // J of every evaluation
            std::vector< uint >                mFetchAfter;    // number of evaluations done at every direction fetch
            Vector< real >                     mLast;          // s phi of the last evaluation

            Mock_Interface( const Matrix< DDRMat >& aM, const Vector< real >& aStart, const Vector< real >& aA, Mode aMode, real aSign = 1.0 )
                    : mM( aM )
                    , mStart( aStart )
                    , mA( aA )
                    , mMode( aMode )
                    , mSign( aSign )
            {
            }

            real
            dot( const Vector< real >& aFirst, const Vector< real >& aSecond ) const
            {
                real tSum = 0.0;
                for ( uint i = 0; i < aFirst.size(); i++ )
                {
                    for ( uint j = 0; j < aSecond.size(); j++ )
                    {
                        tSum += aFirst( i ) * mM( i, j ) * aSecond( j );
                    }
                }
                return tSum;
            }

            real
            mnorm( const Vector< real >& aVector ) const
            {
                return std::sqrt( this->dot( aVector, aVector ) );
            }

            void
            initialize( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds, Matrix< IdMat >& ) override
            {
                if ( par_rank() == 0 )    // the optimizer's rank holds the design
                {
                    aADVs        = mStart;
                    aLowerBounds = Vector< real >( mStart.size(), -10.0 );
                    aUpperBounds = Vector< real >( mStart.size(), 10.0 );
                }
            }

            real
            objective( const Vector< real >& aPhi ) const
            {
                if ( mMode == Mode::RISING )
                {
                    Vector< real > tDiff( aPhi.size() );
                    real           tNorm  = this->mnorm( aPhi );
                    real           tNorm0 = this->mnorm( mStart );
                    for ( uint i = 0; i < aPhi.size(); i++ )
                    {
                        tDiff( i ) = aPhi( i ) / tNorm - mSign * mStart( i ) / tNorm0;    // mStart is s phi_0
                    }
                    return this->mnorm( tDiff );
                }
                return -this->dot( aPhi, mA ) / ( this->mnorm( aPhi ) * this->mnorm( mA ) );
            }

            Vector< real >
            direction( const Vector< real >& aPhi ) const
            {
                uint           tN = aPhi.size();
                Vector< real > tG( tN, 0.0 );
                if ( mMode == Mode::FIXED_POINT || mMode == Mode::ANTIPODE )
                {
                    real tS = mMode == Mode::FIXED_POINT ? 2.5 : -2.5;
                    for ( uint i = 0; i < tN; i++ ) tG( i ) = tS * aPhi( i );
                    return tG;
                }
                real tNormA = this->mnorm( mA );
                if ( mMode == Mode::RISING )
                {
                    for ( uint i = 0; i < tN; i++ ) tG( i ) = mA( i ) / tNormA;
                    return tG;
                }
                // OVERSHOOT: in the plane of a and phi, at angle -1.5 alpha from a (the far side of a), |g~|_M = 3.7
                real tNormPhi = this->mnorm( aPhi );
                real tCos     = std::max( -1.0, std::min( 1.0, this->dot( aPhi, mA ) / ( tNormPhi * tNormA ) ) );
                real tAlpha   = std::acos( tCos );
                if ( std::sin( tAlpha ) < 1.0e-14 )
                {
                    for ( uint i = 0; i < tN; i++ ) tG( i ) = 3.7 * mA( i ) / tNormA;
                    return tG;
                }
                for ( uint i = 0; i < tN; i++ )
                {
                    real tU = ( aPhi( i ) / tNormPhi - tCos * mA( i ) / tNormA ) / std::sin( tAlpha );
                    tG( i ) = 3.7 * ( std::cos( 1.5 * tAlpha ) * mA( i ) / tNormA - std::sin( 1.5 * tAlpha ) * tU );
                }
                return tG;
            }

            Vector< real >
            perform( Vector< real >& aNewADVs ) override
            {
                if ( aNewADVs.size() == 0 )
                {
                    return Vector< real >( 1, 0.0 );    // ranks > 0 (dummy solve)
                }
                std::vector< real > tCopy( aNewADVs.size() );
                mLast = aNewADVs;
                for ( uint i = 0; i < aNewADVs.size(); i++ )
                {
                    tCopy[ i ] = aNewADVs( i );
                    mLast( i ) = mSign * aNewADVs( i );
                }
                mEvaluated.push_back( tCopy );
                real tJ = this->objective( mLast );
                mObjectives.push_back( tJ );
                return Vector< real >( 1, tJ );
            }

            Matrix< DDRMat >
            compute_dcriteria_dadv() override
            {
                return Matrix< DDRMat >( 1, mStart.size(), 0.0 );
            }

            Vector< real >
            get_design_direction() override
            {
                mFetchAfter.push_back( mEvaluated.size() );
                Vector< real > tG = this->direction( mLast );
                for ( uint i = 0; i < tG.size(); i++ ) tG( i ) = mSign * tG( i );
                return tG;
            }

            real
            mass_inner_product( const Vector< real >& aFirst, const Vector< real >& aSecond ) override
            {
                return this->dot( aFirst, aSecond );
            }
        };

        Matrix< DDSMat >
        constraint_types()
        {
            return Matrix< DDSMat >( 0, 1 );
        }

        Matrix< DDRMat >
        objectives( const Vector< real >&, const Vector< real >& aCriteria )
        {
            return Matrix< DDRMat >( 1, 1, aCriteria( 0 ) );
        }

        Matrix< DDRMat >
        constraints( const Vector< real >&, const Vector< real >& )
        {
            return Matrix< DDRMat >( 0, 1 );
        }

        Matrix< DDRMat >
        zeros( const Vector< real >& aADVs, const Vector< real >& )
        {
            return Matrix< DDRMat >( 1, aADVs.size(), 0.0 );
        }

        struct Result
        {
            Algorithm_TD_Slerp::Exit mExit        = Algorithm_TD_Slerp::Exit::NONE;
            uint                     mAccepted    = 0;
            uint                     mEvaluations = 0;
            std::vector< real >      mObjectives;
            std::vector< real >      mThetas;
            std::vector< real >      mKappas;
            real                     mFinalObjective = 0.0;
        };

        Result
        run( const std::shared_ptr< Mock_Interface >& aInterface, Parameter_List aAlgorithmParameterList )
        {
            Parameter_List             tProblemParameterList = moris::prm::create_opt_problem_parameter_list();
            std::shared_ptr< Problem > tProblem              = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList, aInterface, &constraint_types, &objectives, &constraints, &zeros, &zeros, &zeros, &zeros );
            tProblem->initialize();

            Algorithm_TD_Slerp tAlgorithm( aAlgorithmParameterList );
            Result             tResult;
            {
                Tracer tTracer( "OPT", "Manager", "Perform" );    // solve() books its iteration on the Manager's trace
                tAlgorithm.solve( 0, tProblem );
            }

            if ( par_rank() == 0 )
            {
                tResult.mExit           = tAlgorithm.get_exit();
                tResult.mAccepted       = tAlgorithm.get_num_accepted();
                tResult.mEvaluations    = tAlgorithm.get_num_evaluations();
                tResult.mObjectives     = tAlgorithm.get_objective_history();
                tResult.mThetas         = tAlgorithm.get_theta_history();
                tResult.mKappas         = tAlgorithm.get_kappa_history();
                tResult.mFinalObjective = tProblem->get_objectives()( 0 );
            }
            return tResult;
        }

        // the two mass matrices: identity, then a non-diagonal SPD one (1D linear B-spline mass pattern)
        Matrix< DDRMat >
        mass_matrix( bool aIdentity )
        {
            Matrix< DDRMat > tM( 4, 4, 0.0 );
            for ( uint i = 0; i < 4; i++ )
            {
                tM( i, i ) = aIdentity ? 1.0 : 4.0 / 6.0;
                if ( !aIdentity && i > 0 )
                {
                    tM( i, i - 1 ) = 1.0 / 6.0;
                    tM( i - 1, i ) = 1.0 / 6.0;
                }
            }
            return tM;
        }

        std::vector< real >
        read_restart_advs( uint aIter )
        {
            hid_t  tFileID = open_hdf5_file( "ADV_Alg_0_Iter_" + std::to_string( aIter ) + ".hdf5", false, true );
            herr_t tStatus = 0;

            std::vector< real > tADVs;
            load_vector_from_hdf5_file( tFileID, "ADVs", tADVs, tStatus );
            close_hdf5_file( tFileID );
            return tADVs;
        }

        const Vector< real > gStart = { -0.7, 0.2, 0.9, -0.3 };    // not on the unit sphere
        const Vector< real > gA     = { -0.2, 0.9, 0.6, 0.3 };     // about 55 deg from gStart (both M)
    }    // namespace td_slerp

    TEST_CASE( "TD_SLERP parameter list and factory", "[optimization][td_slerp][U12]" )
    {
        Parameter_List tParameterList = moris::prm::create_optimization_algorithm_parameter_list( Optimization_Algorithm_Type::TD_SLERP );
        CHECK( tParameterList.get< std::string >( "algorithm" ) == "td_slerp" );
        CHECK( tParameterList.get< real >( "kappa_min" ) == std::pow( 2.0, -10 ) );
        CHECK( tParameterList.get< real >( "theta_tol" ) == 1.0 );
        CHECK( tParameterList.get< std::string >( "kappa_policy" ) == "reset" );

        std::shared_ptr< Algorithm > tAlgorithm = create_algorithm( tParameterList );
        CHECK( std::dynamic_pointer_cast< Algorithm_TD_Slerp >( tAlgorithm ) != nullptr );

        tParameterList.set( "kappa_policy", "other" );
        CHECK_THROWS( create_algorithm( tParameterList ) );
    }

    TEST_CASE( "TD_SLERP keeps every evaluated design on the unit M sphere and decreases J", "[optimization][td_slerp][U13]" )
    {
        for ( bool tIdentity : { true, false } )
        {
            std::shared_ptr< td_slerp::Mock_Interface > tInterface = std::make_shared< td_slerp::Mock_Interface >(
                    td_slerp::mass_matrix( tIdentity ), td_slerp::gStart, td_slerp::gA, td_slerp::Mode::OVERSHOOT );
            td_slerp::Result tResult = td_slerp::run( tInterface, moris::prm::create_td_slerp_parameter_list() );

            if ( par_rank() == 0 )
            {
                // the overshooting direction: k = 1 is rejected, k = 1/2 accepted every iteration, then THETA
                CHECK( tResult.mExit == Algorithm_TD_Slerp::Exit::THETA );
                REQUIRE( tResult.mAccepted >= 3 );
                CHECK( tResult.mKappas.size() == tResult.mAccepted );
                for ( real tKappa : tResult.mKappas )
                {
                    CHECK( tKappa == 0.5 );
                }
                // the start as given, the normalized start, two trials per iteration
                CHECK( tResult.mEvaluations == 2 + 2 * tResult.mAccepted );
                CHECK( tResult.mThetas.back() < 1.0 );
                CHECK( tInterface->mEvaluated[ 0 ][ 0 ] == td_slerp::gStart( 0 ) );

                // unit M norm of every evaluated design after the start as given (normalize_start)
                for ( uint k = 1; k < tInterface->mEvaluated.size(); k++ )
                {
                    Vector< real > tV( 4 );
                    for ( uint i = 0; i < 4; i++ ) tV( i ) = tInterface->mEvaluated[ k ][ i ];
                    CHECK( std::abs( tInterface->mnorm( tV ) - 1.0 ) < 1.0e-14 );
                }

                // accepted J non-increasing, approaching the minimum -1 (phi -> a)
                for ( uint i = 1; i < tResult.mObjectives.size(); i++ )
                {
                    CHECK( tResult.mObjectives[ i ] <= tResult.mObjectives[ i - 1 ] );
                }
                CHECK( tResult.mObjectives.back() < -0.9999 );
            }
        }
    }

    TEST_CASE( "TD_SLERP exits THETA at a fixed point and THETA_PI at the antipode before any update", "[optimization][td_slerp][U14][U15]" )
    {
        for ( bool tIdentity : { true, false } )
        {
            for ( td_slerp::Mode tMode : { td_slerp::Mode::FIXED_POINT, td_slerp::Mode::ANTIPODE } )
            {
                std::shared_ptr< td_slerp::Mock_Interface > tInterface = std::make_shared< td_slerp::Mock_Interface >(
                        td_slerp::mass_matrix( tIdentity ), td_slerp::gStart, td_slerp::gA, tMode );
                td_slerp::Result tResult = td_slerp::run( tInterface, moris::prm::create_td_slerp_parameter_list() );

                if ( par_rank() == 0 )
                {
                    CHECK( tResult.mExit == ( tMode == td_slerp::Mode::FIXED_POINT ? Algorithm_TD_Slerp::Exit::THETA : Algorithm_TD_Slerp::Exit::THETA_PI ) );
                    CHECK( tResult.mEvaluations == 2 );    // the start as given and normalized
                    CHECK( tResult.mAccepted == 0 );
                    CHECK( tInterface->mFetchAfter.size() == 1 );
                    CHECK( std::abs( tResult.mThetas[ 0 ] - ( tMode == td_slerp::Mode::FIXED_POINT ? 0.0 : 180.0 ) ) < 1.0e-5 );
                }
            }
        }
    }

    TEST_CASE( "TD_SLERP line search halves kappa to kappa_min, exits LINE_SEARCH and restores phi_i", "[optimization][td_slerp][U16]" )
    {
        for ( bool tIdentity : { true, false } )
        {
            std::shared_ptr< td_slerp::Mock_Interface > tInterface = std::make_shared< td_slerp::Mock_Interface >(
                    td_slerp::mass_matrix( tIdentity ), td_slerp::gStart, td_slerp::gA, td_slerp::Mode::RISING );
            td_slerp::Result tResult = td_slerp::run( tInterface, moris::prm::create_td_slerp_parameter_list() );

            if ( par_rank() == 0 )
            {
                CHECK( tResult.mExit == Algorithm_TD_Slerp::Exit::LINE_SEARCH );
                CHECK( tResult.mAccepted == 0 );

                // the start (as given, normalized) plus the trials k = 1, 1/2, ..., 2^-10 (2^-10 = kappa_min is tried,
                // 2^-11 is not)
                CHECK( tResult.mEvaluations == 13 );
                CHECK( tInterface->mFetchAfter.size() == 1 );

                // the problem reports J_i, the last restart file holds phi_i = the normalized start
                CHECK( tResult.mFinalObjective == tInterface->mObjectives[ 1 ] );
                // the files are numbered by evaluation from 0 in each solve: 0 .. 12, the restored design in 13
                std::vector< real > tRestored = td_slerp::read_restart_advs( tResult.mEvaluations );
                REQUIRE( tRestored.size() == 4 );
                for ( uint i = 0; i < 4; i++ )
                {
                    CHECK( tRestored[ i ] == tInterface->mEvaluated[ 1 ][ i ] );
                }
                real tNorm0 = tInterface->mnorm( td_slerp::gStart );
                for ( uint i = 0; i < 4; i++ )
                {
                    CHECK( std::abs( tRestored[ i ] - td_slerp::gStart( i ) / tNorm0 ) < 1.0e-15 );
                }
            }
        }
    }

    TEST_CASE( "TD_SLERP is sign-flip equivariant and fetches g~ only after accepted evaluations", "[optimization][td_slerp][U17]" )
    {
        for ( bool tIdentity : { true, false } )
        {
            Vector< real > tFlippedStart = td_slerp::gStart;
            for ( uint i = 0; i < 4; i++ ) tFlippedStart( i ) = -tFlippedStart( i );

            std::shared_ptr< td_slerp::Mock_Interface > tPlain = std::make_shared< td_slerp::Mock_Interface >(
                    td_slerp::mass_matrix( tIdentity ), td_slerp::gStart, td_slerp::gA, td_slerp::Mode::OVERSHOOT );
            std::shared_ptr< td_slerp::Mock_Interface > tFlipped = std::make_shared< td_slerp::Mock_Interface >(
                    td_slerp::mass_matrix( tIdentity ), tFlippedStart, td_slerp::gA, td_slerp::Mode::OVERSHOOT, -1.0 );

            td_slerp::Result tResultPlain   = td_slerp::run( tPlain, moris::prm::create_td_slerp_parameter_list() );
            td_slerp::Result tResultFlipped = td_slerp::run( tFlipped, moris::prm::create_td_slerp_parameter_list() );

            if ( par_rank() == 0 )
            {
                // (phi, g~) -> (-phi, -g~): the same decisions and the negated iterates
                CHECK( tResultFlipped.mExit == tResultPlain.mExit );
                CHECK( tResultFlipped.mAccepted == tResultPlain.mAccepted );
                REQUIRE( tFlipped->mEvaluated.size() == tPlain->mEvaluated.size() );
                for ( uint k = 0; k < tPlain->mEvaluated.size(); k++ )
                {
                    CHECK( tFlipped->mObjectives[ k ] == Approx( tPlain->mObjectives[ k ] ).epsilon( 1.0e-13 ) );
                    for ( uint i = 0; i < 4; i++ )
                    {
                        CHECK( std::abs( tFlipped->mEvaluated[ k ][ i ] + tPlain->mEvaluated[ k ][ i ] ) < 1.0e-14 );
                    }
                }

                // g~ fetched once after the start and once after every accepted evaluation, never after a rejected one
                CHECK( tPlain->mFetchAfter.size() == tResultPlain.mAccepted + 1 );
                CHECK( tPlain->mFetchAfter[ 0 ] == 2 );
                real tAcceptedJ = tPlain->mObjectives[ 1 ];
                for ( uint k = 1; k < tPlain->mFetchAfter.size(); k++ )
                {
                    uint tLastEvaluation = tPlain->mFetchAfter[ k ] - 1;
                    CHECK( tPlain->mObjectives[ tLastEvaluation ] <= tAcceptedJ * ( 1.0 + 1.0e-10 ) );
                    tAcceptedJ = tPlain->mObjectives[ tLastEvaluation ];
                }
            }
        }
    }

    TEST_CASE( "TD_SLERP restart_index resumes the accepted-iteration count; a bound hit is fatal", "[optimization][td_slerp][U17]" )
    {
        Parameter_List tParameterList = moris::prm::create_td_slerp_parameter_list();
        tParameterList.set( "restart_index", 2 );
        tParameterList.set( "max_its", 3 );

        // a restart design is on the unit sphere already: evaluated once, not rescaled
        td_slerp::Mock_Interface tScale( td_slerp::mass_matrix( false ), td_slerp::gStart, td_slerp::gA, td_slerp::Mode::OVERSHOOT );
        Vector< real >           tUnitStart = td_slerp::gStart;
        real                     tNorm      = tScale.mnorm( td_slerp::gStart );
        for ( uint i = 0; i < 4; i++ ) tUnitStart( i ) = td_slerp::gStart( i ) / tNorm;

        std::shared_ptr< td_slerp::Mock_Interface > tInterface = std::make_shared< td_slerp::Mock_Interface >(
                td_slerp::mass_matrix( false ), tUnitStart, td_slerp::gA, td_slerp::Mode::OVERSHOOT );
        td_slerp::Result tResult = td_slerp::run( tInterface, tParameterList );

        if ( par_rank() == 0 )
        {
            CHECK( tResult.mExit == Algorithm_TD_Slerp::Exit::N_MAX );
            CHECK( tResult.mAccepted == 3 );
            CHECK( tResult.mKappas.size() == 1 );      // one accepted iteration in this run
            CHECK( tResult.mEvaluations == 1 + 2 );    // the start once, then k = 1 (rejected) and k = 1/2
            for ( uint i = 0; i < 4; i++ )
            {
                CHECK( tInterface->mEvaluated[ 0 ][ i ] == tUnitStart( i ) );
            }
        }

        // bounds +-10 with |phi|_M = 1: a start that is not normalized and beyond a bound
        Parameter_List tNoNormalize = moris::prm::create_td_slerp_parameter_list();
        tNoNormalize.set( "normalize_start", false );
        Vector< real >                              tFar          = { 20.0, 0.0, 0.0, 0.0 };
        std::shared_ptr< td_slerp::Mock_Interface > tFarInterface = std::make_shared< td_slerp::Mock_Interface >(
                td_slerp::mass_matrix( true ), tFar, td_slerp::gA, td_slerp::Mode::OVERSHOOT );
        if ( par_size() == 1 )
        {
            CHECK_THROWS( td_slerp::run( tFarInterface, tNoNormalize ) );
        }
    }
}    // namespace moris::opt
