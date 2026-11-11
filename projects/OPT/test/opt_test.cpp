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

                // objectives live on rank 0 only (Problem::get_objectives errors elsewhere); the callers check on rank 0 (MB-088)
                real tObjective = par_rank() == 0 ? tManager.get_objectives()( 0 ) : 0.0;
                return { tObjective, tManager.get_advs() };
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
            // only processor 0 holds the ADVs; the others evaluate an empty dummy vector (Algorithm::dummy_solve)
            if ( par_rank() != 0 )
            {
                return { {} };
            }

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
            if ( par_rank() != 0 )
            {
                return { {} };
            }

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
        // every rank enters solve (the others serve criteria requests in Algorithm::dummy_solve); only rank 0 checks
        // evaluation 1 is the initial design, 2 the (non-conservative) first trial, 3 the first inner trial
        failed_evaluation::reset( failed_evaluation::Fault::single_evaluation, 3 );
        if ( par_rank() == 0 )
        {
            failed_evaluation::remove_restart_files();
        }
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

        if ( par_rank() == 0 )
        {
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
        // every rank enters solve (the others serve criteria requests in Algorithm::dummy_solve); only rank 0 checks
        // every trial fails, however short the step; only the initial design evaluates
        failed_evaluation::reset( failed_evaluation::Fault::off_initial_design );
        if ( par_rank() == 0 )
        {
            failed_evaluation::remove_restart_files();
        }
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

        if ( par_rank() == 0 )
        {
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

    //--------------------------------------------------------------------------------------------------------------
    // Warm MMA restart (restart_keep_state): min sum_i w_i (x_i - t_i)^2 s.t. sum_i x_i <= 0.4 N, 0 <= x <= 1, for any
    // number of ADVs N. Criteria = the ADVs. The interface can request an optimizer restart at a given evaluation (with
    // NaN criteria, as the workflow's reinitialization schedule does) and re-initializes at the last evaluated design
    // plus a shift (a "re-fit" that moves some ADVs).

    namespace warm_restart
    {
        real
        weight( uint aI )
        {
            return 1.0 + 0.7 * aI;
        }

        real
        target( uint aI )
        {
            return ( aI % 2 == 0 ) ? 0.9 : 0.1 + 0.05 * aI;
        }

        void
        unused_initialize( Vector< real >&, Vector< real >&, Vector< real >& )
        {
        }

        Vector< real >
        get_criteria( const Vector< real >& aADVs )
        {
            return aADVs;
        }

        Matrix< DDRMat >
        get_dcriteria_dadv( const Vector< real >& aADVs )
        {
            Matrix< DDRMat > tIdentity( aADVs.size(), aADVs.size(), 0.0 );
            for ( uint iADV = 0; iADV < aADVs.size(); iADV++ )
            {
                tIdentity( iADV, iADV ) = 1.0;
            }
            return tIdentity;
        }

        Matrix< DDSMat >
        get_constraint_types()
        {
            return Matrix< DDSMat >( 1, 1, 1 );
        }

        Matrix< DDRMat >
        compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            real tObjective = 0.0;
            for ( uint iADV = 0; iADV < aCriteria.size(); iADV++ )
            {
                tObjective += weight( iADV ) * std::pow( aCriteria( iADV ) - target( iADV ), 2 );
            }
            return { { tObjective } };
        }

        Matrix< DDRMat >
        compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            real tSum = 0.0;
            for ( uint iADV = 0; iADV < aCriteria.size(); iADV++ )
            {
                tSum += aCriteria( iADV );
            }
            return { { tSum - 0.4 * aCriteria.size() } };
        }

        Matrix< DDRMat >
        zero_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return Matrix< DDRMat >( 1, aADVs.size(), 0.0 );
        }

        Matrix< DDRMat >
        compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            Matrix< DDRMat > tGradient( 1, aCriteria.size(), 0.0 );
            for ( uint iADV = 0; iADV < aCriteria.size(); iADV++ )
            {
                tGradient( 0, iADV ) = 2.0 * weight( iADV ) * ( aCriteria( iADV ) - target( iADV ) );
            }
            return tGradient;
        }

        Matrix< DDRMat >
        compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            return Matrix< DDRMat >( 1, aCriteria.size(), 1.0 );
        }

        /** restart at evaluation mRestartAtCall (0: never); the restart design is the last evaluated one plus mShift */
        class Interface_Forced_Restart : public Interface_User_Defined
        {
          public:
            uint                               mNumADVs         = 4;
            uint                               mRestartAtCall   = 0;
            uint                               mCalls           = 0;
            uint                               mInitializations = 0;
            std::vector< real >                mShift;
            std::vector< real >                mLast;
            std::vector< std::vector< real > > mEvaluated;    ///< every finite evaluation's design

            Interface_Forced_Restart()
                    : Interface_User_Defined( &warm_restart::unused_initialize, &warm_restart::get_criteria, &warm_restart::get_dcriteria_dadv )
            {
            }

            void
            initialize( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds, Matrix< IdMat >& ) override
            {
                mInitializations++;
                aADVs.resize( mNumADVs );
                aLowerBounds.resize( mNumADVs );
                aUpperBounds.resize( mNumADVs );
                for ( uint iADV = 0; iADV < mNumADVs; iADV++ )
                {
                    real tValue = 0.5;
                    if ( mInitializations > 1 )
                    {
                        tValue = mLast[ iADV ] + ( iADV < mShift.size() ? mShift[ iADV ] : 0.0 );
                    }
                    aADVs( iADV )        = std::min( std::max( tValue, 0.0 ), 1.0 );
                    aLowerBounds( iADV ) = 0.0;
                    aUpperBounds( iADV ) = 1.0;
                }
                mInitializeOptimizationRestart = false;
            }

            Vector< real >
            perform( Vector< real >& aNewADVs ) override
            {
                mCalls++;
                if ( mCalls == mRestartAtCall )
                {
                    mInitializeOptimizationRestart = true;
                    return Vector< real >( aNewADVs.size(), std::numeric_limits< real >::quiet_NaN() );
                }
                mLast.assign( aNewADVs.size(), 0.0 );
                for ( uint iADV = 0; iADV < aNewADVs.size(); iADV++ )
                {
                    mLast[ iADV ] = aNewADVs( iADV );
                }
                mEvaluated.push_back( mLast );
                return Interface_User_Defined::perform( aNewADVs );
            }
        };

        std::shared_ptr< Problem >
        make_problem( const std::shared_ptr< Interface_Forced_Restart >& aInterface )
        {
            Parameter_List tProblemParameterList = moris::prm::create_opt_problem_parameter_list();

            std::shared_ptr< Problem > tProblem = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    aInterface,
                    &get_constraint_types,
                    &compute_objectives,
                    &compute_constraints,
                    &zero_dadv,
                    &compute_dobjective_dcriteria,
                    &zero_dadv,
                    &compute_dconstraint_dcriteria );
            tProblem->initialize();
            return tProblem;
        }

        Parameter_List
        algorithm_parameters( bool aKeepState, sint aMaxIts )
        {
            Parameter_List tParameterList = moris::prm::create_mma_parameter_list();
            tParameterList.set( "max_its", aMaxIts );
            tParameterList.set( "norm_drop", 1.0e-14 );
            tParameterList.set( "step_size", 0.1 );
            tParameterList.set( "move_adapt", true );
            tParameterList.set( "move_adapt_init", 0.05 );
            tParameterList.set( "move_adapt_min", 1.0e-3 );
            tParameterList.set( "restart_keep_state", aKeepState );
            return tParameterList;
        }

        /** the Manager's restart loop (Manager::reinitialize), keeping access to the algorithm */
        void
        solve_with_restarts( Algorithm_MMA& aAlgorithm, const std::shared_ptr< Problem >& aProblem )
        {
            Tracer tTracer( "OPT", "Manager", "Perform" );
            aAlgorithm.solve( 0, aProblem );
            while ( aProblem->restart_optimization() )
            {
                aProblem->initialize();
                aAlgorithm.solve( 0, aProblem );
            }
        }

        void
        remove_state_files()
        {
            for ( const std::filesystem::directory_entry& tEntry : std::filesystem::directory_iterator( "." ) )
            {
                const std::string tName = tEntry.path().filename().string();
                if ( tName.rfind( "MMA_State_Alg_0_Iter_", 0 ) == 0 || tName.rfind( "ADV_Alg_0_Iter_", 0 ) == 0 )
                {
                    std::filesystem::remove( tEntry.path() );
                }
            }
        }
    }    // namespace warm_restart

    TEST_CASE( "Native MMA restart state is re-centred, row-mapped and reset for moved ADVs", "[optimization][restart_keep_state]" )
    {
        Algorithm_MMA::Restart_State tState;
        tState.mValid    = true;
        tState.mCentre   = { 0.5, 0.5, 0.5, 0.5 };
        tState.mLowWidth = { 0.10, 0.20, 0.001, 0.30 };
        tState.mUppWidth = { 0.20, 0.30, 0.40, 0.30 };
        tState.mMoveStep = { 0.01, 0.02, 0.03, 0.04 };

        // identity map; ADV 1 has range 2 (move step 0.02 -> fraction 0.01); ADV 2's low width is clamped to 0.01 range;
        // ADV 3 moved by 0.45 > 0.25 of its range: fresh
        std::vector< double > tX    = { 0.55, 0.5, 0.5, 0.95 };
        std::vector< double > tXMin = { 0.0, -1.0, 0.0, 0.0 };
        std::vector< double > tXMax = { 1.0, 1.0, 1.0, 1.0 };
        std::vector< double > tLow, tUpp, tMove;

        Algorithm_MMA::Restart_Counts tCounts = Algorithm_MMA::apply_restart_state(
                tState, {}, tX, tXMin, tXMax, 0.25, 0.5, true, 0.07, 1.0e-4, 0.1, tLow, tUpp, tMove );

        CHECK( tCounts.mKept == 3 );
        CHECK( tCounts.mMoved == 1 );
        CHECK( tCounts.mFresh == 0 );
        CHECK( tLow[ 0 ] == Approx( 0.45 ) );
        CHECK( tUpp[ 0 ] == Approx( 0.75 ) );
        CHECK( tMove[ 0 ] == Approx( 0.01 ) );
        CHECK( tLow[ 1 ] == Approx( 0.30 ) );
        CHECK( tUpp[ 1 ] == Approx( 0.80 ) );
        CHECK( tMove[ 1 ] == Approx( 0.01 ) );
        CHECK( tLow[ 2 ] == Approx( 0.49 ) );
        CHECK( tUpp[ 2 ] == Approx( 0.90 ) );
        CHECK( tMove[ 2 ] == Approx( 0.03 ) );
        CHECK( tLow[ 3 ] == Approx( 0.45 ) );
        CHECK( tUpp[ 3 ] == Approx( 1.45 ) );
        CHECK( tMove[ 3 ] == Approx( 0.07 ) );

        // a fixed move box (no move_adapt) is left untouched
        std::vector< double > tFixed( 4, 0.1 );
        Algorithm_MMA::apply_restart_state( tState, {}, tX, tXMin, tXMax, 0.25, 0.5, false, 0.07, 1.0e-4, 0.1, tLow, tUpp, tFixed );
        CHECK( tFixed == std::vector< double >( 4, 0.1 ) );

        // row map (3 new ADVs): new 0 <- old 3, new 1 fresh, new 2 <- old 0; with the moved test off, any distance keeps
        std::vector< double > tX3    = { 0.1, 0.5, 0.9 };
        std::vector< double > tXMin3 = { 0.0, 0.0, 0.0 };
        std::vector< double > tXMax3 = { 1.0, 1.0, 1.0 };
        tCounts                      = Algorithm_MMA::apply_restart_state(
                tState, { 3, -1, 0 }, tX3, tXMin3, tXMax3, 0.0, 0.5, true, 0.07, 1.0e-4, 0.1, tLow, tUpp, tMove );
        REQUIRE( tLow.size() == 3 );
        CHECK( tCounts.mKept == 2 );
        CHECK( tCounts.mFresh == 1 );
        CHECK( tLow[ 0 ] == Approx( 0.1 - 0.30 ) );
        CHECK( tUpp[ 0 ] == Approx( 0.1 + 0.30 ) );
        CHECK( tMove[ 0 ] == Approx( 0.04 ) );
        CHECK( tLow[ 1 ] == Approx( 0.0 ) );
        CHECK( tUpp[ 1 ] == Approx( 1.0 ) );
        CHECK( tMove[ 1 ] == Approx( 0.07 ) );
        CHECK( tLow[ 2 ] == Approx( 0.8 ) );
        CHECK( tUpp[ 2 ] == Approx( 1.1 ) );
        CHECK( tMove[ 2 ] == Approx( 0.01 ) );

        // a map entry outside [-1, number of kept ADVs) is an error
        CHECK_THROWS( Algorithm_MMA::apply_restart_state(
                tState, { 4, 0, 1 }, tX3, tXMin3, tXMax3, 0.0, 0.5, true, 0.07, 1.0e-4, 0.1, tLow, tUpp, tMove ) );
    }

    TEST_CASE( "Native MMA keeps its state across an in-process restart", "[optimization][restart_keep_state]" )
    {
        // one processor only: the test drives solve() and its restarts on a single rank, without the other ranks'
        // Algorithm::dummy_solve partner (MB-091)
        if ( par_size() > 1 )
        {
            return;
        }
        using namespace warm_restart;

        // the restart is requested at evaluation 6 (the trial of iteration 5: 4 completed iterations); the re-fit moves
        // ADV 1 by 0.003 (kept) and ADV 3 by 0.4 (> 0.25 of its range: fresh)
        std::shared_ptr< Interface_Forced_Restart > tInterface = std::make_shared< Interface_Forced_Restart >();
        tInterface->mRestartAtCall                             = 6;
        tInterface->mShift                                     = { 0.0, 0.003, 0.0, 0.4 };
        std::shared_ptr< Problem > tProblem                    = make_problem( tInterface );

        Algorithm_MMA                tAlgorithm( algorithm_parameters( true, 40 ) );
        Algorithm_MMA::Restart_State tBefore;
        std::vector< real >          tRestartX;
        {
            Tracer tTracer( "OPT", "Manager", "Perform" );
            tAlgorithm.solve( 0, tProblem );
            REQUIRE( tProblem->restart_optimization() );
            REQUIRE( tAlgorithm.get_num_warm_solves() == 0 );
            tBefore = tAlgorithm.get_kept_state();
            REQUIRE( tBefore.mValid );

            tProblem->initialize();
            for ( uint iADV = 0; iADV < 4; iADV++ )
            {
                tRestartX.push_back( tProblem->get_advs()( iADV ) );
            }

            // one iteration after the restart (Manager::reinitialize sets the restart index the same way)
            tAlgorithm.set_restart_index( 39 );
            tAlgorithm.solve( 0, tProblem );
        }

        REQUIRE( tAlgorithm.get_num_warm_solves() == 1 );
        const Algorithm_MMA::Restart_State&  tStart  = tAlgorithm.get_warm_start_state();
        const Algorithm_MMA::Restart_Counts& tCounts = tAlgorithm.get_warm_start_counts();
        CHECK( tCounts.mKept == 3 );
        CHECK( tCounts.mMoved == 1 );
        CHECK( tCounts.mFresh == 0 );

        // the kept quantities equal the pre-restart ones, the widths re-centred on the restart design
        CHECK( tStart.mRaa0 == tBefore.mRaa0 );
        CHECK( tStart.mRaa == tBefore.mRaa );
        CHECK( tStart.mPenalty == tBefore.mPenalty );
        CHECK( tStart.mObjScale == tBefore.mObjScale );
        CHECK( tStart.mKKTRef == tBefore.mKKTRef );
        bool tMovesDiffer = false;
        for ( uint iADV = 0; iADV < 4; iADV++ )
        {
            CHECK( tStart.mCentre[ iADV ] == tRestartX[ iADV ] );
            if ( iADV == 3 )
            {
                // fresh: asymptotes at x -+ 0.5 range, the initial move limit
                CHECK( tStart.mLowWidth[ iADV ] == Approx( 0.5 ) );
                CHECK( tStart.mUppWidth[ iADV ] == Approx( 0.5 ) );
                CHECK( tStart.mMoveStep[ iADV ] == Approx( 0.05 ) );
                continue;
            }
            CHECK( tStart.mLowWidth[ iADV ] == Approx( tBefore.mLowWidth[ iADV ] ).epsilon( 1e-12 ) );
            CHECK( tStart.mUppWidth[ iADV ] == Approx( tBefore.mUppWidth[ iADV ] ).epsilon( 1e-12 ) );
            CHECK( tStart.mMoveStep[ iADV ] == tBefore.mMoveStep[ iADV ] );
            tMovesDiffer = tMovesDiffer || tBefore.mMoveStep[ iADV ] != 0.05;
            // not the reset of a fresh solve
            CHECK( std::abs( tStart.mLowWidth[ iADV ] - 0.5 ) > 1e-3 );
        }
        CHECK( tMovesDiffer );    // the move limits had adapted before the restart

        // the first iteration after the restart used the kept widths (history branch, factor 1): its own next widths
        // (x1 - x' and x' - x' give factor 1 again) are the kept ones, and its step stays inside the kept move limits
        const Algorithm_MMA::Restart_State& tAfter = tAlgorithm.get_kept_state();
        const std::vector< real >&          tX1    = tInterface->mEvaluated.back();
        for ( uint iADV = 0; iADV < 3; iADV++ )
        {
            CHECK( tAfter.mLowWidth[ iADV ] == Approx( tStart.mLowWidth[ iADV ] ).epsilon( 1e-12 ) );
            CHECK( tAfter.mUppWidth[ iADV ] == Approx( tStart.mUppWidth[ iADV ] ).epsilon( 1e-12 ) );
            CHECK( std::abs( tX1[ iADV ] - tRestartX[ iADV ] ) <= tStart.mMoveStep[ iADV ] * ( 1.0 + 1e-12 ) );
        }
        warm_restart::remove_state_files();
    }

    TEST_CASE( "Native MMA warm restart: key off is unchanged, and the restarted run converges", "[optimization][restart_keep_state]" )
    {
        // one processor only: the test drives solve() and its restarts on a single rank, without the other ranks'
        // Algorithm::dummy_solve partner (MB-091)
        if ( par_size() > 1 )
        {
            return;
        }
        using namespace warm_restart;

        // no restart: restart_keep_state changes nothing (same designs, bitwise)
        std::vector< std::vector< std::vector< real > > > tTrajectories;
        for ( bool tKeep : { false, true } )
        {
            std::shared_ptr< Interface_Forced_Restart > tInterface = std::make_shared< Interface_Forced_Restart >();
            std::shared_ptr< Problem >                  tProblem   = make_problem( tInterface );
            Algorithm_MMA                               tAlgorithm( algorithm_parameters( tKeep, 80 ) );
            solve_with_restarts( tAlgorithm, tProblem );
            tTrajectories.push_back( tInterface->mEvaluated );
        }
        CHECK( tTrajectories[ 0 ] == tTrajectories[ 1 ] );

        // with a restart (and a re-fit shift) the run converges to the same optimum, reset or kept
        std::vector< std::vector< real > > tFinal;
        for ( bool tKeep : { false, true } )
        {
            std::shared_ptr< Interface_Forced_Restart > tInterface = std::make_shared< Interface_Forced_Restart >();
            tInterface->mRestartAtCall                             = 6;
            tInterface->mShift                                     = { 0.0, 0.003, 0.0, 0.4 };
            std::shared_ptr< Problem > tProblem                    = make_problem( tInterface );
            Algorithm_MMA              tAlgorithm( algorithm_parameters( tKeep, 80 ) );
            solve_with_restarts( tAlgorithm, tProblem );
            CHECK( tAlgorithm.get_num_warm_solves() == ( tKeep ? 1 : 0 ) );
            tFinal.push_back( tInterface->mEvaluated.back() );
        }
        // KKT optimum: x_i = clip( t_i - lambda / ( 2 w_i ), 0, 1 ), sum x_i = 0.4 N (bisection on lambda)
        real tLambdaLow  = 0.0;
        real tLambdaHigh = 100.0;
        for ( uint iBisect = 0; iBisect < 200; iBisect++ )
        {
            const real tLambda = 0.5 * ( tLambdaLow + tLambdaHigh );
            real       tSum    = 0.0;
            for ( uint iADV = 0; iADV < 4; iADV++ )
            {
                tSum += std::min( std::max( target( iADV ) - tLambda / ( 2.0 * weight( iADV ) ), 0.0 ), 1.0 );
            }
            ( tSum > 1.6 ? tLambdaLow : tLambdaHigh ) = tLambda;
        }
        for ( uint iADV = 0; iADV < 4; iADV++ )
        {
            const real tOptimum = std::min( std::max( target( iADV ) - tLambdaLow / ( 2.0 * weight( iADV ) ), 0.0 ), 1.0 );
            CHECK( tTrajectories[ 0 ].back()[ iADV ] == Approx( tOptimum ).margin( 1e-4 ) );
            CHECK( tFinal[ 0 ][ iADV ] == Approx( tOptimum ).margin( 1e-4 ) );
            CHECK( tFinal[ 1 ][ iADV ] == Approx( tOptimum ).margin( 1e-4 ) );
        }
        warm_restart::remove_state_files();
    }

    TEST_CASE( "Native MMA restarts a new run from a state file with a row map", "[optimization][restart_keep_state]" )
    {
        // one processor only: the test drives solve() and its restarts on a single rank, without the other ranks'
        // Algorithm::dummy_solve partner (MB-091)
        if ( par_size() > 1 )
        {
            return;
        }
        using namespace warm_restart;
        warm_restart::remove_state_files();

        // run 1 (4 ADVs) writes the state next to every committed design
        std::shared_ptr< Interface_Forced_Restart > tInterface = std::make_shared< Interface_Forced_Restart >();
        std::shared_ptr< Problem >                  tProblem   = make_problem( tInterface );
        Algorithm_MMA                               tAlgorithm( algorithm_parameters( true, 6 ) );
        solve_with_restarts( tAlgorithm, tProblem );

        std::string tLatest;
        uint        tLatestIter = 0;
        uint        tNumFiles   = 0;
        for ( const std::filesystem::directory_entry& tEntry : std::filesystem::directory_iterator( "." ) )
        {
            const std::string tName = tEntry.path().filename().string();
            if ( tName.rfind( "MMA_State_Alg_0_Iter_", 0 ) == 0 )
            {
                tNumFiles++;
                uint tIter = (uint)std::stoul( tName.substr( 21 ) );
                if ( tLatest.empty() || tIter > tLatestIter )
                {
                    tLatest     = tName;
                    tLatestIter = tIter;
                }
            }
        }
        REQUIRE( tNumFiles == 6 );    // one per completed iteration
        CHECK( std::filesystem::exists( "ADV_Alg_0_Iter_" + std::to_string( tLatestIter ) + ".hdf5" ) );

        // the file holds the kept state of the last iteration, and its design is the committed one
        Algorithm_MMA::Restart_State        tFile = Algorithm_MMA::read_restart_state( tLatest );
        const Algorithm_MMA::Restart_State& tKept = tAlgorithm.get_kept_state();
        CHECK( tFile.mCentre == tKept.mCentre );
        CHECK( tFile.mLowWidth == tKept.mLowWidth );
        CHECK( tFile.mUppWidth == tKept.mUppWidth );
        CHECK( tFile.mMoveStep == tKept.mMoveStep );
        CHECK( tFile.mRaa0 == tKept.mRaa0 );
        CHECK( tFile.mRaa == tKept.mRaa );
        CHECK( tFile.mObjScale == tKept.mObjScale );
        CHECK( tFile.mCentre == tInterface->mEvaluated.back() );

        // run 2: 5 ADVs (a "prune and grow"), new ADV i <- old ADV { 3, -1, 0, 1, 2 }[ i ]
        std::shared_ptr< Interface_Forced_Restart > tInterface2 = std::make_shared< Interface_Forced_Restart >();
        tInterface2->mNumADVs                                   = 5;
        std::shared_ptr< Problem > tProblem2                    = make_problem( tInterface2 );
        Parameter_List             tParameters                  = algorithm_parameters( true, 1 );
        tParameters.set( "restart_state_file", tLatest );
        tParameters.set( "restart_state_row_map", "3, -1, 0, 1, 2" );
        tParameters.set( "restart_keep_state_moved_fraction", 0.0 );
        Algorithm_MMA tAlgorithm2( tParameters );
        solve_with_restarts( tAlgorithm2, tProblem2 );

        REQUIRE( tAlgorithm2.get_num_warm_solves() == 1 );
        const Algorithm_MMA::Restart_State&  tStart  = tAlgorithm2.get_warm_start_state();
        const Algorithm_MMA::Restart_Counts& tCounts = tAlgorithm2.get_warm_start_counts();
        CHECK( tCounts.mKept == 4 );
        CHECK( tCounts.mFresh == 1 );
        CHECK( tStart.mRaa0 == tFile.mRaa0 );
        CHECK( tStart.mObjScale == tFile.mObjScale );
        const std::vector< sint > tMap = { 3, -1, 0, 1, 2 };
        for ( uint iADV = 0; iADV < 5; iADV++ )
        {
            CHECK( tStart.mCentre[ iADV ] == 0.5 );
            if ( tMap[ iADV ] < 0 )
            {
                CHECK( tStart.mLowWidth[ iADV ] == Approx( 0.5 ) );
                CHECK( tStart.mMoveStep[ iADV ] == Approx( 0.05 ) );
                continue;
            }
            CHECK( tStart.mLowWidth[ iADV ] == Approx( tFile.mLowWidth[ tMap[ iADV ] ] ).epsilon( 1e-12 ) );
            CHECK( tStart.mUppWidth[ iADV ] == Approx( tFile.mUppWidth[ tMap[ iADV ] ] ).epsilon( 1e-12 ) );
            CHECK( tStart.mMoveStep[ iADV ] == tFile.mMoveStep[ tMap[ iADV ] ] );
        }

        // without a row map, a state file of another size is an error
        Parameter_List tWrongSize = algorithm_parameters( true, 1 );
        tWrongSize.set( "restart_state_file", tLatest );
        std::shared_ptr< Interface_Forced_Restart > tInterface3 = std::make_shared< Interface_Forced_Restart >();
        tInterface3->mNumADVs                                   = 5;
        std::shared_ptr< Problem > tProblem3                    = make_problem( tInterface3 );
        Algorithm_MMA              tAlgorithm3( tWrongSize );
        {
            Tracer tTracer( "OPT", "Manager", "Perform" );
            CHECK_THROWS( tAlgorithm3.solve( 0, tProblem3 ) );
        }

        // a state file needs restart_keep_state
        Parameter_List tNoKeep = algorithm_parameters( false, 1 );
        tNoKeep.set( "restart_state_file", tLatest );
        CHECK_THROWS( Algorithm_MMA( tNoKeep ) );

        warm_restart::remove_state_files();
    }

    //--------------------------------------------------------------------------------------------------------------
    // The optimizer state file identifies its ADVs by ID (MB-114): the ADV order depends on the processor count, so a
    // state written in one order must land on the same ADVs in another. The tests emulate the orders of 1, 2 and 4
    // processors by permutations of the ADV IDs (processor by processor, the owned ADVs first); the process count of
    // the test itself does not change the order of a user-defined problem.

    namespace state_ids
    {
        /** ID of the ADV at each position, in the order emulated for aNumProcs processors (8 ADVs, IDs 100 .. 107) */
        std::vector< moris_id >
        order( sint aNumProcs )
        {
            switch ( aNumProcs )
            {
                case 1:
                    return { 100, 101, 102, 103, 104, 105, 106, 107 };
                case 2:
                    return { 104, 105, 106, 107, 100, 101, 102, 103 };
                default:
                    return { 107, 100, 104, 103, 101, 106, 105, 102 };
            }
        }

        Matrix< IdMat >
        id_matrix( const std::vector< moris_id >& aIds )
        {
            Matrix< IdMat > tIds( aIds.size(), 1, gNoID );
            for ( uint iADV = 0; iADV < aIds.size(); iADV++ )
            {
                tIds( iADV ) = aIds[ iADV ];
            }
            return tIds;
        }

        // per-ADV values as functions of the ID, so that a value found at a position names the ADV it belongs to
        real
        centre_of( moris_id aId )
        {
            return 0.1 + 0.1 * ( aId - 100 );
        }

        real
        low_width_of( moris_id aId )
        {
            return 0.02 + 0.01 * ( aId - 100 );
        }

        real
        upp_width_of( moris_id aId )
        {
            return 0.2 + 0.03 * ( aId - 100 );
        }

        real
        move_step_of( moris_id aId )
        {
            return 0.011 + 0.002 * ( aId - 100 );
        }

        /** a state in the given order; with aWithIds = false as a problem without ADV IDs writes it */
        Algorithm_MMA::Restart_State
        make_state( const std::vector< moris_id >& aIds, sint aNumProcs, bool aWithIds )
        {
            Algorithm_MMA::Restart_State tState;
            tState.mValid = true;
            for ( moris_id tId : aIds )
            {
                tState.mCentre.push_back( centre_of( tId ) );
                tState.mLowWidth.push_back( low_width_of( tId ) );
                tState.mUppWidth.push_back( upp_width_of( tId ) );
                tState.mMoveStep.push_back( move_step_of( tId ) );
                tState.mLowerBounds.push_back( 0.0 );
                tState.mUpperBounds.push_back( 1.0 );
            }
            tState.mRaa      = { 1.0e-3 };
            tState.mPenalty  = { 7.0 };
            tState.mRaa0     = 2.0e-3;
            tState.mObjScale = 0.5;
            tState.mKKTRef   = 0.25;
            if ( aWithIds )
            {
                tState.mIjklIds = id_matrix( aIds );
            }
            tState.mNumProcs = aNumProcs;
            return tState;
        }

        /** a state file as written before the IDs and the processor count were stored */
        void
        write_old_format_state( const std::string& aFileName, const Algorithm_MMA::Restart_State& aState )
        {
            hid_t  tFileID = create_hdf5_file( aFileName, false );
            herr_t tStatus = 0;
            save_vector_to_hdf5_file( tFileID, "Centre", aState.mCentre, tStatus );
            save_vector_to_hdf5_file( tFileID, "LowWidth", aState.mLowWidth, tStatus );
            save_vector_to_hdf5_file( tFileID, "UppWidth", aState.mUppWidth, tStatus );
            save_vector_to_hdf5_file( tFileID, "MoveStep", aState.mMoveStep, tStatus );
            save_vector_to_hdf5_file( tFileID, "LowerBounds", aState.mLowerBounds, tStatus );
            save_vector_to_hdf5_file( tFileID, "UpperBounds", aState.mUpperBounds, tStatus );
            save_vector_to_hdf5_file( tFileID, "Raa", aState.mRaa, tStatus );
            save_vector_to_hdf5_file( tFileID, "Penalty", aState.mPenalty, tStatus );
            save_scalar_to_hdf5_file( tFileID, "Raa0", aState.mRaa0, tStatus );
            save_scalar_to_hdf5_file( tFileID, "ObjScale", aState.mObjScale, tStatus );
            save_scalar_to_hdf5_file( tFileID, "KKTRef", aState.mKKTRef, tStatus );
            save_scalar_to_hdf5_file( tFileID, "OptIteration", (sint)3, tStatus );
            close_hdf5_file( tFileID );
        }

        /** the design with every ADV at the centre the state holds for its ID, in the given order */
        std::vector< double >
        design_at_centres( const std::vector< moris_id >& aIds )
        {
            std::vector< double > tX;
            for ( moris_id tId : aIds )
            {
                tX.push_back( centre_of( tId ) );
            }
            return tX;
        }

        // ---- a problem whose ADVs are identified by their IDs: min sum_i w( id_i ) ( x_i - t( id_i ) )^2 ----

        /** ID of the ADV at each position of the problem that is being solved */
        std::vector< moris_id > gSolveOrder;

        Matrix< DDRMat >
        compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            real tObjective = 0.0;
            for ( uint iADV = 0; iADV < aCriteria.size(); iADV++ )
            {
                const uint tVariable = gSolveOrder[ iADV ] - 100;
                tObjective += warm_restart::weight( tVariable ) * std::pow( aCriteria( iADV ) - warm_restart::target( tVariable ), 2 );
            }
            return { { tObjective } };
        }

        Matrix< DDRMat >
        compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
        {
            Matrix< DDRMat > tGradient( 1, aCriteria.size(), 0.0 );
            for ( uint iADV = 0; iADV < aCriteria.size(); iADV++ )
            {
                const uint tVariable = gSolveOrder[ iADV ] - 100;
                tGradient( 0, iADV ) = 2.0 * warm_restart::weight( tVariable ) * ( aCriteria( iADV ) - warm_restart::target( tVariable ) );
            }
            return tGradient;
        }

        /** ADVs in the order mIds, with those IDs; the first design is mStart (by position), 0.5 if empty */
        class Interface_Ordered : public warm_restart::Interface_Forced_Restart
        {
          public:
            std::vector< moris_id > mIds;
            std::vector< real >     mStart;

            void
            initialize( Vector< real >& aADVs, Vector< real >& aLowerBounds, Vector< real >& aUpperBounds, Matrix< IdMat >& aIjklIds ) override
            {
                mNumADVs = mIds.size();
                aADVs.resize( mNumADVs );
                aLowerBounds.resize( mNumADVs );
                aUpperBounds.resize( mNumADVs );
                for ( uint iADV = 0; iADV < mNumADVs; iADV++ )
                {
                    aADVs( iADV )        = mStart.empty() ? 0.5 : mStart[ iADV ];
                    aLowerBounds( iADV ) = 0.0;
                    aUpperBounds( iADV ) = 1.0;
                }
                aIjklIds                       = id_matrix( mIds );
                mInitializeOptimizationRestart = false;
            }
        };

        /** position of each ID of aIds: tPosition[ id - 100 ] */
        std::vector< uint >
        positions( const std::vector< moris_id >& aIds )
        {
            std::vector< uint > tPosition( aIds.size(), 0 );
            for ( uint iADV = 0; iADV < aIds.size(); iADV++ )
            {
                tPosition[ aIds[ iADV ] - 100 ] = iADV;
            }
            return tPosition;
        }

        /** aValues (in the order aIds) sorted by ID */
        std::vector< real >
        by_id( const std::vector< real >& aValues, const std::vector< moris_id >& aIds )
        {
            std::vector< real > tById( aValues.size(), 0.0 );
            for ( uint iADV = 0; iADV < aIds.size(); iADV++ )
            {
                tById[ aIds[ iADV ] - 100 ] = aValues[ iADV ];
            }
            return tById;
        }

        /** what one run leaves behind, per-ADV values sorted by ID */
        struct Run
        {
            std::vector< std::vector< real > > mDesigns;          ///< every evaluated design
            std::vector< real >                mStartLowWidth;    ///< state the warm solve started with
            std::vector< real >                mStartUppWidth;
            std::vector< real >                mStartMoveStep;
            Algorithm_MMA::Restart_Counts      mCounts;
            int                                mNumWarmSolves = 0;
            std::string                        mStateFile;    ///< copy of the state file of the last iteration
        };

        /**
         * aMaxIts iterations on the problem in the order aIds, from the design aStartById (empty: 0.5) and the state
         * file aStateFile (empty: none). The state file of the last iteration is kept as aKeepAs.
         */
        Run
        run( const std::vector< moris_id >& aIds,
                const std::vector< real >&  aStartById,
                const std::string&          aStateFile,
                sint                        aMaxIts,
                const std::string&          aKeepAs )
        {
            warm_restart::remove_state_files();
            gSolveOrder = aIds;

            std::shared_ptr< Interface_Ordered > tInterface = std::make_shared< Interface_Ordered >();
            tInterface->mIds                                = aIds;
            for ( uint iADV = 0; iADV < aIds.size() && !aStartById.empty(); iADV++ )
            {
                tInterface->mStart.push_back( aStartById[ aIds[ iADV ] - 100 ] );
            }

            Parameter_List             tProblemParameterList = moris::prm::create_opt_problem_parameter_list();
            std::shared_ptr< Problem > tProblem              = std::make_shared< Problem_User_Defined >(
                    tProblemParameterList,
                    tInterface,
                    &warm_restart::get_constraint_types,
                    &compute_objectives,
                    &warm_restart::compute_constraints,
                    &warm_restart::zero_dadv,
                    &compute_dobjective_dcriteria,
                    &warm_restart::zero_dadv,
                    &warm_restart::compute_dconstraint_dcriteria );
            tProblem->initialize();

            Parameter_List tParameters = warm_restart::algorithm_parameters( true, aMaxIts );
            tParameters.set( "restart_state_file", aStateFile );
            Algorithm_MMA tAlgorithm( tParameters );
            warm_restart::solve_with_restarts( tAlgorithm, tProblem );

            Run tRun;
            for ( const std::vector< real >& tDesign : tInterface->mEvaluated )
            {
                tRun.mDesigns.push_back( by_id( tDesign, aIds ) );
            }
            tRun.mNumWarmSolves = tAlgorithm.get_num_warm_solves();
            if ( tRun.mNumWarmSolves > 0 )
            {
                tRun.mStartLowWidth = by_id( tAlgorithm.get_warm_start_state().mLowWidth, aIds );
                tRun.mStartUppWidth = by_id( tAlgorithm.get_warm_start_state().mUppWidth, aIds );
                tRun.mStartMoveStep = by_id( tAlgorithm.get_warm_start_state().mMoveStep, aIds );
                tRun.mCounts        = tAlgorithm.get_warm_start_counts();
            }

            // the state file of the last iteration
            std::string tLatest;
            uint        tLatestIter = 0;
            for ( const std::filesystem::directory_entry& tEntry : std::filesystem::directory_iterator( "." ) )
            {
                const std::string tName = tEntry.path().filename().string();
                if ( tName.rfind( "MMA_State_Alg_0_Iter_", 0 ) == 0 )
                {
                    uint tIter = (uint)std::stoul( tName.substr( 21 ) );
                    if ( tLatest.empty() || tIter > tLatestIter )
                    {
                        tLatest     = tName;
                        tLatestIter = tIter;
                    }
                }
            }
            REQUIRE( !tLatest.empty() );
            std::filesystem::copy_file( tLatest, aKeepAs, std::filesystem::copy_options::overwrite_existing );
            tRun.mStateFile = aKeepAs;
            warm_restart::remove_state_files();
            return tRun;
        }
    }    // namespace state_ids

    TEST_CASE( "Native MMA state file is mapped to the ADV order by ID", "[optimization][restart_keep_state][restart_state_ids]" )
    {
        using namespace state_ids;

        // each rank writes its own files
        const std::string tFileName    = "MMA_State_ids_test_" + std::to_string( par_rank() ) + ".hdf5";
        const std::string tOldFileName = "MMA_State_old_test_" + std::to_string( par_rank() ) + ".hdf5";
        const uint        tNumVar      = 8;

        const std::vector< double > tXMin( tNumVar, 0.0 );
        const std::vector< double > tXMax( tNumVar, 1.0 );

        // written in the order of 1, 2 or 4 processors, read in each of the other orders
        const std::vector< sint > tNumProcs = { 1, 2, 4 };
        for ( sint tWriteProcs : tNumProcs )
        {
            const std::vector< moris_id > tWriteIds = order( tWriteProcs );
            Algorithm_MMA::write_restart_state( tFileName, make_state( tWriteIds, tWriteProcs, true ), 3 );
            Algorithm_MMA::Restart_State tFile = Algorithm_MMA::read_restart_state( tFileName );

            // the file carries the order
            REQUIRE( tFile.mIjklIds.numel() == tNumVar );
            CHECK( tFile.mNumProcs == tWriteProcs );
            for ( uint iADV = 0; iADV < tNumVar; iADV++ )
            {
                CHECK( tFile.mIjklIds( iADV ) == tWriteIds[ iADV ] );
                CHECK( tFile.mCentre[ iADV ] == centre_of( tWriteIds[ iADV ] ) );
            }

            for ( sint tReadProcs : tNumProcs )
            {
                const std::vector< moris_id > tReadIds = order( tReadProcs );
                const std::vector< sint >     tMap     = Algorithm_MMA::restart_state_id_map(
                        tFile, id_matrix( tReadIds ), tNumVar, tReadProcs, -1, tFileName );

                // the same order needs no map
                if ( tReadProcs == tWriteProcs )
                {
                    CHECK( tMap.empty() );
                    continue;
                }
                REQUIRE( tMap.size() == tNumVar );

                // the restart design is the state's centre, ADV by ADV: with the map every ADV keeps its own state
                const std::vector< double >   tX = design_at_centres( tReadIds );
                std::vector< double >         tLow, tUpp, tMove;
                Algorithm_MMA::Restart_Counts tCounts = Algorithm_MMA::apply_restart_state(
                        tFile, tMap, tX, tXMin, tXMax, 0.25, 0.5, true, 0.07, 1.0e-4, 0.1, tLow, tUpp, tMove );
                CHECK( tCounts.mKept == (int)tNumVar );
                CHECK( tCounts.mMoved == 0 );
                CHECK( tCounts.mFresh == 0 );
                for ( uint iADV = 0; iADV < tNumVar; iADV++ )
                {
                    const moris_id tId = tReadIds[ iADV ];
                    CHECK( tFile.mIjklIds( tMap[ iADV ] ) == tId );
                    CHECK( tLow[ iADV ] == Approx( centre_of( tId ) - low_width_of( tId ) ).epsilon( 1e-14 ) );
                    CHECK( tUpp[ iADV ] == Approx( centre_of( tId ) + upp_width_of( tId ) ).epsilon( 1e-14 ) );
                    CHECK( tMove[ iADV ] == Approx( move_step_of( tId ) ).epsilon( 1e-14 ) );
                }

                // read by position (the defect): ADVs count as moved, the others carry another ADV's state
                tCounts = Algorithm_MMA::apply_restart_state(
                        tFile, {}, tX, tXMin, tXMax, 0.25, 0.5, true, 0.07, 1.0e-4, 0.1, tLow, tUpp, tMove );
                CHECK( tCounts.mMoved > 0 );
            }
        }

        // an ADV without an ID keeps its position; the others are mapped
        {
            std::vector< moris_id > tWriteIds = order( 1 );
            std::vector< moris_id > tReadIds  = { 107, 106, gNoID, 104, 103, 105, 101, 100 };
            tWriteIds[ 2 ]                    = gNoID;
            Algorithm_MMA::write_restart_state( tFileName, make_state( tWriteIds, 1, true ), 3 );
            Algorithm_MMA::Restart_State tFile = Algorithm_MMA::read_restart_state( tFileName );
            const std::vector< sint >    tMap  = Algorithm_MMA::restart_state_id_map( tFile, id_matrix( tReadIds ), tNumVar, 4, -1, tFileName );
            CHECK( tMap == std::vector< sint >{ 7, 6, 2, 4, 3, 5, 1, 0 } );

            // ... but not if the file holds an ID at that position
            tReadIds = { 107, gNoID, 102, 104, 103, 105, 101, 100 };
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tFile, id_matrix( tReadIds ), tNumVar, 4, -1, tFileName ) );
        }

        // IDs that do not identify the ADVs are an error: an ID the file does not hold, an ID twice in the file or in
        // the current design, another number of ADVs
        {
            Algorithm_MMA::Restart_State tFile = make_state( order( 1 ), 1, true );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tFile, id_matrix( { 107, 100, 104, 103, 101, 106, 105, 999 } ), tNumVar, 4, -1, tFileName ) );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tFile, id_matrix( { 107, 100, 104, 103, 101, 106, 105, 105 } ), tNumVar, 4, -1, tFileName ) );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tFile, id_matrix( { 107, 100, 104, 103, 101, 106, 105 } ), 7, 4, -1, tFileName ) );

            Algorithm_MMA::Restart_State tTwice = make_state( { 100, 101, 102, 103, 104, 105, 106, 106 }, 1, true );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tTwice, id_matrix( order( 4 ) ), tNumVar, 4, -1, tFileName ) );

            // the same IDs at the same positions are the same order, whatever they are
            CHECK( Algorithm_MMA::restart_state_id_map( tTwice, tTwice.mIjklIds, tNumVar, 4, -1, tFileName ).empty() );
        }

        // a state without IDs (ADVs that carry none) is read in file order on its processor count, and nowhere else
        {
            Algorithm_MMA::write_restart_state( tFileName, make_state( order( 2 ), 2, false ), 3 );
            Algorithm_MMA::Restart_State tFile = Algorithm_MMA::read_restart_state( tFileName );
            CHECK( tFile.mIjklIds.numel() == 0 );
            CHECK( tFile.mNumProcs == 2 );
            CHECK( Algorithm_MMA::restart_state_id_map( tFile, Matrix< IdMat >(), tNumVar, 2, -1, tFileName ).empty() );
            CHECK( Algorithm_MMA::restart_state_id_map( tFile, id_matrix( order( 2 ) ), tNumVar, 2, -1, tFileName ).empty() );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tFile, Matrix< IdMat >(), tNumVar, 4, -1, tFileName ) );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tFile, id_matrix( order( 4 ) ), tNumVar, 4, -1, tFileName ) );

            // a state with IDs, a problem without: the same rule
            Algorithm_MMA::Restart_State tWithIds = make_state( order( 2 ), 2, true );
            CHECK( Algorithm_MMA::restart_state_id_map( tWithIds, Matrix< IdMat >(), tNumVar, 2, -1, tFileName ).empty() );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tWithIds, Matrix< IdMat >(), tNumVar, 1, -1, tFileName ) );
        }

        // a file written before the order was stored: read in file order (as before) on the processor count of the ADV
        // restart file read with it; an error on another count, or if there is no ADV restart file to tell
        {
            write_old_format_state( tOldFileName, make_state( order( 1 ), 1, false ) );
            Algorithm_MMA::Restart_State tOld = Algorithm_MMA::read_restart_state( tOldFileName );
            CHECK( tOld.mIjklIds.numel() == 0 );
            CHECK( tOld.mNumProcs == -1 );
            CHECK( tOld.mCentre == make_state( order( 1 ), 1, false ).mCentre );
            CHECK( tOld.mMoveStep == make_state( order( 1 ), 1, false ).mMoveStep );
            CHECK( Algorithm_MMA::restart_state_id_map( tOld, id_matrix( order( 1 ) ), tNumVar, 1, 1, tOldFileName ).empty() );
            CHECK( Algorithm_MMA::restart_state_id_map( tOld, Matrix< IdMat >(), tNumVar, 4, 4, tOldFileName ).empty() );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tOld, id_matrix( order( 4 ) ), tNumVar, 4, 1, tOldFileName ) );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tOld, id_matrix( order( 1 ) ), tNumVar, 1, 4, tOldFileName ) );
            CHECK_THROWS( Algorithm_MMA::restart_state_id_map( tOld, id_matrix( order( 1 ) ), tNumVar, 1, -1, tOldFileName ) );
        }

        std::filesystem::remove( tFileName );
        std::filesystem::remove( tOldFileName );
    }

    TEST_CASE( "Native MMA restarts from a state file written in another ADV order", "[optimization][restart_keep_state][restart_state_ids]" )
    {
        // one processor only: the test drives solve() on a single rank (MB-091)
        if ( par_size() > 1 )
        {
            return;
        }
        using namespace state_ids;
        const uint tNumVar = 8;

        // each pair: a run in the writer's order, then the same continuation in the writer's order (the reference) and
        // in the reader's order, both from the writer's state file
        const std::vector< std::pair< sint, sint > > tPairs = { { 1, 2 }, { 1, 4 }, { 2, 1 }, { 4, 1 }, { 4, 2 } };
        for ( const std::pair< sint, sint >& tPair : tPairs )
        {
            const std::vector< moris_id > tWriteIds = order( tPair.first );
            const std::vector< moris_id > tReadIds  = order( tPair.second );

            Run tFirst = run( tWriteIds, {}, "", 5, "MMA_State_ids_first.hdf5" );
            REQUIRE( tFirst.mNumWarmSolves == 0 );
            const std::vector< real > tDesign = tFirst.mDesigns.back();

            // the file holds the writer's IDs and the committed design
            Algorithm_MMA::Restart_State tFile = Algorithm_MMA::read_restart_state( tFirst.mStateFile );
            REQUIRE( tFile.mIjklIds.numel() == tNumVar );
            CHECK( tFile.mNumProcs == par_size() );
            for ( uint iADV = 0; iADV < tNumVar; iADV++ )
            {
                CHECK( tFile.mIjklIds( iADV ) == tWriteIds[ iADV ] );
                CHECK( tFile.mCentre[ iADV ] == tDesign[ tWriteIds[ iADV ] - 100 ] );
            }

            Run tReference = run( tWriteIds, tDesign, tFirst.mStateFile, 4, "MMA_State_ids_reference.hdf5" );
            Run tReordered = run( tReadIds, tDesign, tFirst.mStateFile, 4, "MMA_State_ids_reordered.hdf5" );

            // every ADV starts from its own state
            REQUIRE( tReference.mNumWarmSolves == 1 );
            REQUIRE( tReordered.mNumWarmSolves == 1 );
            CHECK( tReference.mCounts.mKept == (int)tNumVar );
            CHECK( tReordered.mCounts.mKept == (int)tNumVar );
            CHECK( tReordered.mCounts.mMoved == 0 );
            CHECK( tReordered.mCounts.mFresh == 0 );
            const std::vector< uint > tFilePosition = positions( tWriteIds );
            for ( uint iVariable = 0; iVariable < tNumVar; iVariable++ )
            {
                CHECK( tReordered.mStartLowWidth[ iVariable ] == tReference.mStartLowWidth[ iVariable ] );
                CHECK( tReordered.mStartUppWidth[ iVariable ] == tReference.mStartUppWidth[ iVariable ] );
                CHECK( tReordered.mStartMoveStep[ iVariable ] == tReference.mStartMoveStep[ iVariable ] );
                CHECK( tReordered.mStartMoveStep[ iVariable ] == tFile.mMoveStep[ tFilePosition[ iVariable ] ] );
            }

            // ... and the two continuations visit the same designs (the sums of the subproblem run in another order)
            REQUIRE( tReordered.mDesigns.size() == tReference.mDesigns.size() );
            REQUIRE( tReference.mDesigns.size() > 3 );
            real tMaxDifference = 0.0;
            real tMaxChange     = 0.0;
            for ( uint iDesign = 0; iDesign < tReference.mDesigns.size(); iDesign++ )
            {
                for ( uint iVariable = 0; iVariable < tNumVar; iVariable++ )
                {
                    tMaxDifference = std::max( tMaxDifference, std::abs( tReordered.mDesigns[ iDesign ][ iVariable ] - tReference.mDesigns[ iDesign ][ iVariable ] ) );
                    tMaxChange     = std::max( tMaxChange, std::abs( tReference.mDesigns[ iDesign ][ iVariable ] - tDesign[ iVariable ] ) );
                }
            }
            CHECK( tMaxDifference < 1.0e-10 );
            CHECK( tMaxChange > 1.0e-3 );    // the continuation moves the design

            // the reordered run writes its own order
            Algorithm_MMA::Restart_State tWritten = Algorithm_MMA::read_restart_state( tReordered.mStateFile );
            REQUIRE( tWritten.mIjklIds.numel() == tNumVar );
            for ( uint iADV = 0; iADV < tNumVar; iADV++ )
            {
                CHECK( tWritten.mIjklIds( iADV ) == tReadIds[ iADV ] );
            }
        }

        std::filesystem::remove( "MMA_State_ids_first.hdf5" );
        std::filesystem::remove( "MMA_State_ids_reference.hdf5" );
        std::filesystem::remove( "MMA_State_ids_reordered.hdf5" );
    }
}    // namespace moris::opt
