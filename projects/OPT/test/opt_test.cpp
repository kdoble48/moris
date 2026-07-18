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

#include <fstream>
#include <cstdio>

#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_OPT_Manager.hpp"
#include "fn_OPT_create_interface.hpp"
#include "cl_OPT_Problem_User_Defined.hpp"
#include "cl_OPT_Interface_User_Defined.hpp"
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
