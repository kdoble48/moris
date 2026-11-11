/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Linear_Solver_Status_Test.cpp
 *
 * A failed or unconverged linear solve must reach the linear solver manager: the algorithms return a non-zero status
 * (MB-101), the manager restarts with the configured fallback algorithms without reading past its list (MB-100), and
 * honours DLA_hard_break and DLA_rebuild_lin_solver_on_fail.
 */

#include <memory>
#include <string>

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_DLA_Linear_Solver.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Preconditioner_Trilinos.hpp"

#ifdef MORIS_HAVE_PETSC
#include "cl_DLA_Preconditioner_PETSc.hpp"
#endif

#include "cl_Solver_Interface_Proxy.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

namespace moris::dla
{
    namespace
    {
        /**
         * A linear solver algorithm that returns a fixed status and records how it was called.
         */
        class Status_Solver_Algorithm : public Linear_Solver_Algorithm
        {
          public:
            moris::sint mStatus;
            moris::uint mNumCalls = 0;
            moris::sint mLastIter = 0;

            explicit Status_Solver_Algorithm( const moris::sint aStatus )
                    : Linear_Solver_Algorithm( Parameter_List( "Status Solver" ) )
                    , mStatus( aStatus )
            {
            }

            moris::sint
            solve_linear_system() override
            {
                return mStatus;
            }

            moris::sint
            solve_linear_system( Linear_Problem* aLinearSystem, const moris::sint aIter ) override
            {
                mNumCalls++;
                mLastIter = aIter;
                return mStatus;
            }
        };

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Parameters of an Ifpack ILU(0) preconditioner. The preconditioner keeps a reference to the list.
         */
        Parameter_List
        create_ilu_parameter_list()
        {
            Parameter_List tParamList = prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK );
            tParamList.set( "ifpack_prec_type", std::string( "ILU" ) );
            tParamList.set( "fact: level-of-fill", 0 );
            tParamList.set( "overlap-level", 0 );
            return tParamList;
        }

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Belos GMRES capped at aMaxIterations iterations, tolerance 1e-14.
         */
        Parameter_List
        create_capped_belos_parameter_list( const moris::sint aMaxIterations )
        {
            Parameter_List tParamList = prm::create_linear_algorithm_parameter_list_belos();
            tParamList.set( "Maximum Iterations", aMaxIterations );
            tParamList.set( "Convergence Tolerance", 1e-14 );
            return tParamList;
        }
    }    // namespace

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Linear solver restarts and hard break (MB-100, MB-101)", "[DLA_Linear_Solver_Status],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( 1 );

            Solver_Factory tSolFactory;

            Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );

            Parameter_List tManagerParameterList = prm::create_linear_solver_parameter_list();

            std::shared_ptr< Status_Solver_Algorithm > tFailing    = std::make_shared< Status_Solver_Algorithm >( 1 );
            std::shared_ptr< Status_Solver_Algorithm > tSucceeding = std::make_shared< Status_Solver_Algorithm >( 0 );

            SECTION( "one algorithm, one restart: the restart reuses it and does not read past the list" )
            {
                tManagerParameterList.set( "DLA_max_lin_solver_restarts", 1 );
                tManagerParameterList.set( "DLA_hard_break", false );

                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 0, tFailing );

                tLinSolver.solver_linear_system( tLinProblem, 3 );

                CHECK( tFailing->mNumCalls == 2 );
                CHECK( tFailing->mLastIter == 3 );
            }

            SECTION( "fallback list 0,1 and one restart: the restart runs algorithm 1" )
            {
                tManagerParameterList.set( "DLA_max_lin_solver_restarts", 1 );

                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 0, tFailing );
                tLinSolver.set_linear_algorithm( 1, tSucceeding );

                // the fallback succeeds, so the hard break (default true) does not fire
                CHECK_NOTHROW( tLinSolver.solver_linear_system( tLinProblem, 1 ) );

                CHECK( tFailing->mNumCalls == 1 );
                CHECK( tSucceeding->mNumCalls == 1 );
            }

            SECTION( "fallback list 0,1 and three restarts: the restarts after the list stay on its last entry" )
            {
                tManagerParameterList.set( "DLA_max_lin_solver_restarts", 3 );
                tManagerParameterList.set( "DLA_hard_break", false );

                std::shared_ptr< Status_Solver_Algorithm > tAlsoFailing = std::make_shared< Status_Solver_Algorithm >( 1 );

                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 0, tFailing );
                tLinSolver.set_linear_algorithm( 1, tAlsoFailing );

                tLinSolver.solver_linear_system( tLinProblem, 1 );

                CHECK( tFailing->mNumCalls == 1 );
                CHECK( tAlsoFailing->mNumCalls == 3 );
            }

            SECTION( "setting entry 0 after entry 1 keeps the fallback" )
            {
                tManagerParameterList.set( "DLA_max_lin_solver_restarts", 1 );

                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 1, tSucceeding );
                tLinSolver.set_linear_algorithm( 0, tFailing );

                CHECK_NOTHROW( tLinSolver.solver_linear_system( tLinProblem, 1 ) );

                CHECK( tFailing->mNumCalls == 1 );
                CHECK( tSucceeding->mNumCalls == 1 );
            }

            SECTION( "no restarts: the fallback is not used" )
            {
                tManagerParameterList.set( "DLA_hard_break", false );

                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 0, tFailing );
                tLinSolver.set_linear_algorithm( 1, tSucceeding );

                tLinSolver.solver_linear_system( tLinProblem, 1 );

                CHECK( tFailing->mNumCalls == 1 );
                CHECK( tSucceeding->mNumCalls == 0 );
            }

            SECTION( "DLA_hard_break (default true) aborts on a failed solve" )
            {
                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 0, tFailing );

                CHECK_THROWS( tLinSolver.solver_linear_system( tLinProblem, 1 ) );
            }

            SECTION( "DLA_hard_break = false continues after a failed solve" )
            {
                tManagerParameterList.set( "DLA_hard_break", false );

                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 0, tFailing );

                CHECK_NOTHROW( tLinSolver.solver_linear_system( tLinProblem, 1 ) );
            }

            SECTION( "DLA_rebuild_lin_solver_on_fail: the restart rebuilds the preconditioner (iteration index 1)" )
            {
                tManagerParameterList.set( "DLA_max_lin_solver_restarts", 1 );
                tManagerParameterList.set( "DLA_rebuild_lin_solver_on_fail", true );

                Linear_Solver tLinSolver( tManagerParameterList );
                tLinSolver.set_linear_algorithm( 0, tFailing );
                tLinSolver.set_linear_algorithm( 1, tSucceeding );

                tLinSolver.solver_linear_system( tLinProblem, 4 );

                CHECK( tFailing->mLastIter == 4 );
                CHECK( tSucceeding->mLastIter == 1 );
            }

            delete tLinProblem;
            delete tSolverInterface;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Belos reports an unconverged solve (MB-101)", "[DLA_Linear_Solver_Status],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( 1 );

            Solver_Factory tSolFactory;

            Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );

            tLinProblem->assemble_jacobian();
            tLinProblem->assemble_residual();

            Parameter_List          tPrecParameterList = create_ilu_parameter_list();
            Preconditioner_Trilinos tPreconditioner( tPrecParameterList );

            SECTION( "capped at one iteration: non-zero status" )
            {
                std::shared_ptr< Linear_Solver_Algorithm > tBelos = tSolFactory.create_solver( create_capped_belos_parameter_list( 1 ) );
                tBelos->set_preconditioner( &tPreconditioner );

                CHECK( tBelos->solve_linear_system( tLinProblem ) != 0 );
            }

            SECTION( "default iteration limit: converges, status 0" )
            {
                Parameter_List tBelosParameterList = prm::create_linear_algorithm_parameter_list_belos();

                std::shared_ptr< Linear_Solver_Algorithm > tBelos = tSolFactory.create_solver( tBelosParameterList );
                tBelos->set_preconditioner( &tPreconditioner );

                CHECK( tBelos->solve_linear_system( tLinProblem ) == 0 );

                Matrix< DDRMat > tSol;
                tLinProblem->get_solution( tSol );
                CHECK( equal_to( tSol( 5, 0 ), -0.0138889, 1.0e+08 ) );
            }

            SECTION( "capped at one iteration under the manager: hard break, or the fallback" )
            {
                std::shared_ptr< Linear_Solver_Algorithm > tBelosCapped = tSolFactory.create_solver( create_capped_belos_parameter_list( 1 ) );
                tBelosCapped->set_preconditioner( &tPreconditioner );

                Parameter_List tManagerParameterList = prm::create_linear_solver_parameter_list();

                // without a fallback the default DLA_hard_break = true aborts
                Linear_Solver tLinSolverNoFallback( tManagerParameterList );
                tLinSolverNoFallback.set_linear_algorithm( 0, tBelosCapped );
                CHECK_THROWS( tLinSolverNoFallback.solver_linear_system( tLinProblem, 1 ) );

                // with an uncapped Belos as fallback the restart converges
                Parameter_List tBelosParameterList = prm::create_linear_algorithm_parameter_list_belos();

                std::shared_ptr< Linear_Solver_Algorithm > tBelos = tSolFactory.create_solver( tBelosParameterList );
                tBelos->set_preconditioner( &tPreconditioner );

                tManagerParameterList.set( "DLA_max_lin_solver_restarts", 1 );
                Linear_Solver tLinSolverFallback( tManagerParameterList );
                tLinSolverFallback.set_linear_algorithm( 0, tBelosCapped );
                tLinSolverFallback.set_linear_algorithm( 1, tBelos );
                CHECK_NOTHROW( tLinSolverFallback.solver_linear_system( tLinProblem, 1 ) );

                Matrix< DDRMat > tSol;
                tLinProblem->get_solution( tSol );
                CHECK( equal_to( tSol( 5, 0 ), -0.0138889, 1.0e+08 ) );
            }

            delete tLinProblem;
            delete tSolverInterface;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "ML stand-alone reports an unconverged solve (MB-101)", "[DLA_Linear_Solver_Status],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( 1 );

            Solver_Factory tSolFactory;

            Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );

            tLinProblem->assemble_jacobian();
            tLinProblem->assemble_residual();

            // a one-level hierarchy whose "coarse solve" is a single Jacobi sweep: one cycle cannot reach 1e-12
            Parameter_List tPrecParameterList = prm::create_preconditioner_parameter_list( sol::PreconditionerType::ML );
            tPrecParameterList.set( "ml_prec_type", std::string( "SA" ) );
            tPrecParameterList.set( "max levels", 1 );
            tPrecParameterList.set( "coarse: type", std::string( "Jacobi" ) );
            tPrecParameterList.set( "coarse: sweeps", 1 );

            Preconditioner_Trilinos tPreconditioner( tPrecParameterList );

            Parameter_List tMLParameterList = prm::create_linear_algorithm_parameter_list( sol::SolverType::ML );
            tMLParameterList.set( "Max_Iter", 1 );
            tMLParameterList.set( "Convergence_Tolerance", 1e-12 );

            std::shared_ptr< Linear_Solver_Algorithm > tML = tSolFactory.create_solver( tMLParameterList );
            tML->set_preconditioner( &tPreconditioner );

            CHECK( tML->solve_linear_system( tLinProblem ) != 0 );

            delete tLinProblem;
            delete tSolverInterface;
        }
    }

    //------------------------------------------------------------------------------------------------------------------

#ifdef MORIS_HAVE_PETSC
    TEST_CASE( "PETSc KSP reports an unconverged solve (MB-101)", "[DLA_Linear_Solver_Status],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( 1 );

            Solver_Factory tSolFactory;

            Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Petsc, true );

            tLinProblem->assemble_residual_and_jacobian();

            // the preconditioner keeps a reference to its parameter list
            Parameter_List tPrecParameterList = prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC );
            tPrecParameterList.set( "PCType", std::string( "none" ) );
            Preconditioner_PETSc tPreconditioner( tPrecParameterList );

            SECTION( "capped at one iteration: non-zero status" )
            {
                Parameter_List tKSPParameterList = prm::create_linear_algorithm_parameter_list_petsc();
                tKSPParameterList.set( "KSPType", std::string( "gmres" ) );
                tKSPParameterList.set( "KSPMaxits", 1 );
                tKSPParameterList.set( "KSPTol", 1e-14 );

                std::shared_ptr< Linear_Solver_Algorithm > tKSP = tSolFactory.create_solver( tKSPParameterList );
                tKSP->set_preconditioner( &tPreconditioner );

                CHECK( tKSP->solve_linear_system( tLinProblem ) != 0 );
            }

            SECTION( "default iteration limit: converges, status 0" )
            {
                Parameter_List tKSPParameterList = prm::create_linear_algorithm_parameter_list_petsc();
                tKSPParameterList.set( "KSPType", std::string( "gmres" ) );

                std::shared_ptr< Linear_Solver_Algorithm > tKSP = tSolFactory.create_solver( tKSPParameterList );
                tKSP->set_preconditioner( &tPreconditioner );

                CHECK( tKSP->solve_linear_system( tLinProblem ) == 0 );
            }

            delete tLinProblem;
            delete tSolverInterface;
        }
    }
#endif

}    // namespace moris::dla
