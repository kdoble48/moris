/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_DLA_Preconditioner_Test.cpp
 *
 */

#include <string>
#include <variant>

#include "catch.hpp"
#include "fn_equal_to.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

#include "cl_Communication_Tools.hpp"

#include "cl_DLA_Solver_Factory.hpp"
#include "cl_DLA_Linear_Problem.hpp"
#include "cl_DLA_Linear_Solver_Algorithm.hpp"
#include "cl_DLA_Preconditioner_Trilinos.hpp"
#include "cl_SOL_Dist_Matrix.hpp"
#include "cl_SOL_Dist_Vector.hpp"

#ifdef MORIS_HAVE_PETSC
#include <petscksp.h>
#include "cl_DLA_Preconditioner_PETSc.hpp"
#endif

#include "cl_Solver_Interface_Proxy.hpp"

#include "fn_PRM_SOL_Parameters.hpp"

namespace moris::dla
{
    namespace
    {
        /**
         * Parameters of the Ifpack preconditioner aIfpackType, all other parameters at their PRM defaults.
         */
        Parameter_List
        create_ifpack_parameter_list( const std::string& aIfpackType )
        {
            Parameter_List tPrecParameterList = prm::create_preconditioner_parameter_list( sol::PreconditionerType::IFPACK );
            tPrecParameterList.set( "ifpack_prec_type", aIfpackType );
            return tPrecParameterList;
        }

        //--------------------------------------------------------------------------------------------------------------

        /**
         * Solve the serial proxy problem with aNumRHS right-hand sides with Belos GMRES, preconditioned by the Ifpack
         * preconditioner described by aPrecParameterList, and check the solution.
         */
        void
        solve_proxy_with_ifpack_preconditioner(
                const Parameter_List& aPrecParameterList,
                const uint            aNumRHS )
        {
            Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( aNumRHS );

            Solver_Factory tSolFactory;

            Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Epetra );

            Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_belos();

            std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( tLinearSolverParameterList );

            // the preconditioner keeps a reference to its parameter list, which outlives it here
            Preconditioner_Trilinos tPreconditioner( aPrecParameterList );
            tLinSolver->set_preconditioner( &tPreconditioner );

            tLinProblem->assemble_jacobian();
            tLinProblem->assemble_residual();

            // building the preconditioner used to throw std::bad_variant_access (MB-102)
            moris::sint tStatus = tLinSolver->solve_linear_system( tLinProblem );
            CHECK( tStatus == 0 );

            Matrix< DDRMat > tSol;
            tLinProblem->get_solution( tSol );

            for ( uint iRHS = 0; iRHS < aNumRHS; iRHS++ )
            {
                CHECK( equal_to( tSol( 5, iRHS ), -0.0138889, 1.0e+08 ) );
                CHECK( equal_to( tSol( 12, iRHS ), -0.00694444, 1.0e+08 ) );
            }

            delete tSolverInterface;
            delete tLinProblem;
        }
    }    // namespace

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Ifpack point relaxation preconditioner (MB-102)", "[DLA_Ifpack_Parameter_Types],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            solve_proxy_with_ifpack_preconditioner( create_ifpack_parameter_list( "point relaxation" ), 2 );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Ifpack Krylov preconditioner (MB-102)", "[DLA_Ifpack_Parameter_Types],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            // an inner solve to round-off, so that the preconditioner is a fixed linear operator for the (non-flexible)
            // outer GMRES; this also checks that the integer and real Krylov parameters reach Ifpack
            Parameter_List tPrecParameterList = create_ifpack_parameter_list( "Krylov" );
            tPrecParameterList.set( "krylov: iterations", 100 );
            tPrecParameterList.set( "krylov: tolerance", 1e-14 );

            // one right-hand side: the inner AztecOO solve of Ifpack_Krylov handles one vector at a time
            solve_proxy_with_ifpack_preconditioner( tPrecParameterList, 1 );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Ifpack SPARSKIT preconditioner parameters (MB-102)", "[DLA_Ifpack_Parameter_Types],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            // Ifpack is often built without SPARSKIT; then the build must fail with a MORIS error naming the type,
            // never with std::bad_variant_access from reading a SPARSKIT parameter with the wrong type
            bool tBadVariantAccess = false;
            try
            {
                solve_proxy_with_ifpack_preconditioner( create_ifpack_parameter_list( "SPARSKIT" ), 1 );
            } catch ( const std::bad_variant_access& )
            {
                tBadVariantAccess = true;
            } catch ( ... )
            {
                // Ifpack without SPARSKIT support
            }
            CHECK( not tBadVariantAccess );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

#ifdef MORIS_HAVE_PETSC
    TEST_CASE( "PETSc mat preconditioner applies the block inverse (MB-103)", "[DLA_PETSc_Mat_Preconditioner],[Linear Solver],[DistLinAlg]" )
    {
        if ( par_size() == 1 )
        {
            Solver_Interface* tSolverInterface = new Solver_Interface_Proxy( 1 );

            Solver_Factory tSolFactory;

            Linear_Problem* tLinProblem = tSolFactory.create_linear_system( tSolverInterface, sol::MapType::Petsc, true );

            tLinProblem->assemble_residual_and_jacobian();

            Mat tSystemMat = tLinProblem->get_matrix()->get_petsc_matrix();

            // the preconditioner keeps a reference to its parameter list
            Parameter_List tPrecParameterList = prm::create_preconditioner_parameter_list( sol::PreconditionerType::PETSC );
            tPrecParameterList.set( "PCType", std::string( "mat" ) );

            Preconditioner_PETSc tPreconditioner( tPrecParameterList );

            SECTION( "the PC applies the block inverse, not the system matrix" )
            {
                // build twice: the second build replaces the first block inverse
                for ( uint iBuild = 0; iBuild < 2; iBuild++ )
                {
                    KSP tKSP;
                    KSPCreate( PETSC_COMM_WORLD, &tKSP );
                    KSPSetOperators( tKSP, tSystemMat, tSystemMat );

                    tPreconditioner.build_preconditioner( tLinProblem, tKSP );

                    PC tPC;
                    KSPGetPC( tKSP, &tPC );

                    PCType tPCType;
                    PCGetType( tPC, &tPCType );
                    CHECK( std::string( tPCType ) == std::string( PCMAT ) );

                    // the PC's Pmat must be the block inverse (MB-103: it was the system matrix itself)
                    Mat tAmat;
                    Mat tPmat;
                    PCGetOperators( tPC, &tAmat, &tPmat );
                    CHECK( tAmat == tSystemMat );
                    REQUIRE( tPmat != tSystemMat );

                    // the proxy has no bulk sets, so every block is one dof: Pmat = diag( A )^-1
                    Vec tDiagA;
                    Vec tDiagP;
                    MatCreateVecs( tSystemMat, &tDiagA, nullptr );
                    MatCreateVecs( tSystemMat, &tDiagP, nullptr );
                    MatGetDiagonal( tSystemMat, tDiagA );
                    MatGetDiagonal( tPmat, tDiagP );

                    PetscInt tFirstRow = 0;
                    PetscInt tEndRow   = 0;
                    VecGetOwnershipRange( tDiagA, &tFirstRow, &tEndRow );

                    const PetscScalar* tDiagAValues = nullptr;
                    const PetscScalar* tDiagPValues = nullptr;
                    VecGetArrayRead( tDiagA, &tDiagAValues );
                    VecGetArrayRead( tDiagP, &tDiagPValues );

                    for ( PetscInt iRow = 0; iRow < tEndRow - tFirstRow; iRow++ )
                    {
                        if ( tDiagAValues[ iRow ] != 0.0 )
                        {
                            CHECK( equal_to( tDiagPValues[ iRow ] * tDiagAValues[ iRow ], 1.0 ) );
                        }
                    }

                    VecRestoreArrayRead( tDiagA, &tDiagAValues );
                    VecRestoreArrayRead( tDiagP, &tDiagPValues );
                    VecDestroy( &tDiagA );
                    VecDestroy( &tDiagP );

                    KSPDestroy( &tKSP );
                }
            }

            SECTION( "GMRES preconditioned with the block inverse solves the system" )
            {
                Parameter_List tLinearSolverParameterList = prm::create_linear_algorithm_parameter_list_petsc();
                tLinearSolverParameterList.set( "KSPType", std::string( "gmres" ) );
                tLinearSolverParameterList.set( "PCType", std::string( "mat" ) );

                std::shared_ptr< Linear_Solver_Algorithm > tLinSolver = tSolFactory.create_solver( tLinearSolverParameterList );
                tLinSolver->set_preconditioner( &tPreconditioner );

                // two solves: the block inverse of the first solve is released when the second one builds its own
                for ( uint iSolve = 0; iSolve < 2; iSolve++ )
                {
                    tLinSolver->solve_linear_system( tLinProblem );

                    // solution vector includes the 3 constrained dofs
                    Matrix< DDRMat > tSol;
                    tLinProblem->get_solution( tSol );

                    CHECK( equal_to( tSol( 5 + 3, 0 ), -0.0138889, 1.0e+08 ) );
                    CHECK( equal_to( tSol( 12 + 3, 0 ), -0.00694444, 1.0e+08 ) );
                }
            }

            delete tLinProblem;
            delete tSolverInterface;
        }
    }
#endif

}    // namespace moris::dla
