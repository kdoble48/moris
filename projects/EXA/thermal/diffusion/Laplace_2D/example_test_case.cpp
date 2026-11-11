/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 */

#include <catch.hpp>

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"     // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see doc/internal/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_laplace_2d
{

void
check_results(
        const std::string &aExoFileName,
        uint               aTestCaseIndex,
        bool               aSolvewithPetsc = false )
{
    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file (set verbose to true to get basic mesh information)
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, true );

    // define reference node IDs
    Vector< uint > tReferenceNodeId = { 10, 140 };

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << aTestCaseIndex << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';

        // coordinates of reference point
        moris::print( tExoIO.get_nodal_coordinate( tReferenceNodeId( aTestCaseIndex ) ), "Coordinates of reference point" );

        // time value for reference time step
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << '\n';

        // solution of reference point at reference time step
        std::cout << "Temperature at reference point: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_nodal_field_value( tReferenceNodeId( aTestCaseIndex ), 2, 0 ) << '\n';

        // value of IQI at reference time step
        std::cout << "IQI value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_global_variable( 0, 0 ) << '\n';

        return;
    }

    // define reference values for dimension, number of nodes and number of elements
    Vector< uint > tReferenceNumDims  = { 2, 2 };
    Vector< uint > tReferenceNumNodes = { 126, 25 };
    Vector< uint > tReferenceNumElems = { 52, 4 };

    // check dimension, number of nodes and number of elements
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    MORIS_LOG_INFO( "Check number of dimensions: reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumDims( aTestCaseIndex ),
            tNumDims,
            std::abs( ( tNumDims - tReferenceNumDims( aTestCaseIndex ) ) / tReferenceNumDims( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of nodes:      reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumNodes( aTestCaseIndex ),
            tNumNodes,
            std::abs( ( tNumNodes - tReferenceNumNodes( aTestCaseIndex ) ) / tReferenceNumNodes( aTestCaseIndex ) * 100.0 ) );
    MORIS_LOG_INFO( "Check number of elements:   reference %12d, actual %12d, percent error %12.5e.",
            tReferenceNumElems( aTestCaseIndex ),
            tNumElems,
            std::abs( ( tNumElems - tReferenceNumElems( aTestCaseIndex ) ) / tReferenceNumElems( aTestCaseIndex ) * 100.0 ) );

    CHECK( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
    CHECK( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
    CHECK( tNumElems == tReferenceNumElems( aTestCaseIndex ) );

    real tL2ErrorActual       = tExoIO.get_global_variable( 0, 0 );
    real tDiffusiveFluxActual = tExoIO.get_global_variable( 1, 0 );

    Vector< real > tL2ErrorReference       = { 6.65582e-07, 6.65582e-07 };
    Vector< real > tDiffusiveFluxReference = { 0.0504293, 0.0504293 };
    Vector< real > tEigenValueRefrence     = { 1.371841, 1.371841 };

    real tRelL2Difference   = ( tL2ErrorActual - tL2ErrorReference( aTestCaseIndex ) ) / tL2ErrorReference( aTestCaseIndex );
    real tRelFluxDifference = ( tDiffusiveFluxActual - tDiffusiveFluxReference( aTestCaseIndex ) ) / tDiffusiveFluxReference( aTestCaseIndex );

    MORIS_LOG_INFO( "Check nodal temperature:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tL2ErrorReference( aTestCaseIndex ),
            tL2ErrorActual,
            tRelL2Difference * 100.0 );

    MORIS_LOG_INFO( "Check nodal temperature:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            tDiffusiveFluxReference( aTestCaseIndex ),
            tDiffusiveFluxActual,
            tRelFluxDifference * 100.0 );

    CHECK( tRelL2Difference < 1.0e-6 );
    CHECK( tRelFluxDifference < 1.0e-6 );
    if ( aSolvewithPetsc and aTestCaseIndex == 0 )
    {
        real tEigenValueActual   = tExoIO.get_nodal_field_value( 0, 7, 0 );
        real tRelEigenDifference = ( tEigenValueActual - tEigenValueRefrence( aTestCaseIndex ) ) / tEigenValueRefrence( aTestCaseIndex );
        MORIS_LOG_INFO( "Check eigenvalue:  reference %12.5e, actual %12.5e, percent error %12.5e.",
                tEigenValueRefrence( aTestCaseIndex ),
                tEigenValueActual,
                tRelEigenDifference * 100.0 );

        CHECK( tRelEigenDifference < 1.0e-6 );
    }
}

// free function to check the results fully
void
check_results_serial_parallel( bool aSolvewithPetsc = false )
{
    // check results
    switch ( par_size() )
    {
        // Test Case 0
        case 1:
        {
            // perform check
            check_results( "Laplace.exo.e-s.0000", 0, aSolvewithPetsc );
            break;
        }
        // Test Case 1
        case 4:
        {
            if ( par_rank() == 1 )
            {
                // set screen output processor
                gLogger.set_screen_output_rank( 1 );

                // perform check
                check_results( "Laplace.exo.e-s.0000", 1, aSolvewithPetsc );

                // reset screen output processor
                gLogger.set_screen_output_rank( 0 );
            }
            break;
        }
        default:
        {
            MORIS_ERROR( false, "Example problem not configured for %d processors.", par_size() );
        }
    }
}

TEST_CASE( "Laplace_Anasazi",
        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Anasazi],[EXA_Laplace_2D]" )
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    gTestIndex            = 0;
    gInterpolationOrder   = 1;
    gPrintReferenceValues = false;

    char tString1[] = "";
    char tString2[] = "./Laplace.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check
    REQUIRE( tRet == 0 );

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing HeatConduction: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    check_results_serial_parallel( false );
}

#ifdef MORIS_HAVE_SLEPC
TEST_CASE( "Laplace_Slepc",
        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Slepc],[EXA_Laplace_2D]" )
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    gTestIndex            = 1;
    gInterpolationOrder   = 1;
    gPrintReferenceValues = false;

    char tString1[] = "";
    char tString2[] = "./Laplace.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check
    REQUIRE( tRet == 0 );

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing HeatConduction: Interpolation order 1 - %i Processors.", par_size() );
    MORIS_LOG_INFO( " " );

    check_results_serial_parallel( true );
}
#endif

//---------------------------------------------------------------
// Preconditioner cases: the same temperature problem as above, solved with the preconditioners that could not be
// used before MB-102 (Ifpack point relaxation and Krylov) and MB-103 (PETSc block inverse, PCType "mat"), and with
// Ifpack ILU for comparison. Laplace.cpp describes the cases and when each preconditioner is the right choice.

// check the solution against the direct solve of the same problem, from both sides
void check_matches_direct_solve( const std::string &aExoFileName )
{
    moris::mtk::Exodus_IO_Helper tExoIO( aExoFileName.c_str(), 0, false, true );

    real tL2ErrorActual       = tExoIO.get_global_variable( 0, 0 );
    real tDiffusiveFluxActual = tExoIO.get_global_variable( 1, 0 );

    // values of the Laplace_Anasazi case, whose one-level ML preconditioner is a direct solve (same on 1 and 4 processors)
    real tL2ErrorReference       = 6.65582452961e-07;
    real tDiffusiveFluxReference = 5.04292797169e-02;

    real tRelL2Difference   = std::abs( tL2ErrorActual - tL2ErrorReference ) / tL2ErrorReference;
    real tRelFluxDifference = std::abs( tDiffusiveFluxActual - tDiffusiveFluxReference ) / tDiffusiveFluxReference;

    MORIS_LOG_INFO( "Check L2 error against direct solve:       reference %19.12e, actual %19.12e, relative difference %12.5e.",
            tL2ErrorReference,
            tL2ErrorActual,
            tRelL2Difference );
    MORIS_LOG_INFO( "Check diffusive flux against direct solve: reference %19.12e, actual %19.12e, relative difference %12.5e.",
            tDiffusiveFluxReference,
            tDiffusiveFluxActual,
            tRelFluxDifference );

    CHECK( tRelL2Difference < 1.0e-6 );
    CHECK( tRelFluxDifference < 1.0e-6 );
}

// run the deck with the linear solver case aTestIndex and check the results
void run_preconditioner_case(
        uint aTestIndex,
        bool aSolvewithPetsc )
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    gTestIndex            = aTestIndex;
    gInterpolationOrder   = 1;
    gPrintReferenceValues = false;

    char tString1[] = "";
    char tString2[] = "./Laplace.so";

    char *argv[ 2 ] = { tString1, tString2 };

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check
    REQUIRE( tRet == 0 );

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing HeatConduction: linear solver case %i - %i Processors.", aTestIndex, par_size() );
    MORIS_LOG_INFO( " " );

    check_results_serial_parallel( aSolvewithPetsc );

    // the output is written and read by the same processor as in check_results_serial_parallel
    if ( par_size() == 1 or par_rank() == 1 )
    {
        check_matches_direct_solve( "Laplace.exo.e-s.0000" );
    }
}

TEST_CASE( "Laplace_Belos_Ifpack_ILU",
        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Preconditioners],[EXA_Laplace_2D]" )
{
    run_preconditioner_case( 2, false );
}

TEST_CASE( "Laplace_Belos_Ifpack_Jacobi",
        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Preconditioners],[EXA_Laplace_2D]" )
{
    run_preconditioner_case( 3, false );
}

TEST_CASE( "Laplace_Belos_Ifpack_Symmetric_Gauss_Seidel",
        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Preconditioners],[EXA_Laplace_2D]" )
{
    run_preconditioner_case( 4, false );
}

TEST_CASE( "Laplace_Belos_Ifpack_Krylov",
        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Preconditioners],[EXA_Laplace_2D]" )
{
    run_preconditioner_case( 5, false );
}

#ifdef MORIS_HAVE_SLEPC
// the PETSc deck also sets up the SLEPc eigenvalue solves, hence the same guard as Laplace_Slepc
TEST_CASE( "Laplace_Petsc_Block_Inverse",
        "[moris],[example],[thermal],[Laplace_2D],[Laplace_Preconditioners],[EXA_Laplace_2D]" )
{
    // serial only: the preconditioner reads the rows of every block from the system matrix, and in parallel a block
    // can contain rows of another processor; MORIS stops with an error if PCType "mat" is built on more than one
    // processor (MB-111)
    if ( par_size() == 1 )
    {
        run_preconditioner_case( 6, true );
    }
}
#endif

}    // namespace exa_laplace_2d
