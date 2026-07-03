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

#include "cl_Logger.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"
#include "HDF5_Tools.hpp"

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see doc/internal/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_solver_examples
{

TEST_CASE("Standard_Monolithic",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Standard_Monolithic],[EXA_Solver_Examples]")
{
    // set all globals this test case or its deck consumes (rule R4)
    gHaveStaggeredFA = false;
    gHaveStaggeredSA = false;
    gUseMixedTimeElements = false;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

// FIXME: This solver configuration returns an error in the matrix assembly due to the mixed time elements
// TEST_CASE("Monolithic_Mixed_Time_Elements",
//         "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Monolithic_Mixed_Time_Elements],[EXA_Solver_Examples]")
// {
//     // set all globals this test case or its deck consumes (rule R4)
//     gHaveStaggeredFA = false;
//     gHaveStaggeredSA = false;
//     gUseMixedTimeElements = true;
//     gUseBelosWithILUT = false;
//
//     // define command line call
//     int argc = 2;
//
//     char tString1[] = "";
//     char tString2[] = "Solver_Examples.so";
//
//     char * argv[2] = {tString1,tString2};
//
//     // call to performance manager main interface
//     int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );
//
//     // catch test statements should follow
//     REQUIRE( tRet ==  0 );
// }

//---------------------------------------------------------------

TEST_CASE("Staggered_FA_Monolithic_SA",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Staggered_FA_Monolithic_SA],[EXA_Solver_Examples]")
{
    // set all globals this test case or its deck consumes (rule R4)
    gHaveStaggeredFA = true;
    gHaveStaggeredSA = false;
    gUseMixedTimeElements = false;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

TEST_CASE("Staggered_FA_and_SA",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Staggered_FA_and_SA],[EXA_Solver_Examples]")
{
    // set all globals this test case or its deck consumes (rule R4)
    gHaveStaggeredFA = true;
    gHaveStaggeredSA = true;
    gUseMixedTimeElements = false;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

TEST_CASE("Staggered_FA_and_SA_Mixed_Time_Elements",
        "[moris],[example],[optimization],[Solver_Examples_Thermo_Elastic],[Staggered_FA_and_SA_Mixed_Time_Elements],[EXA_Solver_Examples]")
{
    // set all globals this test case or its deck consumes (rule R4)
    gHaveStaggeredFA = true;
    gHaveStaggeredSA = true;
    gUseMixedTimeElements = true;
    gUseBelosWithILUT = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Solver_Examples.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

}    // namespace exa_solver_examples

