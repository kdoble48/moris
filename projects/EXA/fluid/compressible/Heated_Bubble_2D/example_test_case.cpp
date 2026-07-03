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

#include "cl_Logger.hpp" // MRS/IOS/src

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see EXA_RUNNER_RFC.md)

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_heated_bubble_2d
{

TEST_CASE("Heated_Bubble_2D",
        "[moris],[example],[fluid],[compressible],[Heated_Bubble_2D],[EXA_Heated_Bubble_2D]")
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    // (legacy file-scope values: gInterpolationOrder zero-initialized, gPrintReferenceValues = false)
    gInterpolationOrder   = 0;
    gPrintReferenceValues = false;

    char tString1[] = "";
    char tString2[] = "./Heated_Bubble_2D.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // check
    REQUIRE( tRet ==  0 );
}

}    // namespace exa_heated_bubble_2d
