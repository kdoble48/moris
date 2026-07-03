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

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see doc/internal/EXA_RUNNER_RFC.md)

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_single_phase_hollow_cylinder_static
{

TEST_CASE("Single_Phase_Hollow_Cylinder_Static",
        "[moris],[example],[thermal],[diffusion],[EXA_Single_Phase_Hollow_Cylinder_Static]")
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    // (former load-bearing file-scope initializer: uint gInterpolationOrder = 1;)
    gInterpolationOrder   = 1;
    gPrintReferenceValues = false;

    char tString1[] = "";
    char tString2[] = "./Single_Phase_Hollow_Cylinder_Static.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );
}

}    // namespace exa_single_phase_hollow_cylinder_static
