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

#include "cl_Logger.hpp"                // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"  // MTK/src
#include "cl_Communication_Tools.hpp"   // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see share/doc/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_couette_flow_static
{

TEST_CASE("Couette_Flow_Static",
        "[moris],[example],[fluid],[laminar],[EXA_Couette_Flow_Static]")
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    gPrintReferenceValues = false;

    char tString1[] = "";
    char tString2[] = "./Couette_Flow_Static.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );
}

}    // namespace exa_couette_flow_static
