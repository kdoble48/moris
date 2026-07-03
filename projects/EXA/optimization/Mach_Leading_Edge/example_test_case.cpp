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

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see share/doc/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char * argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_mach_leading_edge
{

TEST_CASE("Mach_Leading_Edge",
        "[moris],[example],[optimization],[Mach_Leading_Edge],[EXA_Mach_Leading_Edge]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "Mach_Leading_Edge.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

}    // namespace exa_mach_leading_edge

