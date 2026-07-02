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

namespace exa_channel_with_four_cylinders_transient
{

TEST_CASE("Channel_with_Four_Cylinders_Transient",
        "[moris],[example],[thermal],[advection],[EXA_Channel_with_Four_Cylinders_Transient]")
{
    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Channel_with_Four_Cylinders_Transient.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    REQUIRE( tRet ==  0 );
}

}    // namespace exa_channel_with_four_cylinders_transient
