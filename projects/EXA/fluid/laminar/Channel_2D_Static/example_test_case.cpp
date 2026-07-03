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

namespace exa_channel_2d_static
{

TEST_CASE("Channel_2D_Static_Inlet_Velocity",
        "[moris],[example],[fluid],[laminar],[EXA_Channel_2D_Static]")
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    // (legacy file-scope initializers: gInletVelocityBCFlag = true, gInletPressureBCFlag = false)
    gInletVelocityBCFlag = true;
    gInletPressureBCFlag = false;

    char tString1[] = "";
    char tString2[] = "./Channel_2D_Static.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    //REQUIRE( tRet ==  0 );
}

//---------------------------------------------------------------

TEST_CASE("Channel_2D_Static_Inlet_Pressure",
        "[moris],[example],[fluid],[laminar],[EXA_Channel_2D_Static]")
{
    // define command line call
    int argc = 2;

    // set all globals this test case or its deck consumes (rule R4)
    // NOTE: the legacy TU never toggled these flags between cases - despite its name,
    // this case ran with the same file-scope values as the velocity case
    // (gInletVelocityBCFlag = true, gInletPressureBCFlag = false); preserved verbatim.
    gInletVelocityBCFlag = true;
    gInletPressureBCFlag = false;

    char tString1[] = "";
    char tString2[] = "./Channel_2D_Static.so";

    char * argv[2] = {tString1,tString2};

    // call to performance manager main interface
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // catch test statements should follow
    //REQUIRE( tRet ==  0 );
}

}    // namespace exa_channel_2d_static
