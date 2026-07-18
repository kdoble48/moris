/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 * Twin of LevelSet_Boxbeam/example_test_case.cpp for the Deck API v2 port:
 * the same gate (the GCMMA run completes with exit 0).
 */

#include <catch.hpp>

#include "cl_Logger.hpp"

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see doc/internal/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_levelset_boxbeam_newio
{

    TEST_CASE( "LevelSet_Boxbeam_NewIO",
            "[moris],[example],[optimization],[levelset_boxbeam],[EXA_LevelSet_Boxbeam_NewIO]" )
    {
        // define command line call
        int argc = 2;

        char tString1[] = "";
        char tString2[] = "LevelSet_Boxbeam_NewIO.so";

        char* argv[ 2 ] = { tString1, tString2 };

        // call to performance manager main interface
        int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

        // catch test statements should follow
        REQUIRE( tRet == 0 );
    }

}    // namespace exa_levelset_boxbeam_newio
