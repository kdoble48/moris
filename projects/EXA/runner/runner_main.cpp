/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * runner_main.cpp
 *
 * Single entry point for the consolidated EXA-test runner (see EXA_RUNNER_RFC.md).
 * Each ctest entry invokes this binary with a Catch2 tag spec ("[EXA_<Example>]")
 * from that example's working directory, where its input .so lives.
 */

#define CATCH_CONFIG_RUNNER
#include <catch.hpp>
#include "cl_Communication_Manager.hpp"   // COM/src
#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_Performance_Reporter.hpp"    // MRS/IOS/src
#include "banner.hpp"                     // COR/src

moris::Comm_Manager gMorisComm;
moris::Logger       gLogger;

//---------------------------------------------------------------

int main( int argc, char* argv[] )
{
    // initialize MORIS global communication manager
    gMorisComm = moris::Comm_Manager( &argc, &argv );

    // set severity level 0 - all outputs
    gLogger.initialize( 2 );

    // print banner
    moris::print_banner( argc, argv );

    // Every example expects to run in ITS working directory with a Catch2 test
    // spec (ctest supplies both). A bare invocation would run all examples'
    // cases in one directory - refuse instead.
    if ( argc < 2 )
    {
        if ( moris::par_rank() == 0 )
        {
            std::cerr << "EXA-test: pass a Catch2 test spec, e.g.  EXA-test.exe \"[EXA_Laplace_2D]\"\n"
                      << "          (run with --list-tests to see all; ctest supplies the spec and working dir)\n";
        }
        gMorisComm.finalize();
        return 1;
    }

    // Run Tests
    int tRet = Catch::Session().run( argc, argv );

    // emit the consolidated performance report (while MPI is still live)
    gPerfReporter.finalize();

    // finalize MORIS global communication manager
    gMorisComm.finalize();

    return tRet;
}
