/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_IOS_Performance_Reporter.cpp
 *
 */

#include <catch.hpp>

#include <cstdio>
#include <fstream>
#include <string>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include "cl_Communication_Tools.hpp"       // COM/src
#include "cl_Performance_Reporter.hpp"      // MRS/IOS/src
#include "cl_Tracer.hpp"                    // MRS/IOS/src

// ----------------------------------------------------------------------------

TEST_CASE(
        "Performance_Reporter accumulates and writes a JSON report",
        "[moris],[ios],[performance_reporter]" )
{
    const std::string tJsonPath = "ut_perf_report.json";
    std::remove( tJsonPath.c_str() );

    // enable the reporter at full granularity
    gPerfReporter.initialize( 3, tJsonPath );

    // --- direct API: deterministic accumulation across repeated visits ---
    // ModuleA is visited twice: wall 2.0 + 3.0 = 5.0 s, mem 10.0 + (-4.0) = 6.0 MB
    gPerfReporter.record_region( "UT_PERF", "ModuleA", "Compute", 1, 2.0, 10.0 );
    gPerfReporter.record_region( "UT_PERF", "ModuleA", "Compute", 1, 3.0, -4.0 );
    // ModuleB is visited once
    gPerfReporter.record_region( "UT_PERF", "ModuleB", "Assemble", 1, 1.5, 5.0 );

    gPerfReporter.update_peak_rss( 123.0 );

    // --- integration: a real Tracer scope drives gLogger.sign_in/sign_out,
    // which in turn feeds gPerfReporter.record_region (proves the hook is wired) ---
    {
        moris::Tracer tTracer( "UT_PERF", "Traced", "Run" );
    }

    // aggregate across ranks, print the table, and write the JSON report
    gPerfReporter.finalize();

    // only rank 0 writes the report
    if ( moris::par_rank() == 0 )
    {
        std::ifstream tIfs( tJsonPath );
        REQUIRE( tIfs.good() );

        boost::property_tree::ptree tTree;
        boost::property_tree::read_json( tIfs, tTree );

        // run-level fields
        REQUIRE( tTree.get< int >( "level" ) == 3 );
        REQUIRE( tTree.get< double >( "peak_mem_mb" ) >= 123.0 );
        REQUIRE( tTree.count( "modules" ) == 1 );

        // find the ModuleA / ModuleB records and check the accumulated values
        bool tFoundA = false;
        bool tFoundB = false;
        bool tFoundTraced = false;

        for ( const auto& tModule : tTree.get_child( "modules" ) )
        {
            const std::string tEntity = tModule.second.get< std::string >( "entity" );
            const std::string tType   = tModule.second.get< std::string >( "type" );

            if ( tEntity == "UT_PERF" && tType == "ModuleA" )
            {
                tFoundA = true;
                CHECK( tModule.second.get< int >( "visits" ) == 2 );
                CHECK( tModule.second.get< double >( "wall_s" ) == Approx( 5.0 ) );
                CHECK( tModule.second.get< double >( "mem_delta_mb" ) == Approx( 6.0 ) );
                CHECK( tModule.second.get< double >( "mem_delta_peak_mb" ) == Approx( 10.0 ) );
            }
            else if ( tEntity == "UT_PERF" && tType == "ModuleB" )
            {
                tFoundB = true;
                CHECK( tModule.second.get< int >( "visits" ) == 1 );
                CHECK( tModule.second.get< double >( "wall_s" ) == Approx( 1.5 ) );
            }
            else if ( tEntity == "UT_PERF" && tType == "Traced" )
            {
                // the Tracer-driven region: proves Logger::sign_out -> record_region works
                tFoundTraced = true;
                CHECK( tModule.second.get< int >( "visits" ) == 1 );
            }
        }

        REQUIRE( tFoundA );
        REQUIRE( tFoundB );
        REQUIRE( tFoundTraced );

        // clean up the artifact
        std::remove( tJsonPath.c_str() );
    }
}
