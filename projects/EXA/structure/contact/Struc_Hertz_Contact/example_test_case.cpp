/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 * Case matrix for the 2D Hertzian-contact benchmark:
 * a 3 (penetration) x 3 (mesh) sweep, gCaseIndex = load*3 + mesh  (0..8).
 *
 *   load = gCaseIndex / 3 :  0 -> 0.0025, 1 -> 0.005, 2 -> 0.01   (penetration)
 *   mesh = gCaseIndex % 3 :  0 -> 40,     1 -> 80,    2 -> 160     (uniform NxN)
 *
 * Each case is its own TEST_CASE (own Catch2 tag [EXA_Struc_Hertz_Contact_N]) and calls
 * fn_WRK_Workflow_Main_Interface EXACTLY ONCE. ArborX/Kokkos is not re-entrant across a
 * second full workflow in one process, so the EXA runner drives one workflow per process
 * (one tag per invocation); a single TEST_CASE looping over the cases SIGABRTs at case 1.
 * Mirrors ../Struc_Contact_Patch_Test/example_test_case.cpp.
 *
 * check_results() is REPORTING-ONLY (finiteness): the deck writes a FULL_DISCONTINUOUS
 * (unzipped) VIS mesh, on which node ids / a hard-coded nodal-displacement reference table
 * are NOT a valid regression target. It logs the interior sigma_yy range per body and the
 * Hertz benchmark inputs (penetration, R = 1/(2*tFactor)) so a/R can be measured in post.
 */

#include <catch.hpp>

#include <algorithm>
#include <cmath>
#include <limits>

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"     // MRS/COM/src

#include "cl_Matrix.hpp"
#include "fn_norm.hpp"

using namespace moris;

//---------------------------------------------------------------
// Deck-visible globals: the input .so resolves these BY NAME against the executable's
// -rdynamic dynamic symbol table. DECLARED extern here (NOT defined): the EXA single-runner's
// EXA_Globals.cpp DEFINES gInterpolationOrder / gCaseIndex / gPrintReferenceValues. The runner
// provides the storage; we only read/write it. Mirrors Struc_Contact_Patch_Test.

// global variable for interpolation order
extern uint gInterpolationOrder;

// flag to print reference values
extern bool gPrintReferenceValues;

// test case index (gCaseIndex = load*3 + mesh, selects penetration/mesh in the deck)
extern uint gCaseIndex;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_struc_hertz_contact
{

//---------------------------------------------------------------
// Strict nodal-displacement regression is GATED OFF. The FULL_DISCONTINUOUS (unzipped) VIS mesh
// makes node ids / a hard-coded displacement table an invalid regression target (see the
// Struc_Contact_Patch_Test note). Keep this false; the structural + finiteness checks guard the deck.
static const bool gHaveReference = false;

//---------------------------------------------------------------
// Interior sigma_yy sampling band (reporting metric).
//
// The non-conforming interface sits at y ~= tInitialGapMiddle (0.503 in the deck). The cut
// element rows straddle that line and give an unreliable stress readout, so we restrict the
// sampled sigma_yy to the UNCUT interior of each body by excluding a band around the interface.
// sigma_yy is read from the NODAL-projected field STRESSNodal{Top,Bottom}YY (Exodus_IO_Helper
// exposes nodal fields only).
static const real gInterfaceY      = 0.503;    // == tInitialGapMiddle in the deck
static const real gInteriorBandHalf = 0.06;    // exclude +/- band around the interface (cut rows)

// Hertz benchmark inputs (mirror the deck constants) for logging a/R inputs.
static const real gFactor              = 1.0;    // == tFactor  -> R = 1/(2*tFactor) = 0.5
static const real gPenetrationLevels[] = { 0.0025, 0.005, 0.01 };
static const uint gMeshLevels[]        = { 40, 80, 160 };

void check_results()
{
    std::string tExoFileName =
            "Struc_Hertz_Contact_Case_" + std::to_string( gCaseIndex ) + ".e-s.0000";

    // decode the load x mesh matrix (gCaseIndex = load*3 + mesh) for logging
    uint tLoadIndex   = gCaseIndex / 3;
    uint tMeshIndex   = gCaseIndex % 3;
    real tPenetration = gPenetrationLevels[ tLoadIndex ];
    uint tNumElems    = gMeshLevels[ tMeshIndex ];
    real tRadius      = 1.0 / ( 2.0 * gFactor );

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", gCaseIndex, par_size() );
    MORIS_LOG_INFO( "  Hertz inputs: penetration %g, R = 1/(2*factor) %g, mesh %ux%u",
            tPenetration, tRadius, tNumElems, tNumElems );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file
    moris::mtk::Exodus_IO_Helper tExoIO( tExoFileName.c_str(), 0, false, false );

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << gCaseIndex << '\n';
        std::cout << "Number of dimensions: " << tExoIO.get_number_of_dimensions() << '\n';
        std::cout << "Number of nodes     : " << tExoIO.get_number_of_nodes() << '\n';
        std::cout << "Number of elements  : " << tExoIO.get_number_of_elements() << '\n';
        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << '\n';
        return;
    }

    // structural checks (always on): 2D, non-degenerate mesh
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElemsExo = tExoIO.get_number_of_elements();

    REQUIRE( tNumDims == 2 );
    REQUIRE( tNumNodes > 0 );
    REQUIRE( tNumElemsExo > 0 );

    //-----------------------------------------------------------
    // reporting-only metric: interior sigma_yy range per body (from the NODAL sigma_yy fields).
    // Not a pass/fail target (the loading + a/R are measured in post); we only assert finiteness.
    {
        real tTopMin = std::numeric_limits< real >::max();
        real tTopMax = -std::numeric_limits< real >::max();
        real tTopSum = 0.0;
        uint tTopCnt = 0;

        real tBotMin = std::numeric_limits< real >::max();
        real tBotMax = -std::numeric_limits< real >::max();
        real tBotSum = 0.0;
        uint tBotCnt = 0;

        uint tTopFieldIndex = tExoIO.get_field_index_by_name( "STRESSNodalTopYY" );
        uint tBotFieldIndex = tExoIO.get_field_index_by_name( "STRESSNodalBottomYY" );

        // On the FULL_DISCONTINUOUS (unzipped) mesh each per-body field is NaN on vertex
        // duplicates that belong to blocks where its IQI is not computed (the other body,
        // refinement-transition duplicates). Those NaNs are structural, not physics — skip
        // them but count them; a genuine solver blowup leaves ZERO finite samples and fails
        // the count REQUIREs below.
        uint tSkippedCnt = 0;

        // MORIS node ids are 0-based (Exodus_IO_Helper applies the on-disk QA id-offset itself)
        for ( uint tNodeId = 0; tNodeId < tNumNodes; ++tNodeId )
        {
            Matrix< DDRMat > tCoord = tExoIO.get_nodal_coordinate( tNodeId );
            real             tY     = tCoord( 1 );

            if ( tY > gInterfaceY + gInteriorBandHalf )
            {
                real tSig = tExoIO.get_nodal_field_value( tNodeId, tTopFieldIndex, 0 );
                if ( !std::isfinite( tSig ) )
                {
                    ++tSkippedCnt;
                    continue;
                }
                tTopMin = std::min( tTopMin, tSig );
                tTopMax = std::max( tTopMax, tSig );
                tTopSum += tSig;
                ++tTopCnt;
            }
            else if ( tY < gInterfaceY - gInteriorBandHalf )
            {
                real tSig = tExoIO.get_nodal_field_value( tNodeId, tBotFieldIndex, 0 );
                if ( !std::isfinite( tSig ) )
                {
                    ++tSkippedCnt;
                    continue;
                }
                tBotMin = std::min( tBotMin, tSig );
                tBotMax = std::max( tBotMax, tSig );
                tBotSum += tSig;
                ++tBotCnt;
            }
        }

        MORIS_LOG_INFO( "  sampled nodes skipped (field undefined on unzipped duplicate): %d", tSkippedCnt );

        REQUIRE( tTopCnt > 0 );
        REQUIRE( tBotCnt > 0 );

        real tTopMean = tTopSum / static_cast< real >( tTopCnt );
        real tBotMean = tBotSum / static_cast< real >( tBotCnt );

        MORIS_LOG_INFO( "Hertz sigma_yy (interior, case %d):", gCaseIndex );
        MORIS_LOG_INFO( "  top    interior nodes %d, mean %12.5e, min %12.5e, max %12.5e",
                tTopCnt, tTopMean, tTopMin, tTopMax );
        MORIS_LOG_INFO( "  bottom interior nodes %d, mean %12.5e, min %12.5e, max %12.5e",
                tBotCnt, tBotMean, tBotMin, tBotMax );

        // compact, greppable one-liner for the load x mesh sweep:
        MORIS_LOG_INFO( "HERTZ_SWEEP case %d penetration %g R %g mesh %u syy_top_min %12.5e syy_bot_min %12.5e",
                gCaseIndex, tPenetration, tRadius, tNumElems, tTopMin, tBotMin );

        // reporting-only: assert only finiteness (benchmark values measured in post).
        REQUIRE( std::isfinite( tTopMin ) );
        REQUIRE( std::isfinite( tTopMax ) );
        REQUIRE( std::isfinite( tBotMin ) );
        REQUIRE( std::isfinite( tBotMax ) );
    }
    //-----------------------------------------------------------

    if ( !gHaveReference )
    {
        MORIS_LOG_INFO( "Struc_Hertz_Contact case %d: reference values not harvested (FULL_DISCONTINUOUS "
                        "mesh -> nodal regression invalid) - skipping strict nodal check.",
                gCaseIndex );
        return;
    }
}

//---------------------------------------------------------------

static void run_case( uint aCaseIndex )
{
    gInterpolationOrder   = 1;
    gCaseIndex            = aCaseIndex;
    gPrintReferenceValues = false;

    // define command line call
    int argc = 2;

    char tString1[] = "";
    char tString2[] = "./Struc_Hertz_Contact.so";

    char *argv[ 2 ] = { tString1, tString2 };

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Struc_Hertz_Contact case %d - %i Processors.", aCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface (EXACTLY ONCE per process)
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    check_results();
}

//---------------------------------------------------------------
// One TEST_CASE per (load, mesh) matrix entry, gCaseIndex = load*3 + mesh (0..8).
// Each sets gCaseIndex and runs the workflow EXACTLY ONCE. The EXA runner invokes the
// binary once per tag -> one workflow per process (ArborX/Kokkos re-entrancy constraint).

#ifdef MORIS_HAVE_ARBORX
#define HERTZ_CONTACT_RUN( INDEX )                                                          \
    MORIS_ERROR( par_size() == 1, "Contact not implemented for parallel computation yet" ); \
    run_case( INDEX );
#else
#define HERTZ_CONTACT_RUN( INDEX )                                                                   \
    MORIS_LOG_INFO( " " );                                                                           \
    MORIS_LOG_INFO( "Struc_Hertz_Contact case " #INDEX ": Example skipped as Arborx not installed" ); \
    MORIS_LOG_INFO( " " );
#endif

#define HERTZ_CONTACT_TEST_CASE( INDEX )                                                        \
    TEST_CASE( "Struc_Hertz_Contact_" #INDEX,                                                    \
            "[moris],[example],[structure],[contact],[EXA_Struc_Hertz_Contact_" #INDEX "]" )     \
    {                                                                                           \
        HERTZ_CONTACT_RUN( INDEX )                                                              \
    }

// load 0: penetration 0.0025  x mesh {40, 80, 160}
HERTZ_CONTACT_TEST_CASE( 0 )
HERTZ_CONTACT_TEST_CASE( 1 )
HERTZ_CONTACT_TEST_CASE( 2 )
// load 1: penetration 0.005   x mesh {40, 80, 160}
HERTZ_CONTACT_TEST_CASE( 3 )
HERTZ_CONTACT_TEST_CASE( 4 )
HERTZ_CONTACT_TEST_CASE( 5 )
// load 2: penetration 0.01    x mesh {40, 80, 160}
HERTZ_CONTACT_TEST_CASE( 6 )
HERTZ_CONTACT_TEST_CASE( 7 )
HERTZ_CONTACT_TEST_CASE( 8 )

}    // namespace exa_struc_hertz_contact
