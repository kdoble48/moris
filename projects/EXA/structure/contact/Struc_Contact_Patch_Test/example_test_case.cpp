/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 * Case matrix for the small-strain constant-stress contact patch test:
 * a 4 (method) x 4 (gamma) sweep, gCaseIndex = method*4 + gamma  (0..15).
 *
 *   method = gCaseIndex / 4 :  0 nitsche-symmetric, 1 nitsche-unsymmetric,
 *                              2 nitsche-neutral,    3 penalty
 *   gamma  = gCaseIndex % 4 :  0 -> 10, 1 -> 100, 2 -> 1000, 3 -> 10000
 *
 * The gamma sweep is thesis Argument B: post-gap-closure the only term coupling the two
 * bodies' stress is gamma*jump, so larger gamma should collapse the sigma_yy split toward
 * uniform. check_results() logs (method, gamma) + mean sigma_yy per body so the collapse
 * is readable straight from the logs.
 *
 * Each case is its own TEST_CASE (own Catch2 tag [EXA_Struc_Contact_Patch_Test_N]) and
 * calls fn_WRK_Workflow_Main_Interface EXACTLY ONCE. ArborX/Kokkos is not re-entrant
 * across a second full workflow in one process, so the EXA runner drives one workflow
 * per process (one tag per invocation); a single TEST_CASE looping over the cases
 * SIGABRTs at case 1. Mirrors ../Parabolic_Indenter/example_test_case.cpp.
 *
 * Patch-test criterion (reporting-only for now): sigma_yy is spatially UNIFORM across
 * the non-conforming contact interface. See the interior-sigma_yy metric in
 * check_results().
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
// -rdynamic dynamic symbol table. They must stay at global (non-namespaced) scope.
//
// These are DECLARED extern here (NOT defined): the EXA single-runner's EXA_Globals.cpp
// (feature/performance-report) DEFINES gInterpolationOrder / gCaseIndex /
// gPrintReferenceValues. Defining our own copies here would collide at link time in that
// runner. The runner provides the storage; we only read/write it. Mirrors Parabolic_Indenter.

// global variable for interpolation order
extern uint gInterpolationOrder;

// flag to print reference values
extern bool gPrintReferenceValues;

// test case index (selects contact method/bias in the deck)
extern uint gCaseIndex;

//---------------------------------------------------------------

int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_struc_contact_patch_test
{

//---------------------------------------------------------------
// PLACEHOLDER regression. The reference values below are per-case placeholders and MUST
// be harvested by the build/run agent: run each case once with gPrintReferenceValues = true,
// copy the printed node id / coordinates / displacement into the vectors below, then flip
// gHaveReference to true. Until then the strict nodal-displacement REQUIRE is skipped so the
// structural test (file opens, correct dimensionality, non-degenerate mesh) still guards the deck.

static const bool gHaveReference = false;

//---------------------------------------------------------------
// Interior sigma_yy sampling band (patch-test metric).
//
// The non-conforming interface sits at y ~= tInitialGapMiddle (0.503 in the deck). The cut
// element rows (VIS sets HMR_dummy_c_p1 / HMR_dummy_c_p2) straddle that line and give an
// UNRELIABLE stress readout (ELEMENTAL_AVG recovered a wrong-sign +3.0 on the interface cut
// row). We therefore restrict the sampled sigma_yy to the UNCUT interior of each body
// (VIS sets HMR_dummy_n_p1 for top, HMR_dummy_n_p2 for bottom) by excluding a band around
// the interface. sigma_yy is read from a NODAL-projected field (Exodus_IO_Helper exposes
// nodal fields only; there is no public elemental getter, and nodal projection is cleaner on
// cut cells than ELEMENTAL_AVG).
static const real gInterfaceY = 0.503;    // == tInitialGapMiddle in the deck
static const real gInteriorBandHalf = 0.06;    // exclude +/- band around the interface (cut rows)

void check_results()
{
    std::string tExoFileName =
            "Struc_Contact_Patch_Test_Case_" + std::to_string( gCaseIndex ) + ".e-s.0000";

    // decode the method x gamma matrix (gCaseIndex = method*4 + gamma) for logging
    uint               tMethodIndex   = gCaseIndex / 4;
    uint               tGammaIndex    = gCaseIndex % 4;
    const char *const  tMethodNames[] = { "nitsche-symmetric", "nitsche-unsymmetric", "nitsche-neutral", "penalty" };
    const real         tGammaValues[] = { 10.0, 100.0, 1000.0, 10000.0 };
    const char *const  tMethodName    = tMethodNames[ tMethodIndex ];
    const real         tGammaValue    = tGammaValues[ tGammaIndex ];

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", gCaseIndex, par_size() );
    MORIS_LOG_INFO( "  contact method %s, stabilization gamma %g", tMethodName, tGammaValue );
    MORIS_LOG_INFO( " " );

    // open and query exodus output file
    moris::mtk::Exodus_IO_Helper tExoIO( tExoFileName.c_str(), 0, false, false );

    // per-case reference node ids (PLACEHOLDER - harvest with gPrintReferenceValues).
    // 16 entries: one per method x gamma case (gCaseIndex 0..15).
    Vector< uint > tReferenceNodeId = { 294, 294, 294, 294, 294, 294, 294, 294,
        294, 294, 294, 294, 294, 294, 294, 294 };

    if ( gPrintReferenceValues )
    {
        std::cout << "Test case index: " << gCaseIndex << '\n';

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        std::cout << "Number of dimensions: " << tNumDims << '\n';
        std::cout << "Number of nodes     : " << tNumNodes << '\n';
        std::cout << "Number of elements  : " << tNumElems << '\n';

        moris::print( tExoIO.get_nodal_coordinate( tReferenceNodeId( gCaseIndex ) ), "Coordinates of reference point" );

        std::cout << "Time value: " << std::scientific << std::setprecision( 15 ) << tExoIO.get_time_value() << '\n';

        // nodal fields 2,3 are UX,UY (0,1 are the exodus coordinates)
        std::cout << "Displacement at reference point: "
                  << std::scientific << std::setprecision( 15 )
                  << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 2, 0 ) << ","
                  << tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 3, 0 ) << '\n';

        return;
    }

    // structural checks (always on): 2D, non-degenerate mesh
    uint tNumDims  = tExoIO.get_number_of_dimensions();
    uint tNumNodes = tExoIO.get_number_of_nodes();
    uint tNumElems = tExoIO.get_number_of_elements();

    REQUIRE( tNumDims == 2 );
    REQUIRE( tNumNodes > 0 );
    REQUIRE( tNumElems > 0 );

    //-----------------------------------------------------------
    // patch-test metric (reporting-only): interior sigma_yy uniformity.
    //
    // sigma_yy must be spatially UNIFORM (nu = 0 uniaxial compression). We read the
    // NODAL-projected sigma_yy fields "STRESSNodalTopY" / "STRESSNodalBottomY" (IQIs
    // IQIStressTopY / IQIStressBottomY, VIS Field_Type NODAL) over the UNCUT interior of
    // each body and report:
    //   - per-body spread  max(sigma_yy) - min(sigma_yy)   (should -> 0 for a passing patch test)
    //   - cross-body diff  |mean_top - mean_bottom|         (top and bottom carry the same stress)
    // Tolerance is intentionally loose/reporting-only while the loading setup is still being
    // validated (a separate triage owns the loading/BC decision).
    {
        real tTopMin = std::numeric_limits< real >::max();
        real tTopMax = -std::numeric_limits< real >::max();
        real tTopSum = 0.0;
        uint tTopCnt = 0;

        real tBotMin = std::numeric_limits< real >::max();
        real tBotMax = -std::numeric_limits< real >::max();
        real tBotSum = 0.0;
        uint tBotCnt = 0;

        // resolve the nodal sigma_yy field indices by name (added in VIS)
        uint tTopFieldIndex = tExoIO.get_field_index_by_name( "STRESSNodalTopY" );
        uint tBotFieldIndex = tExoIO.get_field_index_by_name( "STRESSNodalBottomY" );

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
                // top-body interior node (HMR_dummy_n_p1), away from the cut rows
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
                // bottom-body interior node (HMR_dummy_n_p2), away from the cut rows
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

        real tTopSpread = tTopMax - tTopMin;
        real tBotSpread = tBotMax - tBotMin;
        real tTopMean   = tTopSum / static_cast< real >( tTopCnt );
        real tBotMean   = tBotSum / static_cast< real >( tBotCnt );
        real tCrossBody = std::abs( tTopMean - tBotMean );

        MORIS_LOG_INFO( "Patch test sigma_yy (interior, case %d):", gCaseIndex );
        MORIS_LOG_INFO( "  top    interior nodes %d, mean %12.5e, spread(max-min) %12.5e",
                tTopCnt, tTopMean, tTopSpread );
        MORIS_LOG_INFO( "  bottom interior nodes %d, mean %12.5e, spread(max-min) %12.5e",
                tBotCnt, tBotMean, tBotSpread );
        MORIS_LOG_INFO( "  cross-body |mean_top - mean_bottom| %12.5e", tCrossBody );

        // compact, greppable one-liner for the gamma-collapse sweep (thesis Argument B):
        // method, gamma, mean sigma_yy per body, and the top-vs-bottom jump.
        MORIS_LOG_INFO( "GAMMA_SWEEP case %2d method %-20s gamma %8g mean_top %12.5e mean_bottom %12.5e jump %12.5e",
                gCaseIndex, tMethodName, tGammaValue, tTopMean, tBotMean, tCrossBody );

        // reporting-only: assert only that the metric is finite (setup still being validated).
        REQUIRE( std::isfinite( tTopSpread ) );
        REQUIRE( std::isfinite( tBotSpread ) );
        REQUIRE( std::isfinite( tCrossBody ) );
    }
    //-----------------------------------------------------------

    if ( !gHaveReference )
    {
        MORIS_LOG_INFO( "Struc_Contact_Patch_Test case %d: reference values not yet harvested - "
                        "skipping strict nodal regression (see check_results TODO).",
                gCaseIndex );
        return;
    }

    // per-case reference displacement (PLACEHOLDER - harvest with gPrintReferenceValues).
    // 16 entries: method x gamma matrix (gCaseIndex = method*4 + gamma).
    Vector< Matrix< DDRMat > > tReferenceDisplacement;
    for ( uint tCase = 0; tCase < 16; ++tCase )
    {
        tReferenceDisplacement.push_back( { { 0.0 }, { 0.0 } } );    // TODO harvest per case
    }

    Matrix< DDRMat > tActualDisplacement = {
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 2, 0 ) },
        { tExoIO.get_nodal_field_value( tReferenceNodeId( gCaseIndex ), 3, 0 ) },
    };

    real tRelDispDifference = norm( tActualDisplacement - tReferenceDisplacement( gCaseIndex ) ) / norm( tReferenceDisplacement( gCaseIndex ) );

    MORIS_LOG_INFO( "Check nodal displacements:  reference %12.5e, actual %12.5e, percent error %12.5e.",
            norm( tReferenceDisplacement( gCaseIndex ) ),
            norm( tActualDisplacement ),
            tRelDispDifference * 100.0 );

    REQUIRE( tRelDispDifference < 1.0e-5 );
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
    char tString2[] = "./Struc_Contact_Patch_Test.so";

    char *argv[ 2 ] = { tString1, tString2 };

    MORIS_LOG_INFO( " " );
    MORIS_LOG_INFO( "Executing Struc_Contact_Patch_Test case %d - %i Processors.", aCaseIndex, par_size() );
    MORIS_LOG_INFO( " " );

    // call to performance manager main interface (EXACTLY ONCE per process)
    fn_WRK_Workflow_Main_Interface( argc, argv );

    // check results
    check_results();
}

//---------------------------------------------------------------
// One TEST_CASE per (method, gamma) matrix entry, gCaseIndex = method*4 + gamma (0..15).
// Each sets gCaseIndex and runs the workflow EXACTLY ONCE. The EXA runner invokes the
// binary once per tag -> one workflow per process (ArborX/Kokkos re-entrancy constraint).
// The macro keeps each entry its own TEST_CASE / process while avoiding 16x boilerplate.

#ifdef MORIS_HAVE_ARBORX
#define CONTACT_PATCH_RUN( INDEX )                                                          \
    MORIS_ERROR( par_size() == 1, "Contact not implemented for parallel computation yet" ); \
    run_case( INDEX );
#else
#define CONTACT_PATCH_RUN( INDEX )                                                                            \
    MORIS_LOG_INFO( " " );                                                                                    \
    MORIS_LOG_INFO( "Struc_Contact_Patch_Test case " #INDEX ": Example skipped as Arborx not installed" );    \
    MORIS_LOG_INFO( " " );
#endif

#define CONTACT_PATCH_TEST_CASE( INDEX )                                                           \
    TEST_CASE( "Struc_Contact_Patch_Test_" #INDEX,                                                 \
            "[moris],[example],[structure],[contact],[EXA_Struc_Contact_Patch_Test_" #INDEX "]" )   \
    {                                                                                              \
        CONTACT_PATCH_RUN( INDEX )                                                                  \
    }

// method 0: nitsche-symmetric  x gamma {10,100,1000,10000}
CONTACT_PATCH_TEST_CASE( 0 )
CONTACT_PATCH_TEST_CASE( 1 )
CONTACT_PATCH_TEST_CASE( 2 )
CONTACT_PATCH_TEST_CASE( 3 )
// method 1: nitsche-unsymmetric x gamma {10,100,1000,10000}
CONTACT_PATCH_TEST_CASE( 4 )
CONTACT_PATCH_TEST_CASE( 5 )
CONTACT_PATCH_TEST_CASE( 6 )
CONTACT_PATCH_TEST_CASE( 7 )
// method 2: nitsche-neutral     x gamma {10,100,1000,10000}
CONTACT_PATCH_TEST_CASE( 8 )
CONTACT_PATCH_TEST_CASE( 9 )
CONTACT_PATCH_TEST_CASE( 10 )
CONTACT_PATCH_TEST_CASE( 11 )
// method 3: penalty             x gamma {10,100,1000,10000}
CONTACT_PATCH_TEST_CASE( 12 )
CONTACT_PATCH_TEST_CASE( 13 )
CONTACT_PATCH_TEST_CASE( 14 )
CONTACT_PATCH_TEST_CASE( 15 )

}    // namespace exa_struc_contact_patch_test
