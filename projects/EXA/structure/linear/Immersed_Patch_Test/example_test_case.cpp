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

#include <cmath>
#include <string>
#include <vector>

#include "exodusII.h"

#include "cl_Logger.hpp"                 // MRS/IOS/src
#include "cl_Communication_Tools.hpp"    // MRS/COM/src

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see doc/internal/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_immersed_patch_test
{
    // exact solution of the deck: u = A x + b with E = 1, nu = 0.3 in every phase
    const real tA[ 3 ][ 3 ] = { { 1.0e-3, 2.0e-4, 0.0 }, { 1.0e-4, -5.0e-4, 3.0e-4 }, { 0.0, 2.0e-4, 4.0e-4 } };

    real
    exact_sigma_xx()
    {
        real tNu     = 0.3;
        real tLambda = tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
        real tMu     = 1.0 / ( 2.0 * ( 1.0 + tNu ) );
        real tTrace  = tA[ 0 ][ 0 ] + tA[ 1 ][ 1 ] + tA[ 2 ][ 2 ];
        return tLambda * tTrace + 2.0 * tMu * tA[ 0 ][ 0 ];
    }

    // index (1-based) of a named variable of an exodus entity type, 0 if absent
    int
    variable_index( int aExoId, ex_entity_type aType, const std::string &aName )
    {
        int tNumVars = 0;
        ex_get_variable_param( aExoId, aType, &tNumVars );

        std::vector< std::vector< char > > tStorage( tNumVars, std::vector< char >( MAX_STR_LENGTH + 1, '\0' ) );
        std::vector< char * >              tNames( tNumVars );
        for ( int iVar = 0; iVar < tNumVars; iVar++ )
        {
            tNames[ iVar ] = tStorage[ iVar ].data();
        }
        if ( tNumVars > 0 )
        {
            ex_get_variable_names( aExoId, aType, tNumVars, tNames.data() );
        }

        for ( int iVar = 0; iVar < tNumVars; iVar++ )
        {
            if ( aName == tNames[ iVar ] )
            {
                return iVar + 1;
            }
        }
        return 0;
    }

    void
    check_results()
    {
        int   tWordSize = sizeof( real );
        int   tIOSize   = 0;
        float tVersion  = 0.0;
        int   tExoId    = ex_open( "Immersed_Patch_Test.exo", EX_READ, &tWordSize, &tIOSize, &tVersion );
        REQUIRE( tExoId >= 0 );

        int tLastStep = static_cast< int >( ex_inquire_int( tExoId, EX_INQ_TIME ) );
        REQUIRE( tLastStep >= 1 );

        // the displacement is exact up to the cells dropped from the residual (relative volume < 1e-6): the L2 error is
        // ~3e-10 for displacements of ~1e-3
        int tL2Index = variable_index( tExoId, EX_GLOBAL, "IQIDispErrorL2" );
        REQUIRE( tL2Index > 0 );
        int tNumGlobals = 0;
        ex_get_variable_param( tExoId, EX_GLOBAL, &tNumGlobals );
        std::vector< real > tGlobals( tNumGlobals, 0.0 );
        ex_get_var( tExoId, tLastStep, EX_GLOBAL, 1, 1, tNumGlobals, tGlobals.data() );
        CHECK( std::sqrt( tGlobals[ tL2Index - 1 ] ) < 1.0e-8 );

        // every output cell carries the exact stress, including the cells dropped from the residual
        int tStressIndex = variable_index( tExoId, EX_ELEM_BLOCK, "IQIStressXX" );
        REQUIRE( tStressIndex > 0 );

        int              tNumBlocks = static_cast< int >( ex_inquire_int( tExoId, EX_INQ_ELEM_BLK ) );
        std::vector< int > tBlockIds( tNumBlocks, 0 );
        ex_get_ids( tExoId, EX_ELEM_BLOCK, tBlockIds.data() );

        real tExact      = exact_sigma_xx();
        uint tNumCells   = 0;
        uint tNumWrong   = 0;
        real tMaxRelDiff = 0.0;

        for ( int iBlock = 0; iBlock < tNumBlocks; iBlock++ )
        {
            char tElemType[ MAX_STR_LENGTH + 1 ];
            int  tNumElems = 0, tNodesPerElem = 0, tEdges = 0, tFaces = 0, tAttributes = 0;
            ex_get_block( tExoId, EX_ELEM_BLOCK, tBlockIds[ iBlock ], tElemType, &tNumElems, &tNodesPerElem, &tEdges, &tFaces, &tAttributes );
            if ( tNumElems == 0 )
            {
                continue;
            }

            std::vector< real > tValues( tNumElems, 0.0 );
            ex_get_var( tExoId, tLastStep, EX_ELEM_BLOCK, tStressIndex, tBlockIds[ iBlock ], tNumElems, tValues.data() );

            for ( real iValue : tValues )
            {
                real tRelDiff = std::abs( iValue - tExact ) / std::abs( tExact );
                tMaxRelDiff   = std::max( tMaxRelDiff, std::isnan( tRelDiff ) ? 1.0 : tRelDiff );
                tNumWrong += ( std::isnan( tRelDiff ) or tRelDiff > 1.0e-4 ) ? 1 : 0;
                tNumCells++;
            }
        }
        ex_close( tExoId );

        MORIS_LOG_INFO( "Immersed_Patch_Test: %u output cells, %u off the exact stress, max relative deviation %e",
                tNumCells,
                tNumWrong,
                tMaxRelDiff );

        CHECK( tNumCells > 0 );
        CHECK( tNumWrong == 0 );
    }
}    // namespace exa_immersed_patch_test

//---------------------------------------------------------------

TEST_CASE( "Immersed_Patch_Test",
        "[moris],[example],[structure],[linear],[EXA_Immersed_Patch_Test]" )
{
    if ( par_size() == 1 )
    {
        int argc = 2;

        char tString1[] = "";
        char tString2[] = "./Immersed_Patch_Test.so";

        char *argv[ 2 ] = { tString1, tString2 };

        fn_WRK_Workflow_Main_Interface( argc, argv );

        exa_immersed_patch_test::check_results();
    }
}
