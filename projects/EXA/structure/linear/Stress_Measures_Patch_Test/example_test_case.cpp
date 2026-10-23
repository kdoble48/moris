/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 * Checks for Stress_Measures_Patch_Test.cpp (read its header first). Every expected value of case 0 is computed here
 * from the deck's data (A, E1, E2, nu, the plane normal), so each check reads: formula, number, and the defect that
 * would break it. One TEST_CASE per case index, each running the workflow once (one workflow per runner process).
 */

#include <catch.hpp>

#include <algorithm>
#include <cmath>
#include <map>
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

namespace exa_stress_measures_patch_test
{
    //------------------------------------------------------------------------------
    // the deck's data (Stress_Measures_Patch_Test.cpp)

    const real tA[ 3 ][ 3 ] = { { 1.0e-3, 2.0e-4, 0.0 }, { 1.0e-4, -5.0e-4, 3.0e-4 }, { 0.0, 2.0e-4, 4.0e-4 } };
    const real tE[ 2 ]      = { 1.0, 10.0 };    // material 1 (phase 0), material 2 (phases 1, 2, 3)
    const real tNu          = 0.3;
    const real tN0[ 3 ]     = { 1.0, 0.3, 0.2 };

    // reference strain energies of cases 1 and 2, per material. Not closed forms (the boundary data is
    // non-polynomial): harvested from case 2 (two CM objects) on the build that fixed MB-047. Case 1 (one CM object
    // named on both sides) must reproduce them; before the fix it gave 4.8231318e-06 / 5.8292285e-05 (+1.9e-4 in
    // material 2) after 6 Newton iterations instead of 2.
    const real tReferenceStrainEnergy[ 2 ] = { 4.823084980889959e-06, 5.828111206571638e-05 };

    using Tensor = std::vector< std::vector< real > >;

    //------------------------------------------------------------------------------
    // exact stresses of case 0

    // sigma = lambda tr( eps ) I + 2 mu eps, eps = sym( grad u )
    Tensor
    stress( real aE, const Tensor &aGrad )
    {
        real   tLambda = aE * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
        real   tMu     = aE / ( 2.0 * ( 1.0 + tNu ) );
        real   tTrace  = aGrad[ 0 ][ 0 ] + aGrad[ 1 ][ 1 ] + aGrad[ 2 ][ 2 ];
        Tensor tSigma( 3, std::vector< real >( 3, 0.0 ) );
        for ( uint i = 0; i < 3; i++ )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                tSigma[ i ][ j ] = 2.0 * tMu * 0.5 * ( aGrad[ i ][ j ] + aGrad[ j ][ i ] ) + ( i == j ? tLambda * tTrace : 0.0 );
            }
        }
        return tSigma;
    }

    // the two uniform stresses: material 1 from A, material 2 from A + a n^T (a: traction continuity, deck header)
    std::vector< Tensor >
    exact_stresses()
    {
        real tNorm = std::sqrt( tN0[ 0 ] * tN0[ 0 ] + tN0[ 1 ] * tN0[ 1 ] + tN0[ 2 ] * tN0[ 2 ] );
        real tN[ 3 ];
        for ( uint i = 0; i < 3; i++ )
        {
            tN[ i ] = tN0[ i ] / tNorm;
        }

        Tensor tGrad1( 3, std::vector< real >( 3 ) );
        for ( uint i = 0; i < 3; i++ )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                tGrad1[ i ][ j ] = tA[ i ][ j ];
            }
        }

        // r = ( C1 - C2 ) eps( A ) n,   a = ( r - ( lambda2 + mu2 ) / ( lambda2 + 2 mu2 ) ( n . r ) n ) / mu2
        Tensor tS1 = stress( tE[ 0 ], tGrad1 );
        Tensor tS2 = stress( tE[ 1 ], tGrad1 );
        real   tR[ 3 ];
        for ( uint i = 0; i < 3; i++ )
        {
            tR[ i ] = 0.0;
            for ( uint j = 0; j < 3; j++ )
            {
                tR[ i ] += ( tS1[ i ][ j ] - tS2[ i ][ j ] ) * tN[ j ];
            }
        }
        real tLambda2 = tE[ 1 ] * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) );
        real tMu2     = tE[ 1 ] / ( 2.0 * ( 1.0 + tNu ) );
        real tNR      = tR[ 0 ] * tN[ 0 ] + tR[ 1 ] * tN[ 1 ] + tR[ 2 ] * tN[ 2 ];

        Tensor tGrad2 = tGrad1;
        for ( uint i = 0; i < 3; i++ )
        {
            real tJump = ( tR[ i ] - ( tLambda2 + tMu2 ) / ( tLambda2 + 2.0 * tMu2 ) * tNR * tN[ i ] ) / tMu2;
            for ( uint j = 0; j < 3; j++ )
            {
                tGrad2[ i ][ j ] += tJump * tN[ j ];
            }
        }
        return { tS1, stress( tE[ 1 ], tGrad2 ) };
    }

    // principal stresses, largest first: roots of det( sigma - lambda I ) = 0, found by Jacobi rotations (independent
    // of the closed-form invariant formula used in MORIS)
    std::vector< real >
    principal_stresses( Tensor aS )
    {
        for ( uint iSweep = 0; iSweep < 50; iSweep++ )
        {
            for ( uint p = 0; p < 2; p++ )
            {
                for ( uint q = p + 1; q < 3; q++ )
                {
                    if ( std::abs( aS[ p ][ q ] ) < 1.0e-300 )
                    {
                        continue;
                    }
                    real tTheta = 0.5 * std::atan2( 2.0 * aS[ p ][ q ], aS[ q ][ q ] - aS[ p ][ p ] );
                    real c = std::cos( tTheta ), s = std::sin( tTheta );
                    for ( uint k = 0; k < 3; k++ )    // rows p, q
                    {
                        real tP = c * aS[ p ][ k ] - s * aS[ q ][ k ], tQ = s * aS[ p ][ k ] + c * aS[ q ][ k ];
                        aS[ p ][ k ] = tP, aS[ q ][ k ] = tQ;
                    }
                    for ( uint k = 0; k < 3; k++ )    // columns p, q
                    {
                        real tP = c * aS[ k ][ p ] - s * aS[ k ][ q ], tQ = s * aS[ k ][ p ] + c * aS[ k ][ q ];
                        aS[ k ][ p ] = tP, aS[ k ][ q ] = tQ;
                    }
                }
            }
        }
        std::vector< real > tL = { aS[ 0 ][ 0 ], aS[ 1 ][ 1 ], aS[ 2 ][ 2 ] };
        std::sort( tL.begin(), tL.end(), std::greater< real >() );
        return tL;
    }

    // von Mises: sqrt( 3/2 s:s ), s the deviator
    real
    von_mises( const Tensor &aS )
    {
        real tMean = ( aS[ 0 ][ 0 ] + aS[ 1 ][ 1 ] + aS[ 2 ][ 2 ] ) / 3.0;
        real tSum  = 0.0;
        for ( uint i = 0; i < 3; i++ )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                real tD = aS[ i ][ j ] - ( i == j ? tMean : 0.0 );
                tSum += tD * tD;
            }
        }
        return std::sqrt( 1.5 * tSum );
    }

    //------------------------------------------------------------------------------
    // exodus access

    int
    open_exodus( const std::string &aFile )
    {
        int   tWordSize = sizeof( real );
        int   tIOSize   = 0;
        float tVersion  = 0.0;
        int   tExoId    = ex_open( aFile.c_str(), EX_READ, &tWordSize, &tIOSize, &tVersion );
        REQUIRE( tExoId >= 0 );
        return tExoId;
    }

    // names of an entity's variables (EX_GLOBAL, EX_ELEM_BLOCK, EX_SIDE_SET) or of its blocks/sets
    std::vector< std::string >
    read_names( int aExoId, ex_entity_type aType, bool aEntityNames, int aCount )
    {
        std::vector< std::vector< char > > tStorage( aCount, std::vector< char >( MAX_STR_LENGTH + 1, '\0' ) );
        std::vector< char * >              tPointers( aCount );
        for ( int i = 0; i < aCount; i++ )
        {
            tPointers[ i ] = tStorage[ i ].data();
        }
        if ( aCount > 0 )
        {
            if ( aEntityNames )
            {
                ex_get_names( aExoId, aType, tPointers.data() );
            }
            else
            {
                ex_get_variable_names( aExoId, aType, aCount, tPointers.data() );
            }
        }
        return std::vector< std::string >( tPointers.begin(), tPointers.end() );
    }

    std::vector< std::string >
    variable_names( int aExoId, ex_entity_type aType )
    {
        int tNum = 0;
        ex_get_variable_param( aExoId, aType, &tNum );
        return read_names( aExoId, aType, false, tNum );
    }

    std::map< std::string, real >
    global_values( int aExoId, int aStep )
    {
        std::vector< std::string > tNames = variable_names( aExoId, EX_GLOBAL );
        std::vector< real >        tValues( tNames.size(), 0.0 );
        if ( !tNames.empty() )
        {
            ex_get_var( aExoId, aStep, EX_GLOBAL, 1, 1, tNames.size(), tValues.data() );
        }
        std::map< std::string, real > tMap;
        for ( uint i = 0; i < tNames.size(); i++ )
        {
            tMap[ tNames[ i ] ] = tValues[ i ];
        }
        return tMap;
    }

    // values of an element variable over the blocks of one material (block names end in _p<phase>)
    std::vector< real >
    element_values( int aExoId, int aStep, const std::string &aVariable, uint aMaterial )
    {
        std::vector< std::string > tVars  = variable_names( aExoId, EX_ELEM_BLOCK );
        int                        tIndex = std::find( tVars.begin(), tVars.end(), aVariable ) - tVars.begin() + 1;
        REQUIRE( tIndex <= static_cast< int >( tVars.size() ) );

        int                        tNumBlocks = static_cast< int >( ex_inquire_int( aExoId, EX_INQ_ELEM_BLK ) );
        std::vector< int >         tIds( tNumBlocks );
        std::vector< std::string > tBlockNames = read_names( aExoId, EX_ELEM_BLOCK, true, tNumBlocks );
        ex_get_ids( aExoId, EX_ELEM_BLOCK, tIds.data() );

        std::vector< real > tAll;
        for ( int iB = 0; iB < tNumBlocks; iB++ )
        {
            int  tPhase    = std::stoi( tBlockNames[ iB ].substr( tBlockNames[ iB ].rfind( "_p" ) + 2 ) );
            uint tMaterial = tPhase == 0 ? 1 : 2;
            char tType[ MAX_STR_LENGTH + 1 ];
            int  tNumElems = 0, tNodes = 0, tEdges = 0, tFaces = 0, tAttributes = 0;
            ex_get_block( aExoId, EX_ELEM_BLOCK, tIds[ iB ], tType, &tNumElems, &tNodes, &tEdges, &tFaces, &tAttributes );
            if ( tMaterial != aMaterial or tNumElems == 0 )
            {
                continue;
            }
            std::vector< real > tValues( tNumElems, 0.0 );
            ex_get_var( aExoId, aStep, EX_ELEM_BLOCK, tIndex, tIds[ iB ], tNumElems, tValues.data() );
            tAll.insert( tAll.end(), tValues.begin(), tValues.end() );
        }
        return tAll;
    }

    // largest magnitude of a side-set variable over all side sets that carry it
    real
    max_abs_side_set_value( int aExoId, int aStep, const std::string &aVariable )
    {
        std::vector< std::string > tVars  = variable_names( aExoId, EX_SIDE_SET );
        int                        tIndex = std::find( tVars.begin(), tVars.end(), aVariable ) - tVars.begin() + 1;
        REQUIRE( tIndex <= static_cast< int >( tVars.size() ) );

        int                tNumSets = static_cast< int >( ex_inquire_int( aExoId, EX_INQ_SIDE_SETS ) );
        std::vector< int > tIds( tNumSets );
        ex_get_ids( aExoId, EX_SIDE_SET, tIds.data() );

        // the side-set variable truth table says which sets carry the variable
        std::vector< int > tTable( tNumSets * tVars.size(), 0 );
        ex_get_truth_table( aExoId, EX_SIDE_SET, tNumSets, tVars.size(), tTable.data() );

        real tMax = 0.0;
        for ( int iS = 0; iS < tNumSets; iS++ )
        {
            int tNumSides = 0, tNumDist = 0;
            ex_get_set_param( aExoId, EX_SIDE_SET, tIds[ iS ], &tNumSides, &tNumDist );
            if ( tNumSides == 0 or tTable[ iS * tVars.size() + tIndex - 1 ] == 0 )
            {
                continue;
            }
            std::vector< real > tValues( tNumSides, 0.0 );
            ex_get_var( aExoId, aStep, EX_SIDE_SET, tIndex, tIds[ iS ], tNumSides, tValues.data() );
            for ( real v : tValues )
            {
                tMax = std::max( tMax, std::isfinite( v ) ? std::abs( v ) : 1.0e300 );
            }
        }
        return tMax;
    }

    //------------------------------------------------------------------------------

    void
    check_case_0()
    {
        std::vector< Tensor > tSigma = exact_stresses();
        real                  tScale = 0.0;
        for ( const Tensor &tS : tSigma )
        {
            for ( const auto &tRow : tS )
            {
                for ( real v : tRow )
                {
                    tScale = std::max( tScale, std::abs( v ) );
                }
            }
        }
        const real tTolerance = 1.0e-10 * tScale;

        int tExoId = open_exodus( "Stress_Measures_Patch_Test_0.exo" );
        int tStep  = static_cast< int >( ex_inquire_int( tExoId, EX_INQ_TIME ) );

        // the displacement is the exact field (the L2 error is integrated squared)
        std::map< std::string, real > tGlobals = global_values( tExoId, tStep );
        CHECK( std::sqrt( tGlobals[ "IQIDispErrorL2" ] ) < 1.0e-12 );

        for ( uint iM = 0; iM < 2; iM++ )
        {
            std::string         tM          = "_M" + std::to_string( iM + 1 );
            std::vector< real > tPrincipal  = principal_stresses( tSigma[ iM ] );
            std::map< std::string, real > tExpected = {
                { "SXX", tSigma[ iM ][ 0 ][ 0 ] },
                { "SVM", von_mises( tSigma[ iM ] ) },
                { "SP1", tPrincipal[ 0 ] },    // MB-046: off by up to 1.9 % when det( sigma ) used s_yz for s_zz
                { "SP2", tPrincipal[ 1 ] },
                { "SP3", tPrincipal[ 2 ] },
            };

            for ( const auto &[ tName, tExact ] : tExpected )
            {
                std::vector< real > tValues = element_values( tExoId, tStep, tName + tM, iM + 1 );
                real                tMaxDev = 0.0;
                for ( real v : tValues )
                {
                    tMaxDev = std::max( tMaxDev, std::isfinite( v ) ? std::abs( v - tExact ) : 1.0e300 );
                }
                MORIS_LOG_INFO( "Stress_Measures_Patch_Test: %s on %zu cells, exact %+.10e, max deviation %.2e",
                        ( tName + tM ).c_str(),
                        tValues.size(),
                        tExact,
                        tMaxDev );
                CHECK( tValues.size() > 0 );
                CHECK( tMaxDev < tTolerance );
            }

        }
        CHECK( std::abs( tGlobals[ "VOL_M1" ] + tGlobals[ "VOL_M2" ] - 8.0 ) < 1.0e-10 );
        ex_close( tExoId );

        // traction is continuous: the jump vanishes on the bimaterial interface and on the material-2 interface,
        // whose follower names the leader's CM (MB-048: the sum t_L( n ) + t_F( n ) was returned, 2 |sigma n|;
        // MB-047: with one CM object both terms came from the same cached traction)
        int tFacetId = open_exodus( "Stress_Measures_Patch_Test_0_faceted.exo" );
        int tFacetStep = static_cast< int >( ex_inquire_int( tFacetId, EX_INQ_TIME ) );
        for ( std::string tJump : { "JT_BI", "JT_M2" } )
        {
            real tMax = max_abs_side_set_value( tFacetId, tFacetStep, tJump );
            MORIS_LOG_INFO( "Stress_Measures_Patch_Test: max |%s| = %.2e", tJump.c_str(), tMax );
            CHECK( tMax < tTolerance );
        }
        ex_close( tFacetId );
    }

    void
    check_case_3()
    {
        std::vector< Tensor > tSigma = exact_stresses();

        int tExoId = open_exodus( "Stress_Measures_Patch_Test_3.exo" );
        int tStep  = static_cast< int >( ex_inquire_int( tExoId, EX_INQ_TIME ) );
        std::map< std::string, real > tGlobals = global_values( tExoId, tStep );
        ex_close( tExoId );

        // MAX_NORMAL_STRESS with ( sigma_xx / 1 - 0 )^1 integrates sigma_xx: the mean is the exact sigma_xx
        // (MB-049: the bounds check aborted every component index, so this case did not run at all)
        for ( uint iM = 0; iM < 2; iM++ )
        {
            std::string tM    = "_M" + std::to_string( iM + 1 );
            real        tMean = tGlobals[ "MAXSXX" + tM ] / tGlobals[ "VOL" + tM ];
            MORIS_LOG_INFO( "Stress_Measures_Patch_Test: MAXSXX%s / VOL%s = %+.10e, exact %+.10e", tM.c_str(), tM.c_str(), tMean, tSigma[ iM ][ 0 ][ 0 ] );
            CHECK( std::abs( tMean - tSigma[ iM ][ 0 ][ 0 ] ) < 1.0e-10 * std::abs( tSigma[ iM ][ 0 ][ 0 ] ) );
        }
    }

    void
    check_case_1_or_2( uint aCase )
    {
        int tExoId = open_exodus( "Stress_Measures_Patch_Test_" + std::to_string( aCase ) + ".exo" );
        int tStep  = static_cast< int >( ex_inquire_int( tExoId, EX_INQ_TIME ) );
        std::map< std::string, real > tGlobals = global_values( tExoId, tStep );
        ex_close( tExoId );

        for ( uint iM = 0; iM < 2; iM++ )
        {
            std::string tName = "SE_M" + std::to_string( iM + 1 );
            real        tRel  = std::abs( tGlobals[ tName ] - tReferenceStrainEnergy[ iM ] ) / tReferenceStrainEnergy[ iM ];
            MORIS_LOG_INFO( "Stress_Measures_Patch_Test case %u: %s = %.15e, reference %.15e, relative difference %.2e",
                    aCase,
                    tName.c_str(),
                    tGlobals[ tName ],
                    tReferenceStrainEnergy[ iM ],
                    tRel );
            CHECK( tRel < 1.0e-9 );
        }
    }

    void
    run_case( uint aCase )
    {
        gCaseIndex = aCase;

        int  argc       = 2;
        char tString1[] = "";
        char tString2[] = "./Stress_Measures_Patch_Test.so";
        char *argv[ 2 ] = { tString1, tString2 };

        MORIS_LOG_INFO( "Executing Stress_Measures_Patch_Test case %u - %i Processors.", aCase, par_size() );

        fn_WRK_Workflow_Main_Interface( argc, argv );

        if ( aCase == 0 )
        {
            check_case_0();
        }
        else if ( aCase == 3 )
        {
            check_case_3();
        }
        else
        {
            check_case_1_or_2( aCase );
        }
    }
}    // namespace exa_stress_measures_patch_test

//---------------------------------------------------------------
// one TEST_CASE per case index; the EXA runner runs one per process (CASES 0 1 2 3 in CMakeLists.txt)

TEST_CASE( "Stress_Measures_Patch_Test_0", "[moris],[example],[structure],[linear],[EXA_Stress_Measures_Patch_Test_0]" )
{
    if ( par_size() == 1 )
    {
        exa_stress_measures_patch_test::run_case( 0 );
    }
}

TEST_CASE( "Stress_Measures_Patch_Test_1", "[moris],[example],[structure],[linear],[EXA_Stress_Measures_Patch_Test_1]" )
{
    if ( par_size() == 1 )
    {
        exa_stress_measures_patch_test::run_case( 1 );
    }
}

TEST_CASE( "Stress_Measures_Patch_Test_2", "[moris],[example],[structure],[linear],[EXA_Stress_Measures_Patch_Test_2]" )
{
    if ( par_size() == 1 )
    {
        exa_stress_measures_patch_test::run_case( 2 );
    }
}

TEST_CASE( "Stress_Measures_Patch_Test_3", "[moris],[example],[structure],[linear],[EXA_Stress_Measures_Patch_Test_3]" )
{
    if ( par_size() == 1 )
    {
        exa_stress_measures_patch_test::run_case( 3 );
    }
}
