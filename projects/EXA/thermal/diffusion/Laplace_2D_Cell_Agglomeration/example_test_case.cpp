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
#include <fstream>
#include <string>
#include <vector>

#include <exodusII.h>

#include "cl_Logger.hpp"                  // MRS/IOS/src
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include "cl_Communication_Tools.hpp"     // MRS/COM/src

#include "cl_Matrix.hpp"

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see doc/internal/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char *argv[] );

//---------------------------------------------------------------
// Everything below is TU-local to this example; names may repeat across examples.

namespace exa_laplace_2d_cell_agglomeration
{

    //---------------------------------------------------------------

    // number of values in a text file written by "SOL_save_operator_to_matlab";
    // for the residual file this is the number of unknowns
    uint
    count_values_in_file( const std::string &aFileName )
    {
        std::ifstream tFile( aFileName );
        REQUIRE( tFile.is_open() );

        uint tNumValues = 0;
        real tValue     = 0.0;
        while ( tFile >> tValue )
        {
            tNumValues++;
        }
        return tNumValues;
    }

    //---------------------------------------------------------------

    // index (0-based) of a name in a list of exodus names; -1 if the name is not in the list
    int
    find_exodus_name(
            int                aExoId,
            ex_entity_type     aEntityType,
            bool               aIsVariableName,
            int                aNumNames,
            const std::string &aName )
    {
        std::vector< std::vector< char > > tNameStorage( aNumNames, std::vector< char >( MAX_STR_LENGTH + 1, '\0' ) );
        std::vector< char * >              tNames( aNumNames );
        for ( int iName = 0; iName < aNumNames; iName++ )
        {
            tNames[ iName ] = tNameStorage[ iName ].data();
        }

        if ( aIsVariableName )
        {
            ex_get_variable_names( aExoId, aEntityType, aNumNames, tNames.data() );
        }
        else
        {
            ex_get_names( aExoId, aEntityType, tNames.data() );
        }

        for ( int iName = 0; iName < aNumNames; iName++ )
        {
            if ( aName == std::string( tNames[ iName ] ) )
            {
                return iName;
            }
        }
        return -1;
    }

    //---------------------------------------------------------------

    // Reads the aggregates that the cell agglomeration built from the interpolation mesh written by XTK
    // ("visualize_cell_association"): for every solid cell, its subphase-group index and the index of its root.
    // A cell whose root is another cell was agglomerated; a root of -1 means that no root was assigned.
    void
    count_agglomerated_cells(
            const std::string &aExoFileName,
            uint              &aNumSolidCells,
            uint              &aNumAgglomeratedCells,
            uint              &aNumCellsWithoutRoot )
    {
        aNumSolidCells        = 0;
        aNumAgglomeratedCells = 0;
        aNumCellsWithoutRoot  = 0;

        int   tCompWordSize = sizeof( real );
        int   tIOWordSize   = 0;
        float tVersion      = 0.0;
        int   tExoId        = ex_open( aExoFileName.c_str(), EX_READ, &tCompWordSize, &tIOWordSize, &tVersion );
        REQUIRE( tExoId >= 0 );

        ex_set_max_name_length( tExoId, MAX_STR_LENGTH );

        // solid block of the interpolation mesh
        int tNumBlocks = ex_inquire_int( tExoId, EX_INQ_ELEM_BLK );
        int tBlock     = find_exodus_name( tExoId, EX_ELEM_BLOCK, false, tNumBlocks, "HMR_dummy_p0" );
        REQUIRE( tBlock >= 0 );

        std::vector< int > tBlockIds( tNumBlocks );
        ex_get_ids( tExoId, EX_ELEM_BLOCK, tBlockIds.data() );

        ex_block tBlockInfo;
        tBlockInfo.id   = tBlockIds[ tBlock ];
        tBlockInfo.type = EX_ELEM_BLOCK;
        ex_get_block_param( tExoId, &tBlockInfo );
        int tNumCells = (int)tBlockInfo.num_entry;

        // elemental fields: subphase-group index and root index with respect to B-spline mesh 0
        int tNumVariables = 0;
        ex_get_variable_param( tExoId, EX_ELEM_BLOCK, &tNumVariables );
        int tSpgVariable  = find_exodus_name( tExoId, EX_ELEM_BLOCK, true, tNumVariables, "SPG_index0" );
        int tRootVariable = find_exodus_name( tExoId, EX_ELEM_BLOCK, true, tNumVariables, "Root_Cell_Index_B0" );
        REQUIRE( tSpgVariable >= 0 );
        REQUIRE( tRootVariable >= 0 );

        std::vector< real > tSpgIndices( tNumCells );
        std::vector< real > tRootIndices( tNumCells );
        ex_get_var( tExoId, 1, EX_ELEM_BLOCK, tSpgVariable + 1, tBlockInfo.id, tNumCells, tSpgIndices.data() );
        ex_get_var( tExoId, 1, EX_ELEM_BLOCK, tRootVariable + 1, tBlockInfo.id, tNumCells, tRootIndices.data() );

        ex_close( tExoId );

        aNumSolidCells = (uint)tNumCells;
        for ( int iCell = 0; iCell < tNumCells; iCell++ )
        {
            if ( tRootIndices[ iCell ] < -0.5 )
            {
                aNumCellsWithoutRoot++;
            }
            else if ( std::abs( tRootIndices[ iCell ] - tSpgIndices[ iCell ] ) > 0.5 )
            {
                aNumAgglomeratedCells++;
            }
        }
    }

    //---------------------------------------------------------------

    void
    check_results(
            const std::string &aOutputName,
            uint               aTestCaseIndex )
    {
        MORIS_LOG_INFO( " " );
        MORIS_LOG_INFO( "Checking Results - Test Case %d on %i processor.", aTestCaseIndex, par_size() );
        MORIS_LOG_INFO( " " );

        // open and query exodus output file (set verbose to true to get basic mesh information)
        moris::mtk::Exodus_IO_Helper tExoIO( ( aOutputName + ".exo" ).c_str(), 0, false, false );

        uint tNumDims  = tExoIO.get_number_of_dimensions();
        uint tNumNodes = tExoIO.get_number_of_nodes();
        uint tNumElems = tExoIO.get_number_of_elements();

        // global variables in the order of the VIS field list: L2ErrSq, H1ErrSq, Volume
        real tL2ErrorSquared = tExoIO.get_global_variable( 0, 0 );
        real tH1ErrorSquared = tExoIO.get_global_variable( 1, 0 );
        real tSolidVolume    = tExoIO.get_global_variable( 2, 0 );

        // number of unknowns
        uint tNumUnknowns = count_values_in_file( aOutputName + ".0.0.res.dat" );

        // aggregates
        uint tNumSolidCells        = 0;
        uint tNumAgglomeratedCells = 0;
        uint tNumCellsWithoutRoot  = 0;
        count_agglomerated_cells( "xtk_temp_ip.exo", tNumSolidCells, tNumAgglomeratedCells, tNumCellsWithoutRoot );

        if ( gPrintReferenceValues )
        {
            std::cout << "Test case index: " << aTestCaseIndex << '\n';
            std::cout << "Number of dimensions: " << tNumDims << '\n';
            std::cout << "Number of nodes     : " << tNumNodes << '\n';
            std::cout << "Number of elements  : " << tNumElems << '\n';
            std::cout << "Number of unknowns  : " << tNumUnknowns << '\n';
            std::cout << "Solid cells         : " << tNumSolidCells << '\n';
            std::cout << "Agglomerated cells  : " << tNumAgglomeratedCells << '\n';
            std::cout << "Cells without root  : " << tNumCellsWithoutRoot << '\n';
            std::cout << "L2 error squared    : " << std::scientific << std::setprecision( 15 ) << tL2ErrorSquared << '\n';
            std::cout << "H1 error squared    : " << std::scientific << std::setprecision( 15 ) << tH1ErrorSquared << '\n';
            std::cout << "Solid volume        : " << std::scientific << std::setprecision( 15 ) << tSolidVolume << '\n';

            return;
        }

        // reference values: test case 0 linear, test case 1 quadratic
        Vector< uint > tReferenceNumDims     = { 2, 2 };
        Vector< uint > tReferenceNumNodes    = { 354, 984 };
        Vector< uint > tReferenceNumElems    = { 280, 280 };
        Vector< uint > tReferenceNumUnknowns = { 231, 264 };

        // of the 220 solid cells, the 20 cells of the cut column are attached to their left neighbour
        uint tReferenceNumSolidCells        = 220;
        uint tReferenceNumAgglomeratedCells = 20;

        Vector< real > tReferenceL2ErrorSquared = { 9.922945117516934e-08, 4.681941283696269e-11 };
        Vector< real > tReferenceH1ErrorSquared = { 3.780860906456499e-03, 1.351853690506223e-06 };
        real           tReferenceSolidVolume    = 0.50005;

        CHECK( tNumDims == tReferenceNumDims( aTestCaseIndex ) );
        CHECK( tNumNodes == tReferenceNumNodes( aTestCaseIndex ) );
        CHECK( tNumElems == tReferenceNumElems( aTestCaseIndex ) );

        // without agglomeration there are 252 (linear) and 286 (quadratic) unknowns
        CHECK( tNumUnknowns == tReferenceNumUnknowns( aTestCaseIndex ) );

        CHECK( tNumSolidCells == tReferenceNumSolidCells );
        CHECK( tNumAgglomeratedCells == tReferenceNumAgglomeratedCells );
        CHECK( tNumCellsWithoutRoot == 0 );

        real tRelL2Difference     = ( tL2ErrorSquared - tReferenceL2ErrorSquared( aTestCaseIndex ) ) / tReferenceL2ErrorSquared( aTestCaseIndex );
        real tRelH1Difference     = ( tH1ErrorSquared - tReferenceH1ErrorSquared( aTestCaseIndex ) ) / tReferenceH1ErrorSquared( aTestCaseIndex );
        real tRelVolumeDifference = ( tSolidVolume - tReferenceSolidVolume ) / tReferenceSolidVolume;

        MORIS_LOG_INFO( "Check L2 error squared:  reference %12.5e, actual %12.5e, percent error %12.5e.",
                tReferenceL2ErrorSquared( aTestCaseIndex ),
                tL2ErrorSquared,
                tRelL2Difference * 100.0 );
        MORIS_LOG_INFO( "Check H1 error squared:  reference %12.5e, actual %12.5e, percent error %12.5e.",
                tReferenceH1ErrorSquared( aTestCaseIndex ),
                tH1ErrorSquared,
                tRelH1Difference * 100.0 );
        MORIS_LOG_INFO( "Check solid volume:      reference %12.5e, actual %12.5e, percent error %12.5e.",
                tReferenceSolidVolume,
                tSolidVolume,
                tRelVolumeDifference * 100.0 );

        CHECK( std::abs( tRelL2Difference ) < 1.0e-6 );
        CHECK( std::abs( tRelH1Difference ) < 1.0e-6 );
        CHECK( std::abs( tRelVolumeDifference ) < 1.0e-10 );
    }

    //---------------------------------------------------------------

    // runs the deck for one variant and checks its results
    void
    run_and_check(
            uint               aTestCaseIndex,
            const std::string &aOutputName )
    {
        // set all globals this test case or its deck consumes (rule R4)
        gTestIndex            = aTestCaseIndex;
        gPrintReferenceValues = false;

        // define command line call
        int  argc       = 2;
        char tString1[] = "";
        char tString2[] = "./Laplace_Cell_Agglomeration.so";

        char *argv[ 2 ] = { tString1, tString2 };

        // call to performance manager main interface
        int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );

        // check
        REQUIRE( tRet == 0 );

        check_results( aOutputName, aTestCaseIndex );
    }

    //---------------------------------------------------------------

    TEST_CASE( "Laplace_Cell_Agglomeration_Linear",
            "[moris],[example],[thermal],[Laplace_2D_Cell_Agglomeration],[Laplace_Cell_Agglomeration_Linear],[EXA_Laplace_2D_Cell_Agglomeration]" )
    {
        MORIS_LOG_INFO( " " );
        MORIS_LOG_INFO( "Executing Laplace_Cell_Agglomeration: linear elements on linear B-splines - %i Processors.", par_size() );
        MORIS_LOG_INFO( " " );

        run_and_check( 0, "Laplace_Cell_Agglomeration_Linear" );
    }

    //---------------------------------------------------------------

    // needs the cell agglomeration for more than one follower basis per node (MB-108)
    TEST_CASE( "Laplace_Cell_Agglomeration_Quadratic",
            "[moris],[example],[thermal],[Laplace_2D_Cell_Agglomeration],[Laplace_Cell_Agglomeration_Quadratic],[EXA_Laplace_2D_Cell_Agglomeration]" )
    {
        MORIS_LOG_INFO( " " );
        MORIS_LOG_INFO( "Executing Laplace_Cell_Agglomeration: quadratic elements on quadratic B-splines - %i Processors.", par_size() );
        MORIS_LOG_INFO( " " );

        run_and_check( 1, "Laplace_Cell_Agglomeration_Quadratic" );
    }

}    // namespace exa_laplace_2d_cell_agglomeration
