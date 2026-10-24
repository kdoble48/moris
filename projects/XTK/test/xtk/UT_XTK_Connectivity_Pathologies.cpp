/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Connectivity_Pathologies.cpp
 *
 * Robustness tests of the three connectivity layers built on top of the conformal cut:
 *
 *   1. IG cells -> subphases (flood fill per background cell, Integration_Mesh_Generator::flood_fill_ig_cell_group),
 *   2. the subphase neighbourhood across background faces (Integration_Mesh_Generator::construct_subphase_neighborhood),
 *   3. subphases -> enrichment levels per basis function (legacy Enrichment::perform_basis_cluster_enrichment, the path
 *      production decks use).
 *
 * Each layer is compared with a reference computed from the raw IG mesh geometry, never from XTK's facet tables:
 *
 *   - subphases: union-find over same-phase IG tets of one background cell that share a triangle (sorted vertex indices);
 *   - neighbourhood: two subphases in neighbouring background cells are adjacent iff they have the same bulk phase and
 *     their footprints on the shared face overlap with positive area (identical triangles on a conforming face, clipped
 *     polygons across a refinement transition, the whole face for an uncut cell);
 *   - enrichment: per basis function, union-find over the subphases in its support with the reference adjacency. XTK's
 *     levels (Enrichment_Data::mSubphaseBGBasisEnrLev) must equal these components as partitions. A level that holds
 *     two reference components couples disconnected material through one coefficient (cross-talk).
 *
 * Multi-rank cases gather canonical keys (background cell id, bulk phase, rounded IG cell centroids) so that the
 * subphase partition, the neighbourhood, the enrichment components and the output sets can be compared between rank
 * counts. Run a sweep on 1, 2 and 4 ranks with MORIS_XTK_CONNECTIVITY_CSV=<file> and compare the files.
 *
 * The flood-fill fuzz ([XTK_CONNECTIVITY_FLOOD_FILL]) is a fast regular test. Everything else is hidden ([.]): the
 * sweeps characterize current behaviour, and some cases demonstrate known defects.
 *
 */

#include "catch.hpp"

#include "cl_XTK_Connectivity_Test_Harness.hpp"

namespace moris::xtk
{
    using namespace connectivity_test;

    // ----------------------------------------------------------------------------------
    // section 2 on the single-hex sweeps

    TEST_CASE( "XTK connectivity sweep: single hex, corner values, three GEN regimes", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_HEX]" )
    {
        if ( par_size() != 1 )
        {
            return;
        }
        Grid_Mesh tHex( { 1, 1, 1 }, { 1.0, 1.0, 1.0 } );

        // MORIS_XTK_CONNECTIVITY_HEX_SETS=signs,near_zero (default both)
        const char* tSetsEnv = std::getenv( "MORIS_XTK_CONNECTIVITY_HEX_SETS" );
        std::string tSets    = tSetsEnv ? tSetsEnv : "signs,near_zero";

        Vector< std::pair< std::string, Vector< real > > > tValueSets;
        if ( tSets.find( "signs" ) != std::string::npos )
        {
            tValueSets.push_back( { "signs {-1,0,1}", { -1.0, 0.0, 1.0 } } );
        }
        if ( tSets.find( "near_zero" ) != std::string::npos )
        {
            tValueSets.push_back( { "near zero d=1e-13", { -1.0, -1.0e-13, 1.0e-13, 1.0 } } );
            tValueSets.push_back( { "near zero d=1e-11", { -1.0, -1.0e-11, 1.0e-11, 1.0 } } );
        }

        uint tTotalFailing = 0;
        uint tTotalThrew   = 0;
        for ( const std::pair< std::string, Vector< real > >& iValues : tValueSets )
        {
            for ( Regime iRegime : selected_regimes( { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } ) )
            {
                Sweep_Totals tTotals;
                uint         tNumValues = iValues.second.size();
                uint         tNumCases  = 1;
                for ( uint iCorner = 0; iCorner < 8; iCorner++ )
                {
                    tNumCases *= tNumValues;
                }
                Run_Options tOptions;
                for ( uint iCase = 0; iCase < tNumCases; iCase++ )
                {
                    std::array< real, 8 > tCorners;
                    uint                  tCode = iCase;
                    for ( uint iCorner = 0; iCorner < 8; iCorner++ )
                    {
                        tCorners[ iCorner ] = iValues.second( tCode % tNumValues );
                        tCode /= tNumValues;
                    }
                    set_hex_corners( tCorners );
                    std::string  tName    = values_string( gNodeGrid.mValues );
                    Case_Outcome tOutcome = run_case( tHex.mMesh, { tName, { node_grid() }, {} }, iRegime, tOptions );
                    tTotals.add( tName, tOutcome );
                    if ( tOutcome.mThrew or tOutcome.mCounts.failures() > 0 )
                    {
                        write_case( "hex " + iValues.first, iRegime, tName, tOutcome );
                    }
                }
                print_totals( "single hex, corner values " + iValues.first + ", " + regime_name( iRegime ), tTotals );
                tTotalFailing += tTotals.mFailing;
                tTotalThrew += tTotals.mThrew;
            }
        }
        CHECK( tTotalFailing == 0 );
        CHECK( tTotalThrew == 0 );
    }

    TEST_CASE( "XTK connectivity sweep: single hex, two geometries in both orders", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_HEX_TWO]" )
    {
        // the second geometries of [XTK_PATHOLOGY_HEX_TWO_GEOMETRIES], as planes; geometry A from the corner sweep
        if ( par_size() != 1 )
        {
            return;
        }
        Grid_Mesh tHex( { 1, 1, 1 }, { 1.0, 1.0, 1.0 } );
        Vector< std::pair< std::string, Shape > > tSeconds = {
            { "x=0.5", plane( 1.0, 0.0, 0.0, -0.5 ) },
            { "x+y=1", plane( 1.0, 1.0, 0.0, -1.0 ) },
            { "x+y+z=1", plane( 1.0, 1.0, 1.0, -1.0 ) },
            { "x+y+z=1.5", plane( 1.0, 1.0, 1.0, -1.5 ) },
            { "z=0", plane( 0.0, 0.0, 1.0, 0.0 ) },
            { "x=0.5+1e-11", plane( 1.0, 0.0, 0.0, -0.5 - 1.0e-11 ) },
            { "sphere r=0.4", sphere( { 0.5, 0.5, 0.5 }, 0.4 ) }
        };
        // two analytic/multilinear geometries hit MB-040 (heap corruption) without the fix of PR #53: STORED by default
        uint tTotalFailing = 0;
        for ( Regime iRegime : selected_regimes( { Regime::STORED } ) )
        {
            Sweep_Totals tTotals;
            Run_Options  tOptions;
            for ( const std::pair< std::string, Shape >& iSecond : tSeconds )
            {
                for ( uint iCase = 0; iCase < 6561; iCase++ )
                {
                    std::array< real, 8 > tCorners;
                    uint                  tCode = iCase;
                    for ( uint iCorner = 0; iCorner < 8; iCorner++ )
                    {
                        tCorners[ iCorner ] = real( (int)( tCode % 3 ) - 1 );
                        tCode /= 3;
                    }
                    set_hex_corners( tCorners );
                    std::string tName = values_string( gNodeGrid.mValues ) + " | " + iSecond.first;
                    for ( bool tGridFirst : { true, false } )
                    {
                        Vector< Shape > tShapes  = tGridFirst ? Vector< Shape >{ node_grid(), iSecond.second } : Vector< Shape >{ iSecond.second, node_grid() };
                        std::string     tFull    = tName + ( tGridFirst ? " AB" : " BA" );
                        Case_Outcome    tOutcome = run_case( tHex.mMesh, { tFull, tShapes, {} }, iRegime, tOptions );
                        tTotals.add( tFull, tOutcome );
                        if ( tOutcome.mThrew or tOutcome.mCounts.failures() > 0 )
                        {
                            write_case( "hex two geometries", iRegime, tFull, tOutcome );
                        }
                    }
                }
            }
            print_totals( "single hex, two geometries, " + regime_name( iRegime ), tTotals );
            tTotalFailing += tTotals.mFailing;
        }
        CHECK( tTotalFailing == 0 );
    }

    // ----------------------------------------------------------------------------------
    // section 6 P1-P5, P7 (+ E4, E5, E9) on any number of ranks

    TEST_CASE( "XTK connectivity: processor-boundary cases on any number of ranks", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_MPI]" )
    {
        Grid_Mesh tBox( { 8, 4, 4 }, { 2.0, 1.0, 1.0 } );
        print_ownership( tBox.mMesh );
        set_checkerboard_grid( { 8, 4, 4 }, { 0.25, 0.25, 0.25 } );

        Run_Options tOptions;
        tOptions.mParallel = true;
        Sweep_Totals tTotals = run_cases( "mpi", tBox.mMesh, parallel_cases(), { Regime::STORED, Regime::ANALYTIC }, tOptions, true );
        print_totals( "processor-boundary cases P1-P5, E4, E5, E9", tTotals );
        if ( par_rank() == 0 )
        {
            CHECK( tTotals.mThrew == 0 );
            CHECK( tTotals.mFailing == 0 );
        }
    }

    // ----------------------------------------------------------------------------------
    // section 3: gap and kissing sweeps, any number of ranks

    TEST_CASE( "XTK connectivity sweep: gaps between same-phase regions", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_GAP]" )
    {
        Grid_Mesh tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 } );
        print_ownership( tBox.mMesh );
        Run_Options tOptions;
        tOptions.mParallel = true;

        Vector< Connectivity_Case > tCases = gap_cases( 1.0, 0.5, false );
        Vector< Connectivity_Case > tSlanted = gap_cases( 1.0, 0.5, true );
        tCases.append( tSlanted );
        Sweep_Totals tTotals = run_cases( "gap", tBox.mMesh, tCases, { Regime::STORED, Regime::ANALYTIC }, tOptions, true );
        print_totals( "gap sweep", tTotals );
        if ( par_rank() == 0 )
        {
            CHECK( tTotals.mThrew == 0 );
            CHECK( tTotals.mFailing == 0 );
        }
    }

    TEST_CASE( "XTK connectivity sweep: same-phase regions kissing at an edge or a vertex", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_KISSING]" )
    {
        Grid_Mesh tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 } );
        print_ownership( tBox.mMesh );
        Run_Options tOptions;
        tOptions.mParallel = true;
        Sweep_Totals tTotals = run_cases( "kissing", tBox.mMesh, kissing_cases( 1.0, 0.5 ), { Regime::STORED, Regime::ANALYTIC }, tOptions, true );
        print_totals( "kissing sweep", tTotals );
        if ( par_rank() == 0 )
        {
            CHECK( tTotals.mThrew == 0 );
            CHECK( tTotals.mFailing == 0 );
        }
    }

    TEST_CASE( "XTK connectivity sweep: random node signs on a multi-cell grid", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_RANDOM]" )
    {
        // MORIS_XTK_CONNECTIVITY_RANDOM="<cells per side>,<number of cases>,<seed>" (default 2,2000,1)
        uint        tCellsPerSide = 2;
        uint        tNumCases     = 2000;
        uint        tSeed         = 1;
        if ( const char* tSpec = std::getenv( "MORIS_XTK_CONNECTIVITY_RANDOM" ) )
        {
            std::stringstream tStream( tSpec );
            std::string       tItem;
            std::getline( tStream, tItem, ',' );
            tCellsPerSide = std::stoi( tItem );
            std::getline( tStream, tItem, ',' );
            tNumCases = std::stoi( tItem );
            std::getline( tStream, tItem, ',' );
            tSeed = std::stoi( tItem );
        }
        Grid_Mesh tGrid( { tCellsPerSide, tCellsPerSide, tCellsPerSide }, { 1.0, 1.0, 1.0 } );
        print_ownership( tGrid.mMesh );
        Run_Options tOptions;
        tOptions.mParallel = par_size() > 1;
        tOptions.mGhost    = true;

        std::mt19937 tRandom( tSeed );
        real         tSpacing = 1.0 / tCellsPerSide;
        uint         tSkipped = 0;
        for ( Regime iRegime : selected_regimes( { Regime::STORED } ) )
        {
            Sweep_Totals tTotals;
            for ( uint iCase = 0; iCase < tNumCases; iCase++ )
            {
                // the first cases are structured patterns (checkerboards, saddles); every rank draws the same numbers
                bool tUsable = set_random_grid( { tCellsPerSide, tCellsPerSide, tCellsPerSide }, { tSpacing, tSpacing, tSpacing }, tRandom, iCase < 3 ? iCase : 3 );
                if ( not tUsable )
                {
                    tSkipped++;
                    continue;
                }
                std::string  tName    = "case " + std::to_string( iCase ) + " [" + values_string( gNodeGrid.mValues ) + "]";
                Case_Outcome tOutcome = run_case( tGrid.mMesh, { tName, { node_grid() }, {} }, iRegime, tOptions );
                tTotals.add( tName, tOutcome );
                if ( tOptions.mParallel or tOutcome.mThrew or tOutcome.mCounts.failures() > 0 or tOutcome.mCounts.mBodyMismatch > 0 )
                {
                    write_case( "random " + std::to_string( tCellsPerSide ) + "^3", iRegime, tName, tOutcome );
                }
            }
            print_totals( "random node signs on " + std::to_string( tCellsPerSide ) + "^3 cells, " + regime_name( iRegime ) + ", skipped (all-zero cell) "
                                  + std::to_string( tSkipped ),
                    tTotals );
            if ( par_rank() == 0 )
            {
                CHECK( tTotals.mThrew == 0 );
                CHECK( tTotals.mFailing == 0 );
            }
        }
    }

    // ----------------------------------------------------------------------------------
    // section 5 E6 (triangulate_all_in_post), E7 (quadratic IG cells), E8 (discretizations)

    TEST_CASE( "XTK connectivity: triangulate_all_in_post and quadratic IG cells (E6, E7)", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_POST]" )
    {
        Grid_Mesh tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 } );
        print_ownership( tBox.mMesh );
        Vector< Connectivity_Case > tCases = kissing_cases( 1.0, 0.5 );
        Vector< Connectivity_Case > tGaps  = gap_cases( 1.0, 0.5, false );
        tCases.append( tGaps );

        uint tThrew   = 0;
        uint tFailing = 0;
        for ( uint tVariant = 0; tVariant < 2; tVariant++ )
        {
            Run_Options tOptions;
            tOptions.mParallel             = true;
            tOptions.mTriangulateAllInPost = tVariant == 0;
            tOptions.mIgOrder              = tVariant == 0 ? 1 : 2;
            std::string  tName             = tVariant == 0 ? "triangulate_all_in_post (E6)" : "ig_element_order 2 (E7)";
            Sweep_Totals tTotals           = run_cases( tVariant == 0 ? "post E6" : "post E7", tBox.mMesh, tCases, { Regime::STORED }, tOptions, true );
            print_totals( tName, tTotals );
            tThrew += tTotals.mThrew;
            tFailing += tTotals.mFailing;
        }
        if ( par_rank() == 0 )
        {
            CHECK( tThrew == 0 );
            CHECK( tFailing == 0 );
        }
    }

    TEST_CASE( "XTK connectivity: quadratic B-splines, a coarser B-spline mesh, two discretizations (E8)", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_DISCRETIZATION]" )
    {
        Vector< Connectivity_Case > tCases = parallel_cases();
        Vector< Connectivity_Case > tKiss  = kissing_cases( 1.0, 0.25 );
        tCases.append( tKiss );
        set_checkerboard_grid( { 8, 4, 4 }, { 0.25, 0.25, 0.25 } );

        uint tThrew   = 0;
        uint tFailing = 0;
        {
            // p = 2 B-splines on the Lagrange mesh: 3^3-cell supports
            Grid_Mesh tBox( { 8, 4, 4 }, { 2.0, 1.0, 1.0 }, MORIS_REAL_MAX, "2" );
            print_ownership( tBox.mMesh );
            Run_Options tOptions;
            tOptions.mParallel   = true;
            Sweep_Totals tTotals = run_cases( "E8 p=2", tBox.mMesh, tCases, { Regime::STORED }, tOptions, true );
            print_totals( "E8: quadratic B-splines", tTotals );
            tThrew += tTotals.mThrew;
            tFailing += tTotals.mFailing;
        }
        {
            // Lagrange mesh refined once (h = 0.125) with B-spline mesh 0 (p = 1) on it and B-spline mesh 1 (p = 2) on the
            // coarse pattern (h = 0.25, extraction); both are enriched and checked
            Grid_Mesh tBox( { 8, 4, 4 }, { 2.0, 1.0, 1.0 }, MORIS_REAL_MAX, "1,2", "0,1", "0,1", 1 );
            print_ownership( tBox.mMesh );
            Run_Options tOptions;
            tOptions.mParallel    = true;
            tOptions.mMeshIndices = { 0, 1 };
            Sweep_Totals tTotals  = run_cases( "E8 two meshes", tBox.mMesh, tCases, { Regime::STORED }, tOptions, true );
            print_totals( "E8: two discretizations, the second coarser and quadratic", tTotals );
            tThrew += tTotals.mThrew;
            tFailing += tTotals.mFailing;
        }
        if ( par_rank() == 0 )
        {
            CHECK( tThrew == 0 );
            CHECK( tFailing == 0 );
        }
    }

    TEST_CASE( "XTK connectivity: sorted enrichment levels with coincident centroids", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_SORT]" )
    {
        // a solid core inside a void shell inside solid, centred on the basis node ( 0.5, 0.5, 0.5 ) of a 4^3 mesh: the
        // basis support holds two solid levels (core, outside) whose volume-weighted centroids both sit at the node
        if ( par_size() != 1 )
        {
            return;
        }
        Grid_Mesh                   tGrid( { 4, 4, 4 }, { 1.0, 1.0, 1.0 } );
        Vector< Connectivity_Case > tCases;
        for ( real tShift : { 0.0, 1.0e-3 } )
        {
            Point tCentre = { 0.5 + tShift, 0.5, 0.5 };
            tCases.push_back( { "core r=0.12 in a void shell to r=0.2, centre shifted by " + number_string( tShift ),
                    { sphere( tCentre, 0.12 ), sphere( tCentre, 0.2 ) },
                    phase_table( 2, &solid_except_shell ) } );
        }
        for ( bool tSort : { false, true } )
        {
            Run_Options tOptions;
            tOptions.mSortLevels = tSort;
            Sweep_Totals tTotals = run_cases( tSort ? "sorted levels" : "unsorted levels", tGrid.mMesh, tCases, { Regime::STORED, Regime::ANALYTIC }, tOptions, true );
            print_totals( std::string( "concentric spheres, sort_basis_enrichment_levels " ) + ( tSort ? "on" : "off" ), tTotals );
        }
    }

    TEST_CASE( "XTK connectivity: background facet indices at a refinement transition", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_TRANSITION_FACETS]" )
    {
        // Diagnostic for the SIGSEGV of [XTK_CONNECTIVITY_TWO_LEVEL] case A2: for every cell with a neighbour on another
        // level, the facet index the subphase neighbourhood uses (get_elements_connected_to_element_and_face_ind_loc_inds,
        // i.e. get_hmr_facet( ordinal )->get_index()) next to the facet index of the same side from the element-to-facet
        // connectivity (what the facet ancestry and compute_bg_facet_to_child_facet_connectivity key by).
        if ( par_size() != 1 )
        {
            return;
        }
        Grid_Mesh                tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 }, 1.5 );
        mtk::Interpolation_Mesh* tMesh     = tBox.mMesh;
        uint                     tNumFacets = tMesh->get_num_entities( tMesh->get_facet_rank() );
        std::cout << "\n=== facet indices at the refinement transition (mesh facets: " << tNumFacets << ") ===\n";
        uint tDiffer    = 0;
        uint tInvalid   = 0;
        uint tUnflagged = 0;
        for ( uint iCell = 0; iCell < tMesh->get_num_elems(); iCell++ )
        {
            Matrix< IndexMat > tByIndex   = tMesh->get_elements_connected_to_element_and_face_ind_loc_inds( iCell );
            Matrix< IndexMat > tByOrdinal = tMesh->get_elements_connected_to_element_and_face_ord_loc_inds( iCell );
            Matrix< IndexMat > tCellFacets = tMesh->get_entity_connected_to_entity_loc_inds( iCell, mtk::EntityRank::ELEMENT, tMesh->get_facet_rank() );
            for ( uint iN = 0; iN < tByIndex.n_cols(); iN++ )
            {
                moris_index tNeighbor = tByIndex( 0, iN );
                if ( tByOrdinal( 3, iN ) == MORIS_INDEX_MAX )
                {
                    // a neighbour on another level reported WITHOUT a transition location: construct_subphase_neighborhood
                    // then takes the same-level branch and collects the neighbour's subphases with this cell's facet
                    if ( tMesh->get_mtk_cell( tNeighbor ).get_level() != tMesh->get_mtk_cell( iCell ).get_level() )
                    {
                        Matrix< IndexMat > tOtherFacets = tMesh->get_entity_connected_to_entity_loc_inds( tNeighbor, mtk::EntityRank::ELEMENT, tMesh->get_facet_rank() );
                        tUnflagged++;
                        std::cout << "cell " << std::setw( 3 ) << iCell << " (level " << tMesh->get_mtk_cell( iCell ).get_level() << ") neighbour " << std::setw( 3 )
                                  << tNeighbor << " (level " << tMesh->get_mtk_cell( tNeighbor ).get_level() << ")  ords " << tByOrdinal( 1, iN ) << "/"
                                  << tByOrdinal( 2, iN ) << "  NO transition location  facet used " << tByIndex( 1, iN ) << "  neighbour's own facet "
                                  << tOtherFacets( tByOrdinal( 2, iN ) ) << "\n";
                    }
                    continue;
                }
                moris_index tUsedFacet       = tByIndex( 1, iN );
                moris_index tMyOrdinal       = tByOrdinal( 1, iN );
                moris_index tNeighborOrdinal = tByOrdinal( 2, iN );
                moris_index tAncestryFacet   = tCellFacets( tMyOrdinal );
                Matrix< IndexMat > tNeighborFacets = tMesh->get_entity_connected_to_entity_loc_inds( tNeighbor, mtk::EntityRank::ELEMENT, tMesh->get_facet_rank() );
                moris_index tNeighborFacet = tNeighborFacets( tNeighborOrdinal );
                tDiffer += tUsedFacet != tAncestryFacet;
                tInvalid += tUsedFacet < 0 or tUsedFacet >= (moris_index)tNumFacets or tAncestryFacet < 0 or tAncestryFacet >= (moris_index)tNumFacets;
                Box tMine  = cell_box( tMesh->get_mtk_cell( iCell ) );
                std::cout << "cell " << std::setw( 3 ) << iCell << " (level " << tMesh->get_mtk_cell( iCell ).get_level() << ", x " << tMine.mMin[ 0 ] << ".."
                          << tMine.mMax[ 0 ] << ") neighbour " << std::setw( 3 ) << tNeighbor << " (level " << tMesh->get_mtk_cell( tNeighbor ).get_level()
                          << ")  ords " << tMyOrdinal << "/" << tNeighborOrdinal << "  transition " << tByOrdinal( 3, iN ) << "  facet used " << tUsedFacet
                          << "  element-to-facet: mine " << tAncestryFacet << ", neighbour's " << tNeighborFacet << "\n";
            }
        }
        std::cout << "transition sides whose used facet differs from the element-to-facet one: " << tDiffer << ", invalid indices: " << tInvalid
                  << ", cross-level sides without a transition location: " << tUnflagged << "\n";
    }

    // ----------------------------------------------------------------------------------
    // the refinement-transition branch of construct_subphase_neighborhood

    TEST_CASE( "XTK connectivity: refinement transition with two same-phase pieces on the coarse face", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_TWO_LEVEL]" )
    {
        // The box [0,2] x [0,1] x [0,1], 4 x 2 x 2 cells, is refined one level for x > 1.5 (fine h = 0.25). The transition
        // branch (:1472-1510) links every same-phase subphase on the coarse face with every same-phase subphase on the fine
        // face without checking that they touch. It over-links only if both cells at a transition are cut and the coarse
        // face carries two same-phase pieces. The level check (check_intersected_background_cell_levels) refuses any
        // decomposition in which a rank holds cut cells of two levels, and every rank holds the aura cells next to its
        // owned cells, so a coarse-cut / fine-cut transition pair is always refused (route B, on any number of ranks).
        // A cut coarse cell next to uncut fine cells has a single-phase face (the fine cells share the face-centre node),
        // and an uncut cell is a single subphase, so the reachable configurations (route A) cannot over-link either.
        Grid_Mesh tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 }, 1.5 );
        print_ownership( tBox.mMesh );

        Vector< Connectivity_Case > tCases;
        // route A: the transition branch runs, one side uncut
        tCases.push_back( { "A1 sphere in a coarse cell next to the fine region", { sphere( { 1.25, 0.25, 0.25 }, 0.15 ) }, {} } );
        // A2 crashes (SIGSEGV in collect_subphases_attached_to_facet_on_cell, cl_XTK_Integration_Mesh_Generator.cpp:1637):
        // HMR reports the fine -> coarse side without a transition location, and the same-level branch reads the coarse
        // cut cell's child facets under the fine cell's facet (null). Run it with MORIS_XTK_CONNECTIVITY_CRASH_CASES=1.
        if ( std::getenv( "MORIS_XTK_CONNECTIVITY_CRASH_CASES" ) != nullptr )
        {
            tCases.push_back( { "A2 plane x=1.3: cut coarse cells next to uncut fine", { plane( 1.0, 0.0, 0.0, -1.3 ) }, {} } );
        }
        else if ( par_rank() == 0 )
        {
            std::cout << "case A2 (cut coarse cells next to uncut fine cells) skipped: known SIGSEGV, set MORIS_XTK_CONNECTIVITY_CRASH_CASES=1\n";
        }
        tCases.push_back( { "A3 sphere in a fine cell next to the coarse region", { sphere( { 1.62, 0.37, 0.37 }, 0.1 ) }, {} } );
        tCases.push_back( { "A4 two fine prongs ending on the transition face",
                { box_union( { { 1.5, 1.9, 0.2, 0.3, 0.2, 0.8 }, { 1.5, 1.9, 0.7, 0.8, 0.2, 0.8 } } ) }, {} } );
        tCases.push_back( { "A5 plane x=1.7: cut fine cells next to uncut coarse", { plane( 1.0, 0.0, 0.0, -1.7 ) }, {} } );
        // route B: cut cells on both levels (refused unless each rank sees one level, MB-044)
        tCases.push_back( { "B1 coarse-cut / fine-cut transition pair (slanted plane across x=1.5)", { plane_through( { 1.0, 0.3, 0.2 }, { 1.5, 0.5, 0.5 }, 0.0 ) }, {} } );
        tCases.push_back( { "B2 prongs crossing the transition face",
                { box_union( { { 1.1, 1.8, 0.2, 0.3, 0.2, 0.8 }, { 1.1, 1.8, 0.7, 0.8, 0.2, 0.8 } } ) }, {} } );
        tCases.push_back( { "B3 x=0.3 coarse and x=1.8 fine (MB-044 layout)", { plane( 1.0, 0.0, 0.0, -0.3 ), plane( 1.0, 0.0, 0.0, -1.8 ) },
                phase_table( 2, &solid_if_all_same_sign ) } );

        Run_Options tOptions;
        tOptions.mParallel = true;
        // ANALYTIC crashes in case A1 by the same mechanism as A2: the sphere inside the coarse cell is detected through the
        // subdivision vertices (MB-031 test), so a coarse cell next to the uncut fine region is cut. Opt in as for A2.
        Vector< Regime > tRegimes = { Regime::STORED };
        if ( std::getenv( "MORIS_XTK_CONNECTIVITY_CRASH_CASES" ) != nullptr )
        {
            tRegimes.push_back( Regime::ANALYTIC );
        }
        Sweep_Totals tTotals = run_cases( "two-level", tBox.mMesh, tCases, tRegimes, tOptions, true );
        print_totals( "refinement transition", tTotals );
        // characterization: no CHECK on the reference counts; the verdict is read from the output
        if ( par_rank() == 0 )
        {
            CHECK( tTotals.mThrew == 0 );
        }
    }

}    // namespace moris::xtk
