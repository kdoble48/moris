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


    // ----------------------------------------------------------------------------------
    // refinement-transition lineup: the background neighbour records (facet, side ordinals, transition location) and the
    // local coordinates of every double-sided cluster side, checked against the physical geometry

    namespace
    {
        Box
        vertex_box( const Vector< mtk::Vertex const * >& aVertices )
        {
            Box tBox;
            for ( mtk::Vertex const * iVertex : aVertices )
            {
                Matrix< DDRMat > tX = iVertex->get_coords();
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tBox.mMin[ iDim ] = std::min( tBox.mMin[ iDim ], tX( iDim ) );
                    tBox.mMax[ iDim ] = std::max( tBox.mMax[ iDim ], tX( iDim ) );
                }
            }
            return tBox;
        }

        bool
        same_box( const Box& aA, const Box& aB, real aTolerance )
        {
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                if ( std::abs( aA.mMin[ iDim ] - aB.mMin[ iDim ] ) > aTolerance or std::abs( aA.mMax[ iDim ] - aB.mMax[ iDim ] ) > aTolerance )
                {
                    return false;
                }
            }
            return true;
        }

        bool
        box_inside( const Box& aInner, const Box& aOuter, real aTolerance )
        {
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                if ( aInner.mMin[ iDim ] < aOuter.mMin[ iDim ] - aTolerance or aInner.mMax[ iDim ] > aOuter.mMax[ iDim ] + aTolerance )
                {
                    return false;
                }
            }
            return true;
        }

        // the axis along which a face box is flat
        uint
        flat_axis( const Box& aFace )
        {
            uint tAxis = 0;
            for ( uint iDim = 1; iDim < 3; iDim++ )
            {
                if ( aFace.size( iDim ) < aFace.size( tAxis ) )
                {
                    tAxis = iDim;
                }
            }
            return tAxis;
        }

        // the column of aCell's neighbour records that points at aNeighbor (MORIS_INDEX_MAX if none)
        moris_index
        record_of_neighbor( const Matrix< IndexMat >& aRecords, moris_index aNeighbor )
        {
            for ( uint iN = 0; iN < aRecords.n_cols(); iN++ )
            {
                if ( aRecords( 0, iN ) == aNeighbor )
                {
                    return (moris_index)iN;
                }
            }
            return MORIS_INDEX_MAX;
        }

        struct Lineup_Counts
        {
            uint mSameLevelSides  = 0;
            uint mCoarseToFine    = 0;
            uint mFineToCoarse    = 0;
            uint mLevelJumps      = 0;    // neighbours more than one level apart
            uint mBadSameLevel    = 0;    // the two faces or the two facet indices differ
            uint mBadCoarseToFine = 0;    // fine face not on the coarse face, no transition location, or facet not the coarse cell's own
            uint mBadLocation     = 0;    // one (side ordinal, transition location) names two different positions on the face
            uint mDuplicateLocation = 0;  // two fine neighbours on one coarse side with the same transition location
            uint mBadFineToCoarse = 0;    // the neighbour ordinal does not name the coarse face that contains mine, or its facet differs
            uint mFacetNotNeighbors = 0;  // fine -> coarse: the facet index HMR gives is not the coarse cell's facet (MB-075)
            uint mClusterSides    = 0;
            uint mCrossLevelSides = 0;
            uint mBadLocalCoords  = 0;    // a side vertex's local coordinate, mapped through the interpolation cell, misses the vertex
            uint mBadPatch        = 0;    // leader and follower sides do not cover the same patch
            uint mCrossLevelGhost = 0;
            std::array< uint, 6 > mCrossLevelGhostByDirection = {};    // -x, +x, -y, +y, -z, +z of the coarse side's face
            Vector< std::string > mNotes;

            void
            note( const std::string& aNote )
            {
                if ( mNotes.size() < 12 )
                {
                    mNotes.push_back( aNote );
                }
            }

            void
            reduce()
            {
                for ( uint* iValue : { &mSameLevelSides, &mCoarseToFine, &mFineToCoarse, &mLevelJumps, &mBadSameLevel, &mBadCoarseToFine, &mBadLocation,
                              &mDuplicateLocation, &mBadFineToCoarse, &mFacetNotNeighbors, &mClusterSides, &mCrossLevelSides, &mBadLocalCoords,
                              &mBadPatch, &mCrossLevelGhost } )
                {
                    *iValue = sum_all( *iValue );
                }
                for ( uint& iValue : mCrossLevelGhostByDirection )
                {
                    iValue = sum_all( iValue );
                }
            }
        };

        // (1) background neighbour records against the geometry of the two cells
        void
        check_background_lineup( mtk::Interpolation_Mesh* aMesh, Lineup_Counts& aCounts )
        {
            // (side ordinal, transition location) -> position of the fine face's centre on the coarse face, in units of 1e-6
            // of the coarse face size; HMR must use one convention for all cells
            std::map< std::pair< moris_index, moris_index >, std::pair< std::int64_t, std::int64_t > > tPositionOfLocation;

            for ( uint iCell = 0; iCell < aMesh->get_num_elems(); iCell++ )
            {
                const mtk::Cell&   tCell      = aMesh->get_mtk_cell( iCell );
                real               tTolerance = 1.0e-9 * cell_box( tCell ).size( 0 );
                Matrix< IndexMat > tByIndex   = aMesh->get_elements_connected_to_element_and_face_ind_loc_inds( iCell );
                Matrix< IndexMat > tByOrdinal = aMesh->get_elements_connected_to_element_and_face_ord_loc_inds( iCell );
                Matrix< IndexMat > tMyFacets  = aMesh->get_entity_connected_to_entity_loc_inds( iCell, mtk::EntityRank::ELEMENT, aMesh->get_facet_rank() );
                std::map< moris_index, std::set< moris_index > > tLocationsOnSide;

                for ( uint iN = 0; iN < tByOrdinal.n_cols(); iN++ )
                {
                    moris_index      tNeighborIndex = tByOrdinal( 0, iN );
                    const mtk::Cell& tNeighbor      = aMesh->get_mtk_cell( tNeighborIndex );
                    moris_index      tMyOrdinal     = tByOrdinal( 1, iN );
                    moris_index      tNbOrdinal     = tByOrdinal( 2, iN );
                    moris_index      tLocation      = tByOrdinal( 3, iN );
                    moris_index      tGivenFacet    = tByIndex( 1, iN );
                    Box              tMyFace        = vertex_box( tCell.get_vertices_on_side_ordinal( tMyOrdinal ) );
                    Box              tNbFace        = vertex_box( tNeighbor.get_vertices_on_side_ordinal( tNbOrdinal ) );
                    uint             tAxis          = flat_axis( tMyFace );
                    bool             tCoplanar      = flat_axis( tNbFace ) == tAxis and std::abs( tMyFace.mMin[ tAxis ] - tNbFace.mMin[ tAxis ] ) <= tTolerance;
                    Matrix< IndexMat > tNbFacets    = aMesh->get_entity_connected_to_entity_loc_inds( tNeighborIndex, mtk::EntityRank::ELEMENT, aMesh->get_facet_rank() );
                    Matrix< IndexMat > tNbByIndex   = aMesh->get_elements_connected_to_element_and_face_ind_loc_inds( tNeighborIndex );
                    moris_index      tBackRecord    = record_of_neighbor( tNbByIndex, (moris_index)iCell );
                    uint             tMyLevel       = tCell.get_level();
                    uint             tNbLevel       = tNeighbor.get_level();
                    std::string      tWhere         = "cell " + std::to_string( iCell ) + " (level " + std::to_string( tMyLevel ) + ") -> " + std::to_string( tNeighborIndex )
                                         + " (level " + std::to_string( tNbLevel ) + "), ords " + std::to_string( tMyOrdinal ) + "/" + std::to_string( tNbOrdinal )
                                         + ", location " + std::to_string( tLocation );

                    if ( tMyLevel == tNbLevel )
                    {
                        aCounts.mSameLevelSides++;
                        bool tGood = tCoplanar and same_box( tMyFace, tNbFace, tTolerance ) and tLocation == MORIS_INDEX_MAX and tGivenFacet == tMyFacets( tMyOrdinal )
                                 and tGivenFacet == tNbFacets( tNbOrdinal ) and tBackRecord != MORIS_INDEX_MAX and tNbByIndex( 1, tBackRecord ) == tGivenFacet;
                        if ( not tGood )
                        {
                            aCounts.mBadSameLevel++;
                            aCounts.note( "same level: " + tWhere );
                        }
                    }
                    else if ( tNbLevel == tMyLevel + 1 )
                    {
                        // coarse -> fine: the fine face lies on my face, at the position the transition location names
                        aCounts.mCoarseToFine++;
                        bool tGood = tCoplanar and box_inside( tNbFace, tMyFace, tTolerance ) and tLocation != MORIS_INDEX_MAX and tGivenFacet == tMyFacets( tMyOrdinal );
                        if ( not tGood )
                        {
                            aCounts.mBadCoarseToFine++;
                            aCounts.note( "coarse -> fine: " + tWhere );
                            continue;
                        }
                        if ( not tLocationsOnSide[ tMyOrdinal ].insert( tLocation ).second )
                        {
                            aCounts.mDuplicateLocation++;
                            aCounts.note( "duplicate transition location: " + tWhere );
                        }
                        uint                                  tAxis1 = ( tAxis + 1 ) % 3;
                        uint                                  tAxis2 = ( tAxis + 2 ) % 3;
                        std::pair< std::int64_t, std::int64_t > tPosition = {
                            std::llround( 1.0e6 * ( 0.5 * ( tNbFace.mMin[ tAxis1 ] + tNbFace.mMax[ tAxis1 ] ) - tMyFace.mMin[ tAxis1 ] ) / tMyFace.size( tAxis1 ) ),
                            std::llround( 1.0e6 * ( 0.5 * ( tNbFace.mMin[ tAxis2 ] + tNbFace.mMax[ tAxis2 ] ) - tMyFace.mMin[ tAxis2 ] ) / tMyFace.size( tAxis2 ) )
                        };
                        auto tKnown = tPositionOfLocation.insert( { { tMyOrdinal, tLocation }, tPosition } );
                        if ( not tKnown.second and tKnown.first->second != tPosition )
                        {
                            aCounts.mBadLocation++;
                            aCounts.note( "transition location names another position: " + tWhere );
                        }
                    }
                    else if ( tMyLevel == tNbLevel + 1 )
                    {
                        // fine -> coarse: the neighbour ordinal names the coarse face that contains mine, and that face's facet is
                        // the one the coarse cell reports for this pair
                        aCounts.mFineToCoarse++;
                        bool tGood = tCoplanar and box_inside( tMyFace, tNbFace, tTolerance ) and tBackRecord != MORIS_INDEX_MAX
                                 and tNbByIndex( 1, tBackRecord ) == tNbFacets( tNbOrdinal );
                        if ( not tGood )
                        {
                            aCounts.mBadFineToCoarse++;
                            aCounts.note( "fine -> coarse: " + tWhere );
                        }
                        if ( tGivenFacet != tNbFacets( tNbOrdinal ) )
                        {
                            aCounts.mFacetNotNeighbors++;
                        }
                    }
                    else
                    {
                        aCounts.mLevelJumps++;
                        aCounts.note( "level jump: " + tWhere );
                    }
                }
            }
        }

        // area and centroid of a planar polygon given by its corner vertices in cyclic order
        void
        add_polygon( const Vector< mtk::Vertex const * >& aCorners, real& aArea, Point& aMoment )
        {
            Matrix< DDRMat > tX0 = aCorners( 0 )->get_coords();
            for ( uint iV = 1; iV + 1 < aCorners.size(); iV++ )
            {
                Matrix< DDRMat > tX1 = aCorners( iV )->get_coords();
                Matrix< DDRMat > tX2 = aCorners( iV + 1 )->get_coords();
                real             tA[ 3 ], tB[ 3 ];
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tA[ iDim ] = tX1( iDim ) - tX0( iDim );
                    tB[ iDim ] = tX2( iDim ) - tX0( iDim );
                }
                real tCross[ 3 ] = { tA[ 1 ] * tB[ 2 ] - tA[ 2 ] * tB[ 1 ], tA[ 2 ] * tB[ 0 ] - tA[ 0 ] * tB[ 2 ], tA[ 0 ] * tB[ 1 ] - tA[ 1 ] * tB[ 0 ] };
                real tArea       = 0.5 * std::sqrt( tCross[ 0 ] * tCross[ 0 ] + tCross[ 1 ] * tCross[ 1 ] + tCross[ 2 ] * tCross[ 2 ] );
                aArea += tArea;
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    aMoment[ iDim ] += tArea * ( tX0( iDim ) + tX1( iDim ) + tX2( iDim ) ) / 3.0;
                }
            }
        }

        // (2) every side of every double-sided cluster: the local coordinates the integration reads, mapped through the
        // interpolation cell, must land on the side's vertices, and the leader and follower sides must cover the same patch
        void
        check_cluster_lineup( Model& aModel, Lineup_Counts& aCounts )
        {
            Enriched_Integration_Mesh& tEnrIgMesh = aModel.get_enriched_integ_mesh( 0 );
            for ( uint iSet = 0; iSet < tEnrIgMesh.get_num_sets(); iSet++ )
            {
                mtk::Set* tSet = tEnrIgMesh.get_set_by_index( iSet );
                if ( tSet->get_set_type() != mtk::SetType::DOUBLE_SIDED_SIDESET )
                {
                    continue;
                }
                bool tGhost = tSet->get_set_name().rfind( "ghost", 0 ) == 0;
                for ( const mtk::Cluster* iCluster : tSet->get_clusters_on_set() )
                {
                    const mtk::Cell& tLeaderIp   = iCluster->get_interpolation_cell( mtk::Leader_Follower::LEADER );
                    const mtk::Cell& tFollowerIp = iCluster->get_interpolation_cell( mtk::Leader_Follower::FOLLOWER );
                    bool             tCross      = tLeaderIp.get_level() != tFollowerIp.get_level();
                    real             tArea[ 2 ]  = { 0.0, 0.0 };
                    Point            tMoment[ 2 ] = { { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 } };
                    real             tSize       = std::min( cell_box( tLeaderIp ).size( 0 ), cell_box( tFollowerIp ).size( 0 ) );
                    std::string      tWhere      = tSet->get_set_name() + " leader cell " + std::to_string( tLeaderIp.get_id() ) + " (level "
                                         + std::to_string( tLeaderIp.get_level() ) + "), follower cell " + std::to_string( tFollowerIp.get_id() ) + " (level "
                                         + std::to_string( tFollowerIp.get_level() ) + ")";

                    for ( uint iSide = 0; iSide < 2; iSide++ )
                    {
                        mtk::Leader_Follower tSide   = iSide == 0 ? mtk::Leader_Follower::LEADER : mtk::Leader_Follower::FOLLOWER;
                        Box                  tIpBox  = cell_box( iCluster->get_interpolation_cell( tSide ) );
                        const Vector< mtk::Cell const * >& tCells = iCluster->get_primary_cells_in_cluster( tSide );
                        for ( uint iCell = 0; iCell < tCells.size(); iCell++ )
                        {
                            aCounts.mClusterSides++;
                            aCounts.mCrossLevelSides += tCross;
                            moris_index                    tOrdinal = iCluster->get_cell_side_ordinal( iCell, tSide );
                            Vector< mtk::Vertex const * > tVertices = tCells( iCell )->get_vertices_on_side_ordinal( tOrdinal );
                            Matrix< DDRMat >               tXi      = iCluster->get_cell_local_coords_on_side_wrt_interp_cell( iCell, tSide );
                            bool                           tGood    = tXi.n_rows() == tVertices.size() and tXi.n_cols() == 3;
                            for ( uint iV = 0; tGood and iV < tVertices.size(); iV++ )
                            {
                                Matrix< DDRMat > tX = tVertices( iV )->get_coords();
                                for ( uint iDim = 0; iDim < 3; iDim++ )
                                {
                                    real tMapped = tIpBox.mMin[ iDim ] + 0.5 * ( tXi( iV, iDim ) + 1.0 ) * tIpBox.size( iDim );
                                    tGood        = tGood and std::abs( tMapped - tX( iDim ) ) <= 1.0e-9 * tSize;
                                }
                            }
                            if ( not tGood )
                            {
                                aCounts.mBadLocalCoords++;
                                aCounts.note( std::string( iSide == 0 ? "leader" : "follower" ) + " local coordinates off the side vertices: " + tWhere );
                            }
                            add_polygon( tCells( iCell )->get_vertices_on_side_ordinal( tOrdinal ), tArea[ iSide ], tMoment[ iSide ] );
                        }
                    }

                    bool tSamePatch = std::abs( tArea[ 0 ] - tArea[ 1 ] ) <= 1.0e-9 * tSize * tSize and tArea[ 0 ] > 0.0;
                    for ( uint iDim = 0; tSamePatch and iDim < 3; iDim++ )
                    {
                        tSamePatch = std::abs( tMoment[ 0 ][ iDim ] / tArea[ 0 ] - tMoment[ 1 ][ iDim ] / tArea[ 1 ] ) <= 1.0e-9 * tSize;
                    }
                    if ( not tSamePatch )
                    {
                        aCounts.mBadPatch++;
                        aCounts.note( "leader and follower sides cover different patches: " + tWhere );
                    }

                    if ( tGhost and tCross )
                    {
                        // direction of the shared face as seen from the coarser side
                        aCounts.mCrossLevelGhost++;
                        const mtk::Cell& tCoarse  = tLeaderIp.get_level() < tFollowerIp.get_level() ? tLeaderIp : tFollowerIp;
                        Box              tCoarseBox = cell_box( tCoarse );
                        Point            tCentre   = { tMoment[ 0 ][ 0 ] / tArea[ 0 ], tMoment[ 0 ][ 1 ] / tArea[ 0 ], tMoment[ 0 ][ 2 ] / tArea[ 0 ] };
                        for ( uint iDim = 0; iDim < 3; iDim++ )
                        {
                            if ( std::abs( tCentre[ iDim ] - tCoarseBox.mMin[ iDim ] ) <= 1.0e-9 * tSize )
                            {
                                aCounts.mCrossLevelGhostByDirection[ 2 * iDim ]++;
                            }
                            else if ( std::abs( tCentre[ iDim ] - tCoarseBox.mMax[ iDim ] ) <= 1.0e-9 * tSize )
                            {
                                aCounts.mCrossLevelGhostByDirection[ 2 * iDim + 1 ]++;
                            }
                        }
                    }
                }
            }
        }

        void
        print_lineup( const std::string& aCase, const Lineup_Counts& aCounts, const Vector< std::string >& aLocalNotes )
        {
            if ( par_rank() == 0 )
            {
                const std::array< uint, 6 >& tDir = aCounts.mCrossLevelGhostByDirection;
                std::cout << "\n=== " << aCase << " (" << par_size() << " rank(s)) ===\n"
                          << "background sides: same level " << aCounts.mSameLevelSides << " (bad " << aCounts.mBadSameLevel << "), coarse->fine "
                          << aCounts.mCoarseToFine << " (bad " << aCounts.mBadCoarseToFine << ", location mismatch " << aCounts.mBadLocation << ", duplicate location "
                          << aCounts.mDuplicateLocation << "), fine->coarse " << aCounts.mFineToCoarse << " (bad " << aCounts.mBadFineToCoarse
                          << "; facet given is not the coarse cell's: " << aCounts.mFacetNotNeighbors << "), level jumps " << aCounts.mLevelJumps << "\n"
                          << "cluster sides: " << aCounts.mClusterSides << " (cross-level " << aCounts.mCrossLevelSides << "), bad local coordinates "
                          << aCounts.mBadLocalCoords << ", leader/follower patch mismatch " << aCounts.mBadPatch << "\n"
                          << "cross-level ghost facets: " << aCounts.mCrossLevelGhost << "  by coarse face -x/+x/-y/+y/-z/+z: " << tDir[ 0 ] << "/" << tDir[ 1 ] << "/"
                          << tDir[ 2 ] << "/" << tDir[ 3 ] << "/" << tDir[ 4 ] << "/" << tDir[ 5 ] << "\n";
            }
            for ( const std::string& iNote : aLocalNotes )
            {
                std::cout << "  [rank " << par_rank() << "] " << iNote << "\n";
            }
        }

        Lineup_Counts
        run_lineup_case( mtk::Interpolation_Mesh* aMesh, const std::string& aCase, const Shape& aShape )
        {
            Lineup_Counts tCounts;
            check_background_lineup( aMesh, tCounts );

            gShapes[ 0 ] = aShape;
            gen::Geometry_Engine_Parameters tParameters;
            tParameters.mGeometries = { create_geometry( 0, Regime::STORED ) };
            gen::Geometry_Engine tGeometryEngine( aMesh, tParameters );
            Model                tModel( 3, aMesh, &tGeometryEngine );
            tModel.mVerbose  = false;
            bool tDecomposed = tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );
            REQUIRE( sum_all( (uint)tDecomposed ) == (uint)par_size() );
            tModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, Matrix< IndexMat >( 1, 1, 0 ), false );
            tModel.construct_face_oriented_ghost_penalization_cells();
            check_cluster_lineup( tModel, tCounts );

            Vector< std::string > tLocalNotes = tCounts.mNotes;
            tCounts.reduce();
            print_lineup( aCase, tCounts, tLocalNotes );
            return tCounts;
        }

        void
        check_lineup_counts( const Lineup_Counts& aCounts )
        {
            CHECK( aCounts.mBadSameLevel == 0 );
            CHECK( aCounts.mBadCoarseToFine == 0 );
            CHECK( aCounts.mBadLocation == 0 );
            CHECK( aCounts.mDuplicateLocation == 0 );
            CHECK( aCounts.mBadFineToCoarse == 0 );
            CHECK( aCounts.mLevelJumps == 0 );
            CHECK( aCounts.mBadLocalCoords == 0 );
            CHECK( aCounts.mBadPatch == 0 );
        }
    }    // namespace

    TEST_CASE( "XTK refinement transition lineup: facets, ordinals, transition locations and ghost local coordinates", "[XTK][XTK_CONNECTIVITY][XTK_TRANSITION_LINEUP]" )
    {
        // A 6 x 6 x 6 box on [0,3]^3 whose central 2 x 2 x 2 block ([1,2]^3) is refined one level, so every face direction is a
        // refinement transition. The geometries cut only fine cells next to the transition faces (a cut coarse cell next to
        // uncut fine cells is MB-075), so legacy ghost builds non-trivial (coarse -> fine) facets on all six sides.
        auto tCentralBlock = []( const Matrix< DDRMat >& aCentre, uint aPass ) -> bool {
            real tLow  = aPass == 0 ? 1.0 : 1.25;
            real tHigh = aPass == 0 ? 2.0 : 1.75;
            return aCentre( 0 ) > tLow and aCentre( 0 ) < tHigh and aCentre( 1 ) > tLow and aCentre( 1 ) < tHigh and aCentre( 2 ) > tLow and aCentre( 2 ) < tHigh;
        };

        SECTION( "one refinement level" )
        {
            Grid_Mesh tBox( { 6, 6, 6 }, { 3.0, 3.0, 3.0 }, tCentralBlock, 1 );
            for ( const std::pair< std::string, Shape >& iCase : Vector< std::pair< std::string, Shape > >{
                          { "sphere r=0.45 in the refined block", sphere( { 1.5, 1.5, 1.5 }, 0.45 ) },
                          { "box [1.1,1.9]^3 in the refined block", box_union( { { 1.1, 1.9, 1.1, 1.9, 1.1, 1.9 } } ) } } )
            {
                Lineup_Counts tCounts = run_lineup_case( tBox.mMesh, iCase.first, iCase.second );
                check_lineup_counts( tCounts );
                CHECK( tCounts.mCoarseToFine > 0 );
                CHECK( tCounts.mFineToCoarse == tCounts.mCoarseToFine );
                for ( uint iDirection = 0; iDirection < 6; iDirection++ )
                {
                    CHECK( tCounts.mCrossLevelGhostByDirection[ iDirection ] > 0 );
                }
            }
        }

        SECTION( "two refinement levels, nested" )
        {
            // second pass: [1.25,1.75]^3 inside the refined block; transitions 0 -> 1 outside and 1 -> 2 inside
            Grid_Mesh     tBox( { 6, 6, 6 }, { 3.0, 3.0, 3.0 }, tCentralBlock, 2 );
            Lineup_Counts tCounts = run_lineup_case( tBox.mMesh, "sphere r=0.2 in the level-2 block", sphere( { 1.5, 1.5, 1.5 }, 0.2 ) );
            check_lineup_counts( tCounts );
            CHECK( tCounts.mCoarseToFine > 0 );
            CHECK( tCounts.mCrossLevelGhost > 0 );
        }
    }

    // ----------------------------------------------------------------------------------
    // the checker as a standalone XTK tool on a model built through the parameter list

    TEST_CASE( "XTK connectivity checker on a model run through the parameter list", "[XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_CHECK]" )
    {
        // a sphere cutting several cells of a 4 x 4 x 4 box; decomposition, enrichment and ghost run through Model::perform()
        Grid_Mesh tBox( { 4, 4, 4 }, { 1.0, 1.0, 1.0 } );

        gShapes[ 0 ] = sphere( { 0.5, 0.5, 0.5 }, 0.3 );
        gen::Geometry_Engine_Parameters tGenParameters;
        tGenParameters.mGeometries = { create_geometry( 0, Regime::STORED ) };
        gen::Geometry_Engine tGeometryEngine( tBox.mMesh, tGenParameters );

        Parameter_List tParameters = prm::create_xtk_parameter_list();
        tParameters.set( "ghost_stab", true );
        Model tModel( tParameters );
        tModel.set_geometry_engine( &tGeometryEngine );
        tModel.set_mtk_background_mesh( tBox.mMesh );
        tModel.set_output_performer( std::make_shared< mtk::Mesh_Manager >() );    // perform_enrichment registers the enriched meshes there
        REQUIRE( tModel.perform() );

        Connectivity_Checker tChecker( tBox.mMesh, &tModel, tGeometryEngine.get_num_bulk_phase(), 0 );
        Connectivity_Counts  tCounts = tChecker.check_and_report( { 0 }, true, true );

        CHECK( tCounts.failures() == 0 );
        CHECK( tCounts.mCutCells > 0 );
        CHECK( tCounts.mLinks > 0 );
        CHECK( tCounts.mLinks == tCounts.mReferenceLinks );
        CHECK( tCounts.mEnrichedBases > 0 );
        CHECK( tCounts.mEnrichedBases == tCounts.mReferenceComponents );
        CHECK( tCounts.mGhostPairs > 0 );
        CHECK( tCounts.sparse_report().rfind( "PASS cut_cells=", 0 ) == 0 );

        // the detailed mode checks the same and reports every count
        Connectivity_Checker tDetailedChecker( tBox.mMesh, &tModel, tGeometryEngine.get_num_bulk_phase(), 0 );
        Connectivity_Counts  tDetailed = tDetailedChecker.check_and_report( { 0 }, true, true, Report_Mode::DETAILED );
        CHECK( tDetailed.failures() == 0 );
        CHECK( tDetailed.mLinks == tCounts.mLinks );
        CHECK( tDetailed.mEnrichedBases == tCounts.mEnrichedBases );
        CHECK( tDetailed.detailed_report().find( "neighbourhood: links" ) != std::string::npos );

        // the sparse line of a failing check: the total and every nonzero failure count, as key=value tokens
        Connectivity_Counts tFailing;
        tFailing.mMissingLinks = 2;
        tFailing.mCrossTalk    = 1;
        CHECK( tFailing.sparse_report() == "FAIL failures=3 missing_links=2 cross_talk=1 (details: Report_Mode::DETAILED)" );
    }
}    // namespace moris::xtk
