/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Cut_Pathologies_MPI.cpp
 *
 * A multi-rank sweep ([XTK_PATHOLOGY_MPI]) decomposes a 4 x 2 x 2 box [0,2] x [0,1] x [0,1] with planes on, near
 * and across the processor boundary (x = 1 on 2 ranks). Besides the global invariants (total volume, negative volume
 * against the exact value, the same result on every number of ranks) it gathers every integration vertex on rank 0
 * and checks that the ranks agree on it: one id per location, and per id the same coordinates, owner and geometric
 * region. Run it on 1, 2 and 4 ranks, e.g. mpirun -n 2 XTK-test.exe "[XTK_PATHOLOGY_MPI]", and compare the CSV files
 * (MORIS_XTK_PATHOLOGY_MPI_CSV=<file>).
 *
 * The sweeps are hidden ([.]) because they characterize current behaviour rather than gate it. Run one with
 *   XTK-test.exe "[XTK_PATHOLOGY_HEX_SIGNS]"
 * Set MORIS_XTK_PATHOLOGY_CSV=<file> to write one row per decomposed case.
 *
 * The helpers are in cl_XTK_Cut_Pathology_MPI_Harness.hpp.
 *
 */

#include "catch.hpp"

#include "cl_XTK_Cut_Pathology_MPI_Harness.hpp"

namespace moris::xtk
{
    using namespace cut_pathology_test;

    TEST_CASE( "XTK pathology sweep: planes at the processor boundary of a box, any number of ranks", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI]" )
    {
        Grid_Mesh                tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 } );
        mtk::Interpolation_Mesh* tMesh = tBox.mMesh;

        // where the processor boundaries are
        Owned_Region tOwned( tMesh );
        real         tOwnedMin[ 3 ] = { MORIS_REAL_MAX, MORIS_REAL_MAX, MORIS_REAL_MAX };
        real         tOwnedMax[ 3 ] = { -MORIS_REAL_MAX, -MORIS_REAL_MAX, -MORIS_REAL_MAX };
        for ( const std::array< real, 6 >& iBox : tOwned.mBoxes )
        {
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tOwnedMin[ iDim ] = std::min( tOwnedMin[ iDim ], iBox[ iDim ] );
                tOwnedMax[ iDim ] = std::max( tOwnedMax[ iDim ], iBox[ 3 + iDim ] );
            }
        }
        std::cout << "rank " << par_rank() << " of " << par_size() << " owns " << tOwned.mBoxes.size() << " background cells in ["
                  << tOwnedMin[ 0 ] << "," << tOwnedMax[ 0 ] << "] x [" << tOwnedMin[ 1 ] << "," << tOwnedMax[ 1 ] << "] x ["
                  << tOwnedMin[ 2 ] << "," << tOwnedMax[ 2 ] << "]\n"
                  << std::flush;

        uint tNumThrew = 0, tNumBadVolume = 0, tNumReference = 0, tNumCellDefects = 0, tNumDisagreements = 0;
        if ( par_rank() == 0 )
        {
            std::cout << "\n=== planes at the processor boundary, " << par_size() << " rank(s) ===\n"
                      << std::setw( 12 ) << "regime" << std::setw( 28 ) << "case" << std::setw( 7 ) << "cells" << std::setw( 6 ) << "subph"
                      << std::setw( 12 ) << "min cell V" << std::setw( 6 ) << "V<=0" << std::setw( 7 ) << "ifacet" << std::setw( 12 ) << "ifacet A"
                      << std::setw( 12 ) << "total V-2" << std::setw( 12 ) << "neg V err" << std::setw( 7 ) << "ids" << std::setw( 7 ) << "shared"
                      << std::setw( 7 ) << "coord" << std::setw( 7 ) << "owner" << std::setw( 7 ) << "label" << std::setw( 7 ) << "dupl" << "\n";
        }

        // MORIS_XTK_PATHOLOGY_MPI_CASE=<n> runs only the n-th case (0-based, in every regime), so a crash can be located
        const char* tOnlyCase = std::getenv( "MORIS_XTK_PATHOLOGY_MPI_CASE" );

        for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
        {
            Vector< Parallel_Case > tCases = parallel_cases();
            for ( uint iCaseIndex = 0; iCaseIndex < tCases.size(); iCaseIndex++ )
            {
                const Parallel_Case& iCase = tCases( iCaseIndex );
                if ( tOnlyCase != nullptr and std::atoi( tOnlyCase ) != (int)iCaseIndex )
                {
                    continue;
                }

                // two analytic/multilinear geometries can hit MB-040 (heap corruption) on develop; STORED only
                if ( iCase.mPlanes.size() > 1 and iRegime != Regime::STORED )
                {
                    continue;
                }

                Parallel_Result tResult = decompose_parallel_case( tMesh, iCase.mPlanes, iRegime );

                Vector< real > tReference;
                real           tReferenceError = 0.0;
                for ( uint iGeom = 0; iGeom < iCase.mPlanes.size(); iGeom++ )
                {
                    tReference.push_back( box_reference_negative_volume( iCase.mPlanes( iGeom ) ) );
                    tReferenceError = std::max( tReferenceError, std::abs( tResult.mNegativeVolume( iGeom ) - tReference( iGeom ) ) );
                }

                if ( par_rank() != 0 )
                {
                    continue;
                }

                tNumThrew += tResult.mThrew;
                if ( not tResult.mThrew )
                {
                    tNumBadVolume += std::abs( tResult.mTotalVolume - gBoxVolume ) > gVolumeTolerance;
                    tNumReference += tReferenceError > gVolumeTolerance;
                    tNumCellDefects += tResult.mNumNonPositive + tResult.mNumNoPhase + tResult.mNumMixedSign > 0;
                    tNumDisagreements += tResult.mNumUnassignedIds + tResult.mNumCoordMismatch + tResult.mNumOwnerMismatch + tResult.mNumLabelMismatch > 0;
                }

                std::cout << std::setw( 12 ) << regime_name( iRegime ) << std::setw( 28 ) << iCase.mName;
                if ( tResult.mThrew )
                {
                    std::cout << "  THREW: " << tResult.mMessage.substr( 0, 120 ) << "\n";
                }
                else
                {
                    std::cout << std::setprecision( 3 ) << std::setw( 7 ) << tResult.mNumCells << std::setw( 6 ) << tResult.mNumSubphases
                              << std::setw( 12 ) << tResult.mMinCellVolume << std::setw( 6 ) << tResult.mNumNonPositive
                              << std::setw( 7 ) << tResult.mNumInterfaceFacets << std::setw( 12 ) << tResult.mInterfaceArea
                              << std::setw( 12 ) << tResult.mTotalVolume - gBoxVolume << std::setw( 12 ) << tReferenceError
                              << std::setw( 7 ) << tResult.mNumVertexIds << std::setw( 7 ) << tResult.mNumSharedIds
                              << std::setw( 7 ) << tResult.mNumCoordMismatch << std::setw( 7 ) << tResult.mNumOwnerMismatch
                              << std::setw( 7 ) << tResult.mNumLabelMismatch << std::setw( 7 ) << tResult.mNumDuplicateLocations
                              << std::setprecision( 6 ) << "\n";
                }

                std::ofstream* tStream = parallel_csv_stream();
                if ( tStream != nullptr )
                {
                    std::string tName = iCase.mName;
                    std::replace( tName.begin(), tName.end(), ',', ' ' );
                    *tStream << std::setprecision( 17 ) << par_size() << "," << regime_name( iRegime ) << "," << tName << "," << tResult.mThrew << ","
                             << tResult.mNumCells << "," << tResult.mNumSubphases << "," << tResult.mTotalVolume << ","
                             << tResult.mNegativeVolume( 0 ) << "," << tReference( 0 ) << ","
                             << ( tReference.size() > 1 ? tResult.mNegativeVolume( 1 ) : 0.0 ) << "," << ( tReference.size() > 1 ? tReference( 1 ) : 0.0 ) << ","
                             << tResult.mNumNonPositive << "," << tResult.mNumNoPhase << "," << tResult.mNumMixedSign << ","
                             << tResult.mNumInterfaceFacets << "," << tResult.mInterfaceArea << "," << tResult.mMinCellVolume << ","
                             << tResult.mMinSubphaseVolume << "," << tResult.mNumVertexIds << "," << tResult.mNumSharedIds << ","
                             << tResult.mNumUnassignedIds << "," << tResult.mNumCoordMismatch << "," << tResult.mMaxCoordMismatch << ","
                             << tResult.mNumOwnerMismatch << "," << tResult.mNumLabelMismatch << "," << tResult.mNumDuplicateLocations << "\n"
                             << std::flush;
                }
            }
        }

        if ( par_rank() == 0 )
        {
            std::cout << "\n  threw " << tNumThrew << " | total volume != 2: " << tNumBadVolume << " | negative volume != exact: " << tNumReference
                      << " | cases with collapsed/no-phase/mixed-sign cells: " << tNumCellDefects
                      << " | cases where the ranks disagree on a vertex: " << tNumDisagreements << "\n";
            CHECK( tNumThrew == 0 );
            CHECK( tNumBadVolume == 0 );
            CHECK( tNumDisagreements == 0 );
        }
    }


    TEST_CASE( "XTK pathology sweep: aura copies, enrichment and output sets agree across ranks", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI_AURA]" )
    {
        Grid_Mesh                tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 } );
        mtk::Interpolation_Mesh* tMesh = tBox.mMesh;

        const char* tOnlyCase = std::getenv( "MORIS_XTK_PATHOLOGY_MPI_CASE" );

        uint tNumThrew = 0, tNumCopyMismatch = 0, tNumDuplicateIds = 0;
        if ( par_rank() == 0 )
        {
            std::cout << "\n=== aura copies, enrichment and output sets, " << par_size() << " rank(s) ===\n";
        }

        for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC } )
        {
            Vector< Parallel_Case > tCases = parallel_cases();
            for ( uint iCaseIndex = 0; iCaseIndex < tCases.size(); iCaseIndex++ )
            {
                const Parallel_Case& iCase = tCases( iCaseIndex );
                if ( ( tOnlyCase != nullptr and std::atoi( tOnlyCase ) != (int)iCaseIndex ) or ( iCase.mPlanes.size() > 1 and iRegime != Regime::STORED ) )
                {
                    continue;
                }

                Aura_Result tResult = aura_case( tMesh, iCase.mPlanes, iRegime );
                if ( par_rank() != 0 )
                {
                    continue;
                }

                tNumThrew += tResult.mThrew;
                uint tMismatch  = tResult.mBgCellMismatch + tResult.mIpVertexMismatch + tResult.mIpCellMismatch + tResult.mClusterVerticesWithoutBasis
                               + tResult.mEmptyOwnedVertices;
                uint tDuplicate = tResult.mDuplicateIgCellIds + tResult.mInvalidIgCellIds + tResult.mDuplicateSubphaseIds + tResult.mInvalidSubphaseIds;
                tNumCopyMismatch += tMismatch > 0;
                tNumDuplicateIds += tDuplicate > 0;

                std::cout << std::setw( 9 ) << regime_name( iRegime ) << std::setw( 28 ) << iCase.mName;
                if ( tResult.mThrew )
                {
                    std::cout << "  THREW: " << tResult.mMessage.substr( 0, 120 ) << "\n";
                    continue;
                }
                std::cout << "  bg " << tResult.mBgCellMismatch << "/" << tResult.mSharedBgCells
                          << "  ipv " << tResult.mIpVertexMismatch << "/" << tResult.mSharedIpVertices
                          << "  ipc " << tResult.mIpCellMismatch << "/" << tResult.mSharedIpCells
                          << "  dup ig " << tResult.mDuplicateIgCellIds << "+" << tResult.mInvalidIgCellIds << "/" << tResult.mOwnedIgCells
                          << "  dup sp " << tResult.mDuplicateSubphaseIds << "+" << tResult.mInvalidSubphaseIds << "/" << tResult.mOwnedSubphases
                          << "  basis " << tResult.mNumEnrichedBasis << "  ipv " << tResult.mNumIpVertices << "  ipc " << tResult.mNumIpCells
                          << "  cluster vertices w/o T " << tResult.mClusterVerticesWithoutBasis << "  empty T owned/aura " << tResult.mEmptyOwnedVertices
                          << "/" << tResult.mEmptyAuraVertices << "  ipc id diff " << tResult.mIpCellIdMismatch
                          << "  ipv with several ids " << tResult.mIpVerticesWithSeveralIds << "  corner T diff " << tResult.mCornerTMatrixMismatch << "/"
                          << tResult.mSharedCorners << "\n";
                if ( not tResult.mCornerTMatrixExample.empty() )
                {
                    std::cout << "      first corner T diff: " << tResult.mCornerTMatrixExample << "\n";
                }
                for ( const std::string& iExample : { tResult.mBgCellExample, tResult.mIpVertexExample, tResult.mIpCellExample } )
                {
                    if ( not iExample.empty() )
                    {
                        std::cout << "      first mismatch: " << iExample << "\n";
                    }
                }

                std::string tName = iCase.mName;
                std::replace( tName.begin(), tName.end(), ',', ' ' );
                if ( std::ofstream* tStream = aura_csv_stream( "MORIS_XTK_PATHOLOGY_AURA_CSV",
                             "num_ranks,regime,case,threw,shared_bg,bg_mismatch,shared_ipv,ipv_mismatch,shared_ipc,ipc_mismatch,"
                             "owned_ig,dup_ig,invalid_ig,owned_sp,dup_sp,invalid_sp,enriched_basis,ip_vertices,ip_cells,cluster_vertices_without_t,"
                             "empty_t_owned,empty_t_aura,ipc_id_mismatch,example" ) )
                {
                    std::string tExample = tResult.mBgCellExample + " " + tResult.mIpVertexExample + " " + tResult.mIpCellExample;
                    std::replace( tExample.begin(), tExample.end(), ',', ' ' );
                    *tStream << par_size() << "," << regime_name( iRegime ) << "," << tName << "," << tResult.mThrew << ","
                             << tResult.mSharedBgCells << "," << tResult.mBgCellMismatch << "," << tResult.mSharedIpVertices << ","
                             << tResult.mIpVertexMismatch << "," << tResult.mSharedIpCells << "," << tResult.mIpCellMismatch << ","
                             << tResult.mOwnedIgCells << "," << tResult.mDuplicateIgCellIds << "," << tResult.mInvalidIgCellIds << ","
                             << tResult.mOwnedSubphases << "," << tResult.mDuplicateSubphaseIds << "," << tResult.mInvalidSubphaseIds << ","
                             << tResult.mNumEnrichedBasis << "," << tResult.mNumIpVertices << "," << tResult.mNumIpCells << ","
                             << tResult.mClusterVerticesWithoutBasis << "," << tResult.mEmptyOwnedVertices << "," << tResult.mEmptyAuraVertices << ","
                             << tResult.mIpCellIdMismatch << "," << tExample << "\n"
                             << std::flush;
                }
                if ( std::ofstream* tStream = aura_csv_stream( "MORIS_XTK_PATHOLOGY_SETS_CSV", "num_ranks,regime,case,set,owned_clusters,measure" ) )
                {
                    for ( const std::pair< const std::string, std::pair< uint, real > >& iSet : tResult.mSets )
                    {
                        *tStream << par_size() << "," << regime_name( iRegime ) << "," << tName << "," << iSet.first << ","
                                 << iSet.second.first << "," << std::setprecision( 17 ) << iSet.second.second << std::setprecision( 6 ) << "\n";
                    }
                    *tStream << std::flush;
                }
            }
        }

        if ( par_rank() == 0 )
        {
            std::cout << "\n  threw " << tNumThrew << " | cases with disagreeing copies: " << tNumCopyMismatch
                      << " | cases with duplicate or invalid owned ids: " << tNumDuplicateIds << "\n";
            CHECK( tNumThrew == 0 );
            CHECK( tNumCopyMismatch == 0 );
            CHECK( tNumDuplicateIds == 0 );
        }
    }

    TEST_CASE( "XTK pathology: a merged enriched IP vertex carries one id on every rank", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI_IPV_IDS]" )
    {
        // Finding 27 (studies/xtk_holistic_review/findings/15_walkthrough_stops_2_3.md §8.7): enriched IP vertex ids are
        // assigned per unzipped copy (one copy per subphase) before merge_duplicate_interpolation_vertices, and each rank
        // keeps the id of its own first local copy. The owner and the other ranks can therefore keep different copies'
        // ids for the same vertex. Checked here: every vertex, identified by its location and enriched T-matrix, has a
        // single id across all ranks. Vertices without a T-matrix (aura) are skipped. Run on 2 ranks:
        //   mpirun -n 2 XTK-test.exe "[XTK_PATHOLOGY_MPI_IPV_IDS]"
        if ( par_size() < 2 )
        {
            WARN( "needs at least 2 ranks" );
            return;
        }

        // linear B-splines on the Lagrange mesh give one-basis T-matrices; quadratic ones (as in our decks) have a
        // 3-cell support that reaches past a rank's cells plus aura, and refinement adds truncated two-level T-matrices.
        // The cuts below lie in unrefined cells, so the same-level guard for intersected cells is not triggered. Refining
        // from x > 1.25 instead puts the x = 1 cut cells next to refined cells and crashes in the subphase neighbourhood
        // (MB-075, findings/18_mb075_transition_segfault_fix.md), so the refined zone starts one cell further out.
        struct Mesh_Config
        {
            std::string           mName;
            std::array< uint, 3 > mElements;
            real                  mRefineAboveX;
            std::string           mBsplineOrders;
        };
        Vector< Mesh_Config > tConfigs = {
            { "linear B-splines 4x2x2", { 4, 2, 2 }, MORIS_REAL_MAX, "1" },
            { "quadratic B-splines 8x4x4", { 8, 4, 4 }, MORIS_REAL_MAX, "2" },
            { "quadratic, refined x>1.5", { 8, 4, 4 }, 1.5, "2" },
        };

        for ( const Mesh_Config& iConfig : tConfigs )
        {
            Grid_Mesh                tBox( iConfig.mElements, { 2.0, 1.0, 1.0 }, iConfig.mRefineAboveX, iConfig.mBsplineOrders );
            mtk::Interpolation_Mesh* tMesh = tBox.mMesh;
            if ( par_rank() == 0 )
            {
                std::cout << "\n=== " << iConfig.mName << ", " << par_size() << " ranks ===\n";
            }

            // report on rank 0 and check: one id per vertex, one vertex per id, one T-matrix per cell corner
            auto tReport = []( const std::string& aName, const Aura_Result& aResult ) {
                if ( par_rank() != 0 )
                {
                    return;
                }
                std::cout << std::setw( 28 ) << aName << "  vertices " << aResult.mDistinctIpVertices << "  with several ids "
                          << aResult.mIpVerticesWithSeveralIds << " (of which on one rank " << aResult.mIpVerticesWithSeveralIdsOnOneRank << ")"
                          << "  same id, different vertex " << aResult.mIpVertexMismatch << "/" << aResult.mSharedIpVertices
                          << "  cell corners with different T-matrices " << aResult.mCornerTMatrixMismatch << "/" << aResult.mSharedCorners
                          << "  bg cells per rank " << aResult.mMinBgCellsOnRank << "-" << aResult.mMaxBgCellsOnRank << "\n";
                if ( aResult.mThrew )
                {
                    std::cout << "      THREW: " << aResult.mMessage.substr( 0, 160 ) << "\n";
                }
                if ( not aResult.mIpVertexIdExample.empty() )
                {
                    std::cout << "      first: " << aResult.mIpVertexIdExample << "\n";
                }
                if ( not aResult.mCornerTMatrixExample.empty() )
                {
                    std::cout << "      first corner: " << aResult.mCornerTMatrixExample << "\n";
                }
                CHECK( not aResult.mThrew );
                CHECK( aResult.mIpVerticesWithSeveralIds == 0 );
                CHECK( aResult.mIpVertexMismatch == 0 );
                CHECK( aResult.mCornerTMatrixMismatch == 0 );
            };

            for ( const Parallel_Case& iCase : parallel_cases() )
            {
                // uncut, cut away from the processor face, cut on it
                if ( iCase.mName != "x=3 (no cut, control)" and iCase.mName != "x=0.3 (away from the face)" and iCase.mName != "x=1 (processor face)" )
                {
                    continue;
                }
                tReport( iCase.mName, aura_case( tMesh, iCase.mPlanes, Regime::STORED ) );
            }

            // thin features need node rows at y = 0.25 and 0.75 (h = 0.25): the unrefined 8 x 4 x 4 mesh only (on the
            // refined one they reach x > 1.5, cut cells sit on two levels and decomposition is refused)
            if ( iConfig.mElements[ 0 ] != 8 or iConfig.mRefineAboveX < MORIS_REAL_MAX )
            {
                continue;
            }
            for ( const Thin_Case& iCase : thin_cases() )
            {
                gBoxes = iCase.mBoxes;
                tReport( iCase.mName, aura_case_fields( tMesh, { &boxes_field }, Regime::STORED ) );
            }
        }
    }

    TEST_CASE( "XTK pathology: interfaces at the outer edge of the aura", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI_AURA_EDGE]" )
    {
        // An interface near the edge of the aura (x = 0.5 and x = 1.5 on 2 ranks, see aura_edge_cases) rather than on
        // the processor face: the cut cell in the outermost aura layer lacks its outer neighbours on that rank, and its
        // outer vertices have no T-matrix there. Every copy of a cell must still be decomposed, enriched and numbered
        // as on its owner. Run on 2 ranks (4 ranks splits the box differently, the edges then lie elsewhere):
        //   mpirun -n 2 XTK-test.exe "[XTK_PATHOLOGY_MPI_AURA_EDGE]"
        if ( par_size() < 2 )
        {
            WARN( "needs at least 2 ranks" );
            return;
        }

        struct Mesh_Config
        {
            std::string           mName;
            std::array< uint, 3 > mElements;
            real                  mRefineAboveX;
            std::string           mBsplineOrders;
            Vector< real >        mEdges;    // aura edges to cut at; the refined mesh only away from its refined zone
        };
        Vector< Mesh_Config > tConfigs = {
            { "linear B-splines 4x2x2", { 4, 2, 2 }, MORIS_REAL_MAX, "1", { 0.5, 1.5 } },
            { "quadratic B-splines 8x4x4", { 8, 4, 4 }, MORIS_REAL_MAX, "2", { 0.5, 1.5 } },
            { "quadratic, refined x>1.5", { 8, 4, 4 }, 1.5, "2", { 0.5 } },
        };

        uint tNumFailedCases = 0;
        for ( const Mesh_Config& iConfig : tConfigs )
        {
            Grid_Mesh                tBox( iConfig.mElements, { 2.0, 1.0, 1.0 }, iConfig.mRefineAboveX, iConfig.mBsplineOrders );
            mtk::Interpolation_Mesh* tMesh = tBox.mMesh;
            if ( par_rank() == 0 )
            {
                std::cout << "\n=== " << iConfig.mName << ", " << par_size() << " ranks ===\n";
            }

            for ( const Parallel_Case& iCase : aura_edge_cases( iConfig.mEdges ) )
            {
                Aura_Result tResult = aura_case( tMesh, iCase.mPlanes, Regime::STORED );
                if ( par_rank() != 0 )
                {
                    continue;
                }

                // copies that disagree with the owner, or data missing on a rank that assembles the cluster
                uint tCopies = tResult.mBgCellMismatch + tResult.mIpVertexMismatch + tResult.mIpCellMismatch;
                uint tMissing = tResult.mClusterVerticesWithoutBasis + tResult.mEmptyOwnedVertices;
                uint tIds     = tResult.mIpVerticesWithSeveralIds + tResult.mDuplicateIgCellIds + tResult.mInvalidIgCellIds
                          + tResult.mDuplicateSubphaseIds + tResult.mInvalidSubphaseIds;
                bool tFailed = tResult.mThrew or tCopies > 0 or tMissing > 0 or tIds > 0 or tResult.mCornerTMatrixMismatch > 0;
                tNumFailedCases += tFailed;

                std::cout << ( tFailed ? "  X " : "    " ) << std::setw( 22 ) << iCase.mName;
                if ( tResult.mThrew )
                {
                    std::cout << "  THREW: " << tResult.mMessage.substr( 0, 160 ) << "\n";
                    continue;
                }
                std::cout << "  bg " << tResult.mBgCellMismatch << "/" << tResult.mSharedBgCells
                          << "  ipc " << tResult.mIpCellMismatch << "/" << tResult.mSharedIpCells
                          << "  ipv " << tResult.mIpVertexMismatch << "/" << tResult.mSharedIpVertices
                          << "  T corners " << tResult.mCornerTMatrixMismatch << "/" << tResult.mSharedCorners
                          << "  ipv several ids " << tResult.mIpVerticesWithSeveralIds
                          << "  cluster vertices w/o T " << tResult.mClusterVerticesWithoutBasis
                          << "  empty T owned/aura " << tResult.mEmptyOwnedVertices << "/" << tResult.mEmptyAuraVertices
                          << "  dup ig/sp " << tResult.mDuplicateIgCellIds + tResult.mInvalidIgCellIds << "/"
                          << tResult.mDuplicateSubphaseIds + tResult.mInvalidSubphaseIds << "\n";
                for ( const std::string& iExample : { tResult.mBgCellExample, tResult.mIpCellExample, tResult.mIpVertexExample,
                              tResult.mCornerTMatrixExample, tResult.mIpVertexIdExample } )
                {
                    if ( not iExample.empty() )
                    {
                        std::cout << "        first: " << iExample << "\n";
                    }
                }
            }
        }

        if ( par_rank() == 0 )
        {
            std::cout << "\n  cases with a disagreement, missing data or a throw: " << tNumFailedCases << "\n";
            CHECK( tNumFailedCases == 0 );
        }
    }

    TEST_CASE( "XTK pathology: intersected cells on different refinement levels on different ranks", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_MPI_LEVELS]" )
    {
        // finding 6: check_intersected_background_cell_levels compares each rank's intersected cells with that rank's
        // first one and combines with all_land, so ranks that are each uniform at different levels pass. The box is
        // refined one level for x > 1.5; the processor face is x = 1 on 2 ranks, and each rank's aura reaches 0.5 past it.
        Grid_Mesh                tBox( { 4, 2, 2 }, { 2.0, 1.0, 1.0 }, 1.5 );
        mtk::Interpolation_Mesh* tMesh = tBox.mMesh;

        struct Level_Case
        {
            std::string     mName;
            Vector< Plane > mPlanes;
        };
        const Point          tXAxis = { 1.0, 0.0, 0.0 };
        Vector< Level_Case > tCases = {
            { "x=1.8 only (fine cells)", { plane_through( tXAxis, { 1.8, 0.5, 0.5 }, 0.0 ) } },
            { "x=0.3 only (coarse cells)", { plane_through( tXAxis, { 0.3, 0.5, 0.5 }, 0.0 ) } },
            { "x=0.3 and x=1.8", { plane_through( tXAxis, { 0.3, 0.5, 0.5 }, 0.0 ), plane_through( tXAxis, { 1.8, 0.5, 0.5 }, 0.0 ) } },
            { "slanted across x=1.5", { plane_through( { 1.0, 0.3, 0.2 }, { 1.5, 0.5, 0.5 }, 0.0 ) } }
        };

        uint tNumInconsistent = 0;
        for ( const Level_Case& iCase : tCases )
        {
            // levels of the cells each plane cuts on this rank (corner signs), owned and aura
            std::set< uint > tLevels;
            for ( uint iCell = 0; iCell < tMesh->get_num_elems(); iCell++ )
            {
                mtk::Cell& tCell = tMesh->get_mtk_cell( iCell );
                for ( const Plane& iPlane : iCase.mPlanes )
                {
                    bool tNegative = false, tPositive = false;
                    for ( mtk::Vertex* iVertex : tCell.get_vertex_pointers() )
                    {
                        Matrix< DDRMat > tX   = iVertex->get_coords();
                        real             tPhi = iPlane[ 0 ] * tX( 0 ) + iPlane[ 1 ] * tX( 1 ) + iPlane[ 2 ] * tX( 2 ) + iPlane[ 3 ];
                        tNegative             = tNegative or tPhi < 0.0;
                        tPositive             = tPositive or tPhi > 0.0;
                    }
                    if ( tNegative and tPositive )
                    {
                        tLevels.insert( tCell.get_level() );
                    }
                }
            }

            // decompose: true if the level check passed on all ranks
            bool tDecomposed = false;
            try
            {
                Vector< std::shared_ptr< gen::Geometry > > tGeometries;
                for ( uint iPlane = 0; iPlane < iCase.mPlanes.size(); iPlane++ )
                {
                    gPlanes[ iPlane ] = iCase.mPlanes( iPlane );
                    tGeometries.push_back( create_geometry( gPlaneFields[ iPlane ], Regime::STORED ) );
                }
                gen::Geometry_Engine_Parameters tParameters;
                tParameters.mGeometries = tGeometries;
                gen::Geometry_Engine tGeometryEngine( tMesh, tParameters );
                Model                tModel( 3, tMesh, &tGeometryEngine );
                tModel.mVerbose = false;
                tDecomposed     = tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );
            }
            catch ( const std::exception& aException )
            {
                std::cout << "rank " << par_rank() << " " << iCase.mName << " THREW: " << std::string( aException.what() ).substr( 0, 160 ) << "\n";
            }

            // the levels every rank sees, and whether the whole problem has mixed levels
            uint tMinLevel = min_all( tLevels.empty() ? (uint)99 : *tLevels.begin() );
            uint tMaxLevel = max_all( tLevels.empty() ? (uint)0 : *tLevels.rbegin() );
            bool tMixed    = tMinLevel < tMaxLevel;

            std::ostringstream tLine;
            tLine << "rank " << par_rank() << " of " << par_size() << "  " << std::setw( 26 ) << iCase.mName << "  cut levels on this rank {";
            for ( uint iLevel : tLevels )
            {
                tLine << " " << iLevel;
            }
            tLine << " }  globally " << ( tMixed ? "MIXED" : "uniform" ) << "  decomposed " << tDecomposed << "\n";
            std::cout << tLine.str() << std::flush;

            // the level check must refuse every mixed-level problem, on any number of ranks
            if ( par_rank() == 0 and tMixed == tDecomposed )
            {
                tNumInconsistent++;
                std::cout << "  ^ " << ( tDecomposed ? "decomposed although the cut cells are on different levels" : "refused although every cut cell is on one level" ) << "\n";
            }
        }
        if ( par_rank() == 0 )
        {
            CHECK( tNumInconsistent == 0 );
        }
    }


}    // namespace moris::xtk
