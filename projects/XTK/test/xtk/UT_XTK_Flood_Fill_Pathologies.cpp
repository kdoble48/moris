/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Flood_Fill_Pathologies.cpp
 *
 * The generic flood fill (fn_mesh_flood_fill.hpp) against union-find on random and adversarial graphs, its
 * active-front bound, and the SortPoint ordering of enrichment-level centroids. Split out of
 * UT_XTK_Connectivity_Pathologies.cpp; the helpers are in cl_XTK_Connectivity_Test_Harness.hpp.
 *
 */

#include "catch.hpp"

#include "cl_XTK_Connectivity_Test_Harness.hpp"

namespace moris::xtk
{
    using namespace connectivity_test;

    TEST_CASE( "XTK connectivity: flood fill against union-find on random and adversarial graphs", "[XTK][XTK_CONNECTIVITY_FLOOD_FILL]" )
    {
        if ( par_rank() != 0 )
        {
            return;
        }
        std::mt19937 tRandom( 20261021 );
        Fuzz_Outcome tOutcome;

        // random graphs: 1 to 60 nodes, mean degree 0 to 8, 1 to 3 phases; the neighbour matrix is either exactly as wide as
        // the largest row (full rows have no dummy end marker, the branch at fn_mesh_flood_fill.hpp:231) or 50 wide
        for ( uint iTrial = 0; iTrial < 3000; iTrial++ )
        {
            uint                  tNumNodes  = 1 + tRandom() % 60;
            real                  tDegree    = 8.0 * ( tRandom() % 1000 ) / 1000.0;
            Graph                 tGraph     = random_graph( tNumNodes, tDegree, tRandom );
            Vector< moris_index > tPhases    = random_phases( tNumNodes, 1 + tRandom() % 3, tRandom );
            uint                  tNumColumns = iTrial % 2 ? std::max( 1u, max_degree( tGraph ) ) : 50;
            if ( max_degree( tGraph ) > tNumColumns )
            {
                continue;
            }
            fuzz_one( tGraph, tPhases, tNumColumns, tOutcome );
        }

        // adversarial graphs, one phase: a single node, isolated nodes, a long path (every inner row full), a star, a
        // ladder of triangles, and duplicate neighbour entries (a neighbour listed twice in a row)
        fuzz_one( Graph( 1 ), Vector< moris_index >( 1, 1 ), 1, tOutcome );
        fuzz_one( Graph( 25 ), Vector< moris_index >( 25, 1 ), 3, tOutcome );
        Graph tPath( 5000 );
        for ( uint iNode = 0; iNode + 1 < 5000; iNode++ )
        {
            add_edge( tPath, iNode, iNode + 1 );
        }
        fuzz_one( tPath, Vector< moris_index >( 5000, 1 ), 2, tOutcome );
        fuzz_one( tPath, Vector< moris_index >( 5000, 1 ), 50, tOutcome );
        Graph tStar( 400 );
        for ( uint iNode = 1; iNode < 400; iNode++ )
        {
            add_edge( tStar, 0, iNode );
        }
        fuzz_one( tStar, Vector< moris_index >( 400, 1 ), 399, tOutcome );
        Graph tLadder( 300 );
        for ( uint iNode = 0; iNode + 2 < 300; iNode++ )
        {
            add_edge( tLadder, iNode, iNode + 1 );
            add_edge( tLadder, iNode, iNode + 2 );
        }
        fuzz_one( tLadder, Vector< moris_index >( 300, 1 ), 4, tOutcome );
        Graph tDuplicates = random_graph( 40, 3.0, tRandom );
        for ( std::vector< moris_index >& iRow : tDuplicates )
        {
            std::vector< moris_index > tCopy = iRow;
            iRow.insert( iRow.end(), tCopy.begin(), tCopy.end() );
        }
        fuzz_one( tDuplicates, Vector< moris_index >( 40, 1 ), std::max( 1u, max_degree( tDuplicates ) ), tOutcome );

        // complete graphs below the measured bound (K_n with n <= 20 stays inside 10 n + 1, see [XTK_CONNECTIVITY_FLOOD_FILL_BOUND])
        for ( uint tN = 2; tN <= 20; tN++ )
        {
            fuzz_one( complete_graph( tN ), Vector< moris_index >( tN, 1 ), tN - 1, tOutcome );
        }

        std::cout << "flood fill fuzz: " << tOutcome.mGraphs << " graphs, " << tOutcome.mRun << " through flood_fill, " << tOutcome.mSkipped
                  << " skipped (front would overflow), replica differs " << tOutcome.mReplicaWrong << ", partition wrong " << tOutcome.mWrong
                  << ", max value wrong " << tOutcome.mBadMax << "\n";
        CHECK( tOutcome.mReplicaWrong == 0 );
        CHECK( tOutcome.mWrong == 0 );
        CHECK( tOutcome.mBadMax == 0 );
        CHECK( tOutcome.mRun > 2500 );

        // sort_points_by_coordinates on well separated points sorts lexicographically
        uint tSortWrong = 0;
        for ( uint iTrial = 0; iTrial < 500; iTrial++ )
        {
            uint                         tNumPoints = 1 + tRandom() % 12;
            Vector< Matrix< DDRMat > >   tPoints;
            std::vector< std::array< real, 3 > > tPlain;
            for ( uint iPoint = 0; iPoint < tNumPoints; iPoint++ )
            {
                std::array< real, 3 > tX = { real( tRandom() % 7 ), real( tRandom() % 7 ), real( iPoint ) };
                tPlain.push_back( tX );
                tPoints.push_back( Matrix< DDRMat >( { { tX[ 0 ] }, { tX[ 1 ] }, { tX[ 2 ] } } ) );
            }
            Vector< moris_index > tOrder = sort_points_by_coordinates( tPoints );
            for ( uint iPoint = 1; iPoint < tNumPoints; iPoint++ )
            {
                tSortWrong += not( tPlain[ tOrder( iPoint - 1 ) ] < tPlain[ tOrder( iPoint ) ] );
            }
        }
        CHECK( tSortWrong == 0 );
    }

    TEST_CASE( "XTK connectivity: active-front bound of the generic flood fill (10N+1 front)", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_FLOOD_FILL_BOUND]" )
    {
        // flood_fill sizes its active front 10 N + 1 (fn_mesh_flood_fill.hpp:80) and pushes a neighbour whenever it is not
        // yet labelled, so a node can sit on the front many times. The overflow check is a MORIS_ASSERT (:222), compiled
        // out unless MORIS_HAVE_DEBUG: an optimized build writes past the allocation. This measures where that happens.
        if ( par_rank() != 0 )
        {
            return;
        }
        std::ofstream tCsv;
        if ( const char* tPath = std::getenv( "MORIS_XTK_CONNECTIVITY_BOUND_CSV" ) )
        {
            tCsv.open( tPath );
            tCsv << "graph,nodes,mean_degree,samples,overflow_samples,max_peak_front,front_size,max_ratio\n";
        }

        std::cout << "\n=== complete graphs K_n: peak active front vs 10 n + 1 ===\n"
                  << std::setw( 6 ) << "n" << std::setw( 10 ) << "peak" << std::setw( 10 ) << "10n+1" << std::setw( 10 ) << "ratio" << "\n";
        uint tFirstOverflow = 0;
        for ( uint tN = 2; tN <= 64; tN++ )
        {
            Graph          tGraph   = complete_graph( tN );
            Replica_Result tReplica = flood_fill_replica( graph_matrix( tGraph, tN - 1 ), Vector< moris_index >( tN, 1 ) );
            if ( tReplica.mOverflow and tFirstOverflow == 0 )
            {
                tFirstOverflow = tN;
            }
            real tRatio = real( tReplica.mPeakFront ) / tReplica.mFrontSize;
            if ( tN <= 30 or tN % 8 == 0 )
            {
                std::cout << std::setw( 6 ) << tN << std::setw( 10 ) << tReplica.mPeakFront << std::setw( 10 ) << tReplica.mFrontSize
                          << std::setw( 10 ) << std::setprecision( 3 ) << tRatio << std::setprecision( 6 ) << ( tReplica.mOverflow ? "  OVERFLOW" : "" ) << "\n";
            }
            if ( tCsv.is_open() )
            {
                tCsv << "complete," << tN << "," << tN - 1 << ",1," << tReplica.mOverflow << "," << tReplica.mPeakFront << "," << tReplica.mFrontSize << ","
                     << tRatio << "\n";
            }
        }
        std::cout << "first complete graph that overflows: K_" << tFirstOverflow << "\n";

        // random graphs: overflow frequency against node count and mean degree
        std::mt19937 tRandom( 7 );
        std::cout << "\n=== random graphs: samples (of 200) whose front overflows, and the largest peak / (10 N + 1) ===\n"
                  << std::setw( 8 ) << "N" << std::setw( 8 ) << "degree" << std::setw( 10 ) << "overflow" << std::setw( 12 ) << "max ratio" << "\n";
        for ( uint tN : { 10u, 20u, 50u, 100u, 400u } )
        {
            for ( real tDegree : { 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 16.0, 20.0, 30.0, 40.0 } )
            {
                if ( tDegree > tN - 1 )
                {
                    continue;
                }
                uint tOverflow = 0;
                uint tMaxPeak  = 0;
                real tMaxRatio = 0.0;
                for ( uint iSample = 0; iSample < 200; iSample++ )
                {
                    Graph          tGraph   = random_graph( tN, tDegree, tRandom );
                    Replica_Result tReplica = flood_fill_replica( graph_matrix( tGraph, std::max( 1u, max_degree( tGraph ) ) ), Vector< moris_index >( tN, 1 ) );
                    tOverflow += tReplica.mOverflow;
                    tMaxPeak  = std::max( tMaxPeak, tReplica.mPeakFront );
                    tMaxRatio = std::max( tMaxRatio, real( tReplica.mPeakFront ) / tReplica.mFrontSize );
                }
                std::cout << std::setw( 8 ) << tN << std::setw( 8 ) << tDegree << std::setw( 10 ) << tOverflow << std::setw( 12 ) << std::setprecision( 3 )
                          << tMaxRatio << std::setprecision( 6 ) << "\n";
                if ( tCsv.is_open() )
                {
                    tCsv << "random," << tN << "," << tDegree << ",200," << tOverflow << "," << tMaxPeak << "," << 10 * tN + 1 << "," << tMaxRatio << "\n";
                }
            }
        }

        // MORIS_XTK_CONNECTIVITY_OVERFLOW_DEMO=<n>: call the real flood_fill on K_n. Above the bound this writes past the
        // heap allocation of the front; run it in its own process (e.g. with MALLOC_CHECK_=3) to see what an optimized
        // build does.
        if ( const char* tDemo = std::getenv( "MORIS_XTK_CONNECTIVITY_OVERFLOW_DEMO" ) )
        {
            uint                  tN     = std::atoi( tDemo );
            Graph                 tGraph = complete_graph( tN );
            moris_index           tMax   = -1;
            Vector< moris_index > tLabels = real_flood_fill( graph_matrix( tGraph, tN - 1 ), Vector< moris_index >( tN, 1 ), tMax );
            std::cout << "overflow demo: real flood_fill on K_" << tN << " returned, max value " << tMax << ", components "
                      << count_distinct( tLabels ) << " (expected 1)\n";
        }
        CHECK( tFirstOverflow > 0 );
    }

    TEST_CASE( "XTK connectivity: SortPoint ordering of enrichment-level centroids", "[.][XTK][XTK_CONNECTIVITY][XTK_CONNECTIVITY_SORTPOINT]" )
    {
        // Enrichment::sort_enrichment_levels_in_basis_support (sort_basis_enrichment_levels, default off) sorts levels by
        // centroid with SortPoint::operator< (fn_sort_points_by_coordinates.hpp:34). Two defects are demonstrated:
        //   (a) coincident centroids stop the run (MORIS_ERROR at :48),
        //   (b) the epsilon tolerance makes the comparison intransitive, which std::stable_sort does not allow.
        if ( par_rank() != 0 )
        {
            return;
        }
        const real tEps = MORIS_REAL_EPS;

        // (a) two symmetric regions, e.g. a frame around a central island, share their centroid
        Vector< Matrix< DDRMat > > tCoincident = { Matrix< DDRMat >( { { 0.5 }, { 0.5 }, { 0.5 } } ), Matrix< DDRMat >( { { 0.5 }, { 0.5 }, { 0.5 } } ) };
        bool        tThrew = false;
        std::string tMessage;
        try
        {
            sort_points_by_coordinates( tCoincident );
        }
        catch ( const std::exception& aException )
        {
            tThrew   = true;
            tMessage = aException.what();
        }
        std::cout << "(a) coincident centroids: threw " << tThrew << ": " << tMessage.substr( 0, 160 ) << "\n";
        CHECK( tThrew );

        // (b) P1 = (0, 1), P2 = (0.6 eps, 0), P3 = (1.2 eps, -1): x differs by <= eps for P1,P2 and P2,P3 but not P1,P3
        SortPoint tP1( Matrix< DDRMat >( { { 0.0 }, { 1.0 } } ) );
        SortPoint tP2( Matrix< DDRMat >( { { 0.6 * tEps }, { 0.0 } } ) );
        SortPoint tP3( Matrix< DDRMat >( { { 1.2 * tEps }, { -1.0 } } ) );
        bool      tP2BeforeP1 = tP2 < tP1;
        bool      tP3BeforeP2 = tP3 < tP2;
        bool      tP1BeforeP3 = tP1 < tP3;
        std::cout << "(b) P2<P1 " << tP2BeforeP1 << ", P3<P2 " << tP3BeforeP2 << ", P1<P3 " << tP1BeforeP3
                  << ( tP2BeforeP1 and tP3BeforeP2 and tP1BeforeP3 ? "  -> cycle: not a strict weak ordering" : "" ) << "\n";
        CHECK( ( tP2BeforeP1 and tP3BeforeP2 and tP1BeforeP3 ) );

        // the sorted order of the three points depends on their input order
        std::set< std::vector< uint > > tOutcomes;
        std::vector< uint >             tPermutation = { 0, 1, 2 };
        Vector< Matrix< DDRMat > >      tBase         = { tP1.mX, tP2.mX, tP3.mX };
        do
        {
            Vector< Matrix< DDRMat > > tPoints;
            for ( uint iPoint : tPermutation )
            {
                tPoints.push_back( tBase( iPoint ) );
            }
            Vector< moris_index > tOrder = sort_points_by_coordinates( tPoints );
            std::vector< uint >   tSorted;
            for ( moris_index iIndex : tOrder )
            {
                tSorted.push_back( tPermutation[ iIndex ] );
            }
            tOutcomes.insert( tSorted );
        } while ( std::next_permutation( tPermutation.begin(), tPermutation.end() ) );
        std::cout << "    distinct sorted orders over the 6 input permutations: " << tOutcomes.size() << " (1 for a valid ordering)\n";
        CHECK( tOutcomes.size() > 1 );
    }

}    // namespace moris::xtk
