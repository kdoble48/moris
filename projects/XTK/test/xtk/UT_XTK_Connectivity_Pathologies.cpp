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

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include "cl_Communication_Tools.hpp"

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "fn_mesh_flood_fill.hpp"
#include "fn_sort_points_by_coordinates.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Database.hpp"
#include "cl_HMR_Background_Mesh_Base.hpp"
#include "cl_HMR_Background_Element_Base.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_User_Defined_Field.hpp"

namespace moris::xtk
{
    namespace
    {
        using Point = std::array< real, 3 >;

        // ----------------------------------------------------------------------------------
        // union-find and partition comparison

        struct Union_Find
        {
            Vector< uint > mParent;

            explicit Union_Find( uint aSize )
                    : mParent( aSize )
            {
                std::iota( mParent.begin(), mParent.end(), 0 );
            }

            uint
            find( uint aElement )
            {
                while ( mParent( aElement ) != aElement )
                {
                    mParent( aElement ) = mParent( mParent( aElement ) );
                    aElement            = mParent( aElement );
                }
                return aElement;
            }

            void
            unite( uint aFirst, uint aSecond )
            {
                uint tFirst  = this->find( aFirst );
                uint tSecond = this->find( aSecond );
                if ( tFirst != tSecond )
                {
                    mParent( std::max( tFirst, tSecond ) ) = std::min( tFirst, tSecond );
                }
            }

            // component label per element (the smallest element of its component)
            Vector< moris_index >
            labels()
            {
                Vector< moris_index > tLabels( mParent.size() );
                for ( uint iElement = 0; iElement < mParent.size(); iElement++ )
                {
                    tLabels( iElement ) = this->find( iElement );
                }
                return tLabels;
            }
        };

        uint
        count_distinct( const Vector< moris_index >& aLabels )
        {
            std::set< moris_index > tDistinct( aLabels.begin(), aLabels.end() );
            return tDistinct.size();
        }

        // two labellings of the same elements: mMerged counts test labels that span several reference labels (cross-talk),
        // mSplit counts reference labels that span several test labels
        struct Partition_Comparison
        {
            uint mMerged = 0;
            uint mSplit  = 0;

            bool
            equal() const
            {
                return mMerged == 0 and mSplit == 0;
            }
        };

        Partition_Comparison
        compare_partitions( const Vector< moris_index >& aReference, const Vector< moris_index >& aTest )
        {
            std::map< moris_index, std::set< moris_index > > tReferenceOfTest;
            std::map< moris_index, std::set< moris_index > > tTestOfReference;
            for ( uint iElement = 0; iElement < aReference.size(); iElement++ )
            {
                tReferenceOfTest[ aTest( iElement ) ].insert( aReference( iElement ) );
                tTestOfReference[ aReference( iElement ) ].insert( aTest( iElement ) );
            }
            Partition_Comparison tComparison;
            for ( const std::pair< const moris_index, std::set< moris_index > >& iTest : tReferenceOfTest )
            {
                tComparison.mMerged += iTest.second.size() > 1;
            }
            for ( const std::pair< const moris_index, std::set< moris_index > >& iReference : tTestOfReference )
            {
                tComparison.mSplit += iReference.second.size() > 1;
            }
            return tComparison;
        }

        // ----------------------------------------------------------------------------------
        // section 1: the generic flood fill (fn_mesh_flood_fill.hpp) on arbitrary graphs

        using Graph = std::vector< std::vector< moris_index > >;

        // neighbour matrix as the enrichment builds it: one row per node, unused columns hold the dummy MORIS_INDEX_MAX
        Matrix< IndexMat >
        graph_matrix( const Graph& aGraph, uint aNumColumns )
        {
            Matrix< IndexMat > tMatrix( aGraph.size(), aNumColumns, MORIS_INDEX_MAX );
            for ( uint iNode = 0; iNode < aGraph.size(); iNode++ )
            {
                MORIS_ERROR( aGraph[ iNode ].size() <= aNumColumns, "graph_matrix() - row longer than the number of columns" );
                for ( uint iNeighbor = 0; iNeighbor < aGraph[ iNode ].size(); iNeighbor++ )
                {
                    tMatrix( iNode, iNeighbor ) = aGraph[ iNode ][ iNeighbor ];
                }
            }
            return tMatrix;
        }

        // connected components of the same-phase subgraph
        Vector< moris_index >
        reference_flood_fill( const Graph& aGraph, const Vector< moris_index >& aPhases )
        {
            Union_Find tUnionFind( aGraph.size() );
            for ( uint iNode = 0; iNode < aGraph.size(); iNode++ )
            {
                for ( moris_index iNeighbor : aGraph[ iNode ] )
                {
                    if ( aPhases( iNode ) == aPhases( iNeighbor ) )
                    {
                        tUnionFind.unite( iNode, iNeighbor );
                    }
                }
            }
            return tUnionFind.labels();
        }

        struct Replica_Result
        {
            Vector< moris_index > mLabels;
            uint                  mFrontSize = 0;        // the allocation of flood_fill: 10 N + 1
            uint                  mPeakFront = 0;        // largest position written + 1
            bool                  mOverflow  = false;    // a write at or past mFrontSize
        };

        void
        write_front( std::vector< moris_index >& aFront, uint aPosition, moris_index aValue, Replica_Result& aResult )
        {
            if ( aPosition >= aFront.size() )
            {
                aFront.resize( 2 * aPosition + 2, MORIS_INDEX_MAX );
            }
            aFront[ aPosition ] = aValue;
            aResult.mPeakFront  = std::max( aResult.mPeakFront, aPosition + 1 );
            aResult.mOverflow   = aResult.mOverflow or aPosition >= aResult.mFrontSize;
        }

        // a statement-by-statement copy of flood_fill (fn_mesh_flood_fill.hpp:38-257, aIncludeAllElements = true, every
        // element included) whose active front grows instead of overflowing, so the peak front can be measured
        Replica_Result
        flood_fill_replica( const Matrix< IndexMat >& aGraph, const Vector< moris_index >& aPhases )
        {
            const moris_index tDummy        = MORIS_INDEX_MAX;
            uint              tNumElements  = aGraph.n_rows();
            uint              tMaxNeighbors = aGraph.n_cols();

            Replica_Result tResult;
            tResult.mFrontSize = 10 * tNumElements + 1;
            tResult.mLabels    = Vector< moris_index >( tNumElements, tDummy );

            std::vector< bool >        tSet( tNumElements, false );
            std::vector< moris_index > tFront( tResult.mFrontSize, 0 );
            uint                       tCount   = 0;
            moris_index                tCurrent = 0;

            for ( uint iE = 0; iE < tNumElements; iE++ )
            {
                if ( tSet[ iE ] )
                {
                    continue;
                }
                moris_index tPhase    = aPhases( iE );
                tResult.mLabels( iE ) = tCurrent;
                tSet[ iE ]            = true;

                for ( uint iN = 0; iN < tMaxNeighbors; iN++ )
                {
                    moris_index tNeighbor = aGraph( iE, iN );
                    if ( tNeighbor == tDummy )
                    {
                        break;
                    }
                    if ( not tSet[ tNeighbor ] and aPhases( tNeighbor ) == tPhase )
                    {
                        write_front( tFront, tCount, tNeighbor, tResult );
                        tCount++;
                    }
                }

                while ( tCount != 0 )
                {
                    moris_index tTop = tFront[ tCount - 1 ];
                    if ( aPhases( tTop ) == tPhase and not tSet[ tTop ] )
                    {
                        tResult.mLabels( tTop ) = tCurrent;
                        tSet[ tTop ]            = true;

                        bool tReplaced = false;
                        for ( uint iN = 0; iN < tMaxNeighbors; iN++ )
                        {
                            moris_index tNeighbor = aGraph( tTop, iN );
                            if ( tNeighbor == tDummy )
                            {
                                if ( not tReplaced )
                                {
                                    tCount--;
                                }
                                break;
                            }
                            if ( not tSet[ tNeighbor ] )
                            {
                                if ( not tReplaced )
                                {
                                    tFront[ tCount - 1 ] = tNeighbor;
                                    tReplaced            = true;
                                }
                                else
                                {
                                    write_front( tFront, tCount, tNeighbor, tResult );
                                    tCount++;
                                }
                            }
                            if ( iN == tMaxNeighbors - 1 and not tReplaced )
                            {
                                tFront[ tCount - 1 ] = tDummy;
                                tCount--;
                            }
                        }
                    }
                    else
                    {
                        tFront[ tCount - 1 ] = tDummy;
                        tCount--;
                    }
                }
                tCurrent++;
            }
            return tResult;
        }

        // the real flood_fill, called as Enrichment::assign_subphase_bin_enrichment_levels_in_basis_support does
        Vector< moris_index >
        real_flood_fill( const Matrix< IndexMat >& aGraph, const Vector< moris_index >& aPhases, moris_index& aMaxValueAssigned )
        {
            uint               tNumElements = aGraph.n_rows();
            Matrix< IndexMat > tActive( 1, tNumElements );
            Matrix< IndexMat > tPhases( 1, tNumElements );
            for ( uint iE = 0; iE < tNumElements; iE++ )
            {
                tActive( 0, iE ) = iE;
                tPhases( 0, iE ) = aPhases( iE );
            }
            Matrix< IndexMat > tIncluded( 1, tNumElements, 1 );
            Matrix< IndexMat > tLabels = flood_fill( aGraph, tPhases, tActive, tIncluded, 3, MORIS_INDEX_MAX, aMaxValueAssigned, true );

            Vector< moris_index > tResult( tNumElements );
            for ( uint iE = 0; iE < tNumElements; iE++ )
            {
                tResult( iE ) = tLabels( iE );
            }
            return tResult;
        }

        void
        add_edge( Graph& aGraph, uint aFirst, uint aSecond )
        {
            aGraph[ aFirst ].push_back( aSecond );
            aGraph[ aSecond ].push_back( aFirst );
        }

        Graph
        complete_graph( uint aNumNodes )
        {
            Graph tGraph( aNumNodes );
            for ( uint iFirst = 0; iFirst < aNumNodes; iFirst++ )
            {
                for ( uint iSecond = iFirst + 1; iSecond < aNumNodes; iSecond++ )
                {
                    add_edge( tGraph, iFirst, iSecond );
                }
            }
            return tGraph;
        }

        // Erdos-Renyi graph with the given mean degree
        Graph
        random_graph( uint aNumNodes, real aMeanDegree, std::mt19937& aRandom )
        {
            Graph                                  tGraph( aNumNodes );
            real                                   tProbability = aNumNodes > 1 ? std::min( 1.0, aMeanDegree / ( aNumNodes - 1 ) ) : 0.0;
            std::uniform_real_distribution< real > tUniform( 0.0, 1.0 );
            for ( uint iFirst = 0; iFirst < aNumNodes; iFirst++ )
            {
                for ( uint iSecond = iFirst + 1; iSecond < aNumNodes; iSecond++ )
                {
                    if ( tUniform( aRandom ) < tProbability )
                    {
                        add_edge( tGraph, iFirst, iSecond );
                    }
                }
            }
            return tGraph;
        }

        uint
        max_degree( const Graph& aGraph )
        {
            uint tMax = 0;
            for ( const std::vector< moris_index >& iRow : aGraph )
            {
                tMax = std::max( tMax, (uint)iRow.size() );
            }
            return tMax;
        }

        bool
        same_labels( const Vector< moris_index >& aFirst, const Vector< moris_index >& aSecond )
        {
            return aFirst.size() == aSecond.size() and std::equal( aFirst.begin(), aFirst.end(), aSecond.begin() );
        }

        struct Fuzz_Outcome
        {
            uint mGraphs       = 0;
            uint mRun          = 0;    // real flood_fill called
            uint mSkipped      = 0;    // replica predicts a front overflow: the real call would write out of bounds
            uint mReplicaWrong = 0;    // replica labels differ from flood_fill (the replica is then not faithful)
            uint mWrong        = 0;    // flood_fill partition differs from the union-find reference
            uint mBadMax       = 0;    // aMaxValueAssigned is not the number of components - 1
        };

        // runs one graph through the replica, and through flood_fill if that is safe; checks against union-find
        void
        fuzz_one( const Graph& aGraph, const Vector< moris_index >& aPhases, uint aNumColumns, Fuzz_Outcome& aOutcome )
        {
            aOutcome.mGraphs++;
            Matrix< IndexMat >    tMatrix    = graph_matrix( aGraph, aNumColumns );
            Replica_Result        tReplica   = flood_fill_replica( tMatrix, aPhases );
            Vector< moris_index > tReference = reference_flood_fill( aGraph, aPhases );
            if ( tReplica.mOverflow )
            {
                aOutcome.mSkipped++;
                aOutcome.mWrong += not compare_partitions( tReference, tReplica.mLabels ).equal();
                return;
            }
            moris_index           tMax    = -1;
            Vector< moris_index > tLabels = real_flood_fill( tMatrix, aPhases, tMax );
            aOutcome.mRun++;
            aOutcome.mReplicaWrong += not same_labels( tLabels, tReplica.mLabels );
            aOutcome.mWrong += not compare_partitions( tReference, tLabels ).equal();
            aOutcome.mBadMax += tMax + 1 != (moris_index)count_distinct( tReference );
        }

        Vector< moris_index >
        random_phases( uint aNumNodes, uint aNumPhases, std::mt19937& aRandom )
        {
            Vector< moris_index > tPhases( aNumNodes );
            for ( uint iNode = 0; iNode < aNumNodes; iNode++ )
            {
                tPhases( iNode ) = aRandom() % aNumPhases;
            }
            return tPhases;
        }
    }    // namespace

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

    namespace
    {
        // ----------------------------------------------------------------------------------
        // geometries: up to six level-set fields, each one of a few shapes

        enum class Regime
        {
            STORED,
            ANALYTIC,
            MULTILINEAR
        };

        std::string
        regime_name( Regime aRegime )
        {
            switch ( aRegime )
            {
                case Regime::STORED:
                    return "stored";
                case Regime::ANALYTIC:
                    return "analytic";
                default:
                    return "multilinear";
            }
        }

        enum class Shape_Type
        {
            PLANE,        // a x + b y + c z + d
            SPHERE,       // | x - c | - r
            BOX_UNION,    // min over boxes of max_i ( | x_i - c_i | - h_i ): negative inside the union
            NODE_GRID     // trilinear interpolant of node values on a regular grid (gNodeGrid)
        };

        struct Shape
        {
            Shape_Type                     mType   = Shape_Type::PLANE;
            std::array< real, 4 >          mPlane  = { 0.0, 0.0, 0.0, 0.0 };
            Point                          mCentre = { 0.0, 0.0, 0.0 };
            real                           mRadius = 0.0;
            Vector< std::array< real, 6 > > mBoxes;    // ( xmin, xmax, ymin, ymax, zmin, zmax )
        };

        struct Node_Grid
        {
            std::array< uint, 3 > mCells   = { 1, 1, 1 };
            Point                 mOrigin  = { 0.0, 0.0, 0.0 };
            Point                 mSpacing = { 1.0, 1.0, 1.0 };
            Vector< real >        mValues;    // ( nx + 1 ) ( ny + 1 ) ( nz + 1 ), x fastest
        };

        const uint               gMaxShapes = 6;
        std::array< Shape, 6 >   gShapes    = {};
        Node_Grid                gNodeGrid;

        real
        node_grid_value( const Matrix< DDRMat >& aPoint )
        {
            std::array< uint, 3 > tIndex;
            Point                 tLocal;
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                real tScaled  = ( aPoint( iDim ) - gNodeGrid.mOrigin[ iDim ] ) / gNodeGrid.mSpacing[ iDim ];
                int  tCell    = std::min( std::max( (int)std::floor( tScaled ), 0 ), (int)gNodeGrid.mCells[ iDim ] - 1 );
                tIndex[ iDim ] = tCell;
                tLocal[ iDim ] = tScaled - tCell;
            }
            real tValue = 0.0;
            for ( uint iCorner = 0; iCorner < 8; iCorner++ )
            {
                real tWeight = 1.0;
                uint tNode   = 0;
                uint tStride = 1;
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    uint tBit = ( iCorner >> iDim ) & 1;
                    tWeight *= tBit ? tLocal[ iDim ] : 1.0 - tLocal[ iDim ];
                    tNode += ( tIndex[ iDim ] + tBit ) * tStride;
                    tStride *= gNodeGrid.mCells[ iDim ] + 1;
                }
                tValue += tWeight * gNodeGrid.mValues( tNode );
            }
            return tValue;
        }

        real
        shape_value( const Shape& aShape, const Matrix< DDRMat >& aPoint )
        {
            switch ( aShape.mType )
            {
                case Shape_Type::PLANE:
                    return aShape.mPlane[ 0 ] * aPoint( 0 ) + aShape.mPlane[ 1 ] * aPoint( 1 ) + aShape.mPlane[ 2 ] * aPoint( 2 ) + aShape.mPlane[ 3 ];
                case Shape_Type::SPHERE:
                {
                    real tSquared = 0.0;
                    for ( uint iDim = 0; iDim < 3; iDim++ )
                    {
                        tSquared += ( aPoint( iDim ) - aShape.mCentre[ iDim ] ) * ( aPoint( iDim ) - aShape.mCentre[ iDim ] );
                    }
                    return std::sqrt( tSquared ) - aShape.mRadius;
                }
                case Shape_Type::BOX_UNION:
                {
                    real tValue = MORIS_REAL_MAX;
                    for ( const std::array< real, 6 >& iBox : aShape.mBoxes )
                    {
                        real tBoxValue = -MORIS_REAL_MAX;
                        for ( uint iDim = 0; iDim < 3; iDim++ )
                        {
                            real tCentre = 0.5 * ( iBox[ 2 * iDim ] + iBox[ 2 * iDim + 1 ] );
                            real tHalf   = 0.5 * ( iBox[ 2 * iDim + 1 ] - iBox[ 2 * iDim ] );
                            tBoxValue    = std::max( tBoxValue, std::abs( aPoint( iDim ) - tCentre ) - tHalf );
                        }
                        tValue = std::min( tValue, tBoxValue );
                    }
                    return tValue;
                }
                default:
                    return node_grid_value( aPoint );
            }
        }

        template< uint N >
        real
        shape_field( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return shape_value( gShapes[ N ], aPoint );
        }

        const gen::Field_Function gShapeFields[ 6 ] = {
            &shape_field< 0 >, &shape_field< 1 >, &shape_field< 2 >, &shape_field< 3 >, &shape_field< 4 >, &shape_field< 5 >
        };

        std::shared_ptr< gen::Geometry >
        create_geometry( uint aSlot, Regime aRegime )
        {
            gen::Level_Set_Parameters tParameters;
            tParameters.mDiscretizationIndex         = aRegime == Regime::STORED ? -1 : -2;
            tParameters.mUseMultilinearInterpolation = aRegime == Regime::MULTILINEAR;

            std::shared_ptr< gen::Field > tField = std::make_shared< gen::User_Defined_Field >( gShapeFields[ aSlot ] );
            return std::make_shared< gen::Level_Set_Geometry >( tField, tParameters );
        }

        Shape
        plane( real aA, real aB, real aC, real aD )
        {
            Shape tShape;
            tShape.mType  = Shape_Type::PLANE;
            tShape.mPlane = { aA, aB, aC, aD };
            return tShape;
        }

        // unit normal through a point, shifted along the normal
        Shape
        plane_through( const Point& aNormal, const Point& aPoint, real aShift )
        {
            real  tNorm  = std::sqrt( aNormal[ 0 ] * aNormal[ 0 ] + aNormal[ 1 ] * aNormal[ 1 ] + aNormal[ 2 ] * aNormal[ 2 ] );
            Shape tShape = plane( 0.0, 0.0, 0.0, -aShift );
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tShape.mPlane[ iDim ] = aNormal[ iDim ] / tNorm;
                tShape.mPlane[ 3 ] -= tShape.mPlane[ iDim ] * aPoint[ iDim ];
            }
            return tShape;
        }

        Shape
        sphere( const Point& aCentre, real aRadius )
        {
            Shape tShape;
            tShape.mType   = Shape_Type::SPHERE;
            tShape.mCentre = aCentre;
            tShape.mRadius = aRadius;
            return tShape;
        }

        Shape
        box_union( const Vector< std::array< real, 6 > >& aBoxes )
        {
            Shape tShape;
            tShape.mType  = Shape_Type::BOX_UNION;
            tShape.mBoxes = aBoxes;
            return tShape;
        }

        Shape
        node_grid()
        {
            Shape tShape;
            tShape.mType = Shape_Type::NODE_GRID;
            return tShape;
        }

        // ----------------------------------------------------------------------------------
        // custom phase tables: two phases, "solid" (0) and "void" (1), from a predicate over the geometry signs.
        // Phase-table entry e: geometry g is positive iff bit ( n - 1 - g ) of e is set (gen::Phase_Table::get_phase_index)

        using Solid_Predicate = bool ( * )( uint aEntry, uint aNumGeometries );

        bool
        is_negative( uint aEntry, uint aGeometry, uint aNumGeometries )
        {
            return not( ( aEntry >> ( aNumGeometries - 1 - aGeometry ) ) & 1 );
        }

        // solid wherever some geometry is negative (union of the negative regions)
        bool
        solid_if_any_negative( uint aEntry, uint aNumGeometries )
        {
            for ( uint iGeom = 0; iGeom < aNumGeometries; iGeom++ )
            {
                if ( is_negative( aEntry, iGeom, aNumGeometries ) )
                {
                    return true;
                }
            }
            return false;
        }

        // solid where all geometries are negative or all are positive (two opposite quadrants / octants)
        bool
        solid_if_all_same_sign( uint aEntry, uint aNumGeometries )
        {
            bool tAllNegative = true;
            bool tAllPositive = true;
            for ( uint iGeom = 0; iGeom < aNumGeometries; iGeom++ )
            {
                tAllNegative = tAllNegative and is_negative( aEntry, iGeom, aNumGeometries );
                tAllPositive = tAllPositive and not is_negative( aEntry, iGeom, aNumGeometries );
            }
            return tAllNegative or tAllPositive;
        }

        // the first half of the geometries all negative, or the second half all positive (two overlapping opposite
        // quadrants / octants: the epsilon-overlap control of the kissing sweep)
        bool
        solid_if_first_half_negative_or_second_half_positive( uint aEntry, uint aNumGeometries )
        {
            uint tHalf        = aNumGeometries / 2;
            bool tAllNegative = true;
            bool tAllPositive = true;
            for ( uint iGeom = 0; iGeom < tHalf; iGeom++ )
            {
                tAllNegative = tAllNegative and is_negative( aEntry, iGeom, aNumGeometries );
                tAllPositive = tAllPositive and not is_negative( aEntry, tHalf + iGeom, aNumGeometries );
            }
            return tAllNegative or tAllPositive;
        }

        // solid except where the first geometry is positive and the second negative (a void shell between two spheres)
        bool
        solid_except_shell( uint aEntry, uint aNumGeometries )
        {
            return not( not is_negative( aEntry, 0, aNumGeometries ) and is_negative( aEntry, 1, aNumGeometries ) );
        }

        Matrix< DDUMat >
        phase_table( uint aNumGeometries, Solid_Predicate aSolid )
        {
            uint             tNumEntries = 1u << aNumGeometries;
            Matrix< DDUMat > tTable( 1, tNumEntries );
            for ( uint iEntry = 0; iEntry < tNumEntries; iEntry++ )
            {
                tTable( iEntry ) = aSolid( iEntry, aNumGeometries ) ? 0 : 1;
            }
            return tTable;
        }

        // ----------------------------------------------------------------------------------
        // background meshes: parameter list + perform() + Mesh_Manager (the setter path segfaults on > 1 rank)

        struct Grid_Mesh
        {
            std::shared_ptr< hmr::HMR >          mHMR;
            std::shared_ptr< mtk::Mesh_Manager > mMeshManager;
            mtk::Interpolation_Mesh*             mMesh = nullptr;

            // aRefineAboveX: background cells whose centre lies at x > aRefineAboveX are refined one level
            // aBsplineOrders / aBsplinePatterns / aLagrangeToBspline: HMR strings; aRefineFirstPattern: initial refinement of pattern 0
            // (pattern 1 unrefined). E.g. "1,2", "0,1", "0,1", 1: linear B-splines on the Lagrange mesh plus quadratic ones one level coarser
            Grid_Mesh(
                    const std::array< uint, 3 >& aElements,
                    const Point&                 aDimensions,
                    real                         aRefineAboveX      = MORIS_REAL_MAX,
                    const std::string&           aBsplineOrders     = "1",
                    const std::string&           aBsplinePatterns   = "0",
                    const std::string&           aLagrangeToBspline = "0",
                    uint                         aRefineFirstPattern = 0 )
            {
                Module_Parameter_Lists tParameters( Module_Type::HMR );
                tParameters( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );
                tParameters( 0 )( 0 ).set( "number_of_elements_per_dimension", aElements[ 0 ], aElements[ 1 ], aElements[ 2 ] );
                tParameters( 0 )( 0 ).set( "domain_dimensions", aDimensions[ 0 ], aDimensions[ 1 ], aDimensions[ 2 ] );
                tParameters( 0 )( 0 ).set( "domain_offset", 0.0, 0.0, 0.0 );
                tParameters( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
                tParameters( 0 )( 0 ).set( "lagrange_orders", "1" );
                tParameters( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
                tParameters( 0 )( 0 ).set( "bspline_orders", aBsplineOrders );
                tParameters( 0 )( 0 ).set( "bspline_pattern", aBsplinePatterns );
                tParameters( 0 )( 0 ).set( "lagrange_to_bspline", aLagrangeToBspline );
                tParameters( 0 )( 0 ).set( "refinement_buffer", 1 );
                tParameters( 0 )( 0 ).set( "staircase_buffer", 1 );
                if ( aRefineFirstPattern > 0 )
                {
                    tParameters( 0 )( 0 ).set( "pattern_initial_refinement", aRefineFirstPattern, 0u );    // pattern 0 refined, pattern 1 not
                }

                mHMR         = std::make_shared< hmr::HMR >( tParameters );
                mMeshManager = std::make_shared< mtk::Mesh_Manager >();
                mHMR->set_performer( mMeshManager );
                mHMR->perform_initial_refinement();
                if ( aRefineAboveX < MORIS_REAL_MAX )
                {
                    // flag as Database::flag_element does, without HMR::flag_element (MB-043: rank-local flag deadlocks)
                    hmr::Background_Mesh_Base* tBackgroundMesh = mHMR->get_database()->get_background_mesh();
                    luint                      tNumElements    = tBackgroundMesh->get_number_of_active_elements_on_proc();
                    Matrix< DDRMat >           tCentre( 3, 1 );
                    uint                       tWorkingPattern = mHMR->get_database()->get_parameters()->get_working_pattern();
                    for ( luint iElement = 0; iElement < tNumElements; iElement++ )
                    {
                        hmr::Background_Element_Base* tElement = tBackgroundMesh->get_element( iElement );
                        tBackgroundMesh->calc_center_of_element( tElement, tCentre );
                        if ( tCentre( 0 ) > aRefineAboveX )
                        {
                            tElement->put_on_refinement_queue();
                            tElement->set_refined_flag( tWorkingPattern );
                        }
                    }
                    mHMR->perform_refinement_based_on_working_pattern( 0 );
                }
                mHMR->perform();
                mMesh = mMeshManager->get_interpolation_mesh( 0 );
            }
        };

        // ----------------------------------------------------------------------------------
        // hashing (64-bit FNV-1a, split into two exactly representable reals for gathering as DDRMat)

        struct Hash
        {
            std::uint64_t mValue = 1469598103934665603ULL;

            void
            add( std::int64_t aValue )
            {
                for ( uint iByte = 0; iByte < 8; iByte++ )
                {
                    mValue ^= ( static_cast< std::uint64_t >( aValue ) >> ( 8 * iByte ) ) & 0xffULL;
                    mValue *= 1099511628211ULL;
                }
            }

            real
            high() const
            {
                return (real)( mValue >> 32 );
            }

            real
            low() const
            {
                return (real)( mValue & 0xffffffffULL );
            }
        };

        using Key = std::pair< std::uint64_t, std::uint64_t >;    // ( high, low ) of a hash, as integers

        Key
        key_of( const Hash& aHash )
        {
            return { aHash.mValue >> 32, aHash.mValue & 0xffffffffULL };
        }

        // rows gathered on rank 0 (empty elsewhere)
        Vector< Matrix< DDRMat > >
        gather_rows( const Vector< Vector< real > >& aRows, uint aNumColumns )
        {
            Matrix< DDRMat > tRows( aRows.size(), aNumColumns );
            for ( uint iRow = 0; iRow < aRows.size(); iRow++ )
            {
                for ( uint iCol = 0; iCol < aNumColumns; iCol++ )
                {
                    tRows( iRow, iCol ) = aRows( iRow )( iCol );
                }
            }
            Vector< Matrix< DDRMat > > tGathered;
            gatherv_mats( tRows, tGathered, 0 );
            return tGathered;
        }

        // ----------------------------------------------------------------------------------
        // geometry helpers for the reference

        struct Box
        {
            Point mMin = { MORIS_REAL_MAX, MORIS_REAL_MAX, MORIS_REAL_MAX };
            Point mMax = { -MORIS_REAL_MAX, -MORIS_REAL_MAX, -MORIS_REAL_MAX };

            real
            size( uint aDim ) const
            {
                return mMax[ aDim ] - mMin[ aDim ];
            }

            real
            volume() const
            {
                return this->size( 0 ) * this->size( 1 ) * this->size( 2 );
            }
        };

        std::array< std::int64_t, 3 >
        rounded_centroid( const Box& aBox )
        {
            std::array< std::int64_t, 3 > tCentroid;
            for ( uint iDim = 0; iDim < 3; iDim++ )
            {
                tCentroid[ iDim ] = std::llround( 0.5 * ( aBox.mMin[ iDim ] + aBox.mMax[ iDim ] ) * 1.0e9 );
            }
            return tCentroid;
        }

        Box
        cell_box( const mtk::Cell& aCell )
        {
            Box tBox;
            for ( mtk::Vertex* iVertex : aCell.get_vertex_pointers() )
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

        // the shared face of two background cells: the axis normal to it, its coordinate, and the overlap rectangle in the
        // two other axes ( ( a + 1 ) % 3, ( a + 2 ) % 3 )
        struct Shared_Face
        {
            uint                  mAxis  = 0;
            real                  mPlane = 0.0;
            std::array< real, 4 > mRectangle = { 0.0, 0.0, 0.0, 0.0 };    // umin, umax, vmin, vmax

            real
            area() const
            {
                return ( mRectangle[ 1 ] - mRectangle[ 0 ] ) * ( mRectangle[ 3 ] - mRectangle[ 2 ] );
            }
        };

        bool
        shared_face( const Box& aFirst, const Box& aSecond, Shared_Face& aFace )
        {
            real tScale = std::min( { aFirst.size( 0 ), aFirst.size( 1 ), aFirst.size( 2 ), aSecond.size( 0 ), aSecond.size( 1 ), aSecond.size( 2 ) } );
            real tTol   = 1.0e-9 * tScale;
            for ( uint iAxis = 0; iAxis < 3; iAxis++ )
            {
                real tPlane;
                if ( std::abs( aFirst.mMax[ iAxis ] - aSecond.mMin[ iAxis ] ) < tTol )
                {
                    tPlane = aFirst.mMax[ iAxis ];
                }
                else if ( std::abs( aSecond.mMax[ iAxis ] - aFirst.mMin[ iAxis ] ) < tTol )
                {
                    tPlane = aFirst.mMin[ iAxis ];
                }
                else
                {
                    continue;
                }
                uint tU = ( iAxis + 1 ) % 3;
                uint tV = ( iAxis + 2 ) % 3;
                std::array< real, 4 > tRectangle = { std::max( aFirst.mMin[ tU ], aSecond.mMin[ tU ] ), std::min( aFirst.mMax[ tU ], aSecond.mMax[ tU ] ),
                    std::max( aFirst.mMin[ tV ], aSecond.mMin[ tV ] ), std::min( aFirst.mMax[ tV ], aSecond.mMax[ tV ] ) };
                if ( tRectangle[ 1 ] - tRectangle[ 0 ] > tTol and tRectangle[ 3 ] - tRectangle[ 2 ] > tTol )
                {
                    aFace.mAxis      = iAxis;
                    aFace.mPlane     = tPlane;
                    aFace.mRectangle = tRectangle;
                    return true;
                }
            }
            return false;
        }

        using Polygon = std::vector< std::array< real, 2 > >;

        real
        polygon_area( const Polygon& aPolygon )
        {
            real tTwiceArea = 0.0;
            for ( uint iV = 0; iV < aPolygon.size(); iV++ )
            {
                const std::array< real, 2 >& tA = aPolygon[ iV ];
                const std::array< real, 2 >& tB = aPolygon[ ( iV + 1 ) % aPolygon.size() ];
                tTwiceArea += tA[ 0 ] * tB[ 1 ] - tA[ 1 ] * tB[ 0 ];
            }
            return 0.5 * std::abs( tTwiceArea );
        }

        // counter-clockwise copy
        Polygon
        counter_clockwise( const Polygon& aPolygon )
        {
            real tTwiceArea = 0.0;
            for ( uint iV = 0; iV < aPolygon.size(); iV++ )
            {
                const std::array< real, 2 >& tA = aPolygon[ iV ];
                const std::array< real, 2 >& tB = aPolygon[ ( iV + 1 ) % aPolygon.size() ];
                tTwiceArea += tA[ 0 ] * tB[ 1 ] - tA[ 1 ] * tB[ 0 ];
            }
            Polygon tResult = aPolygon;
            if ( tTwiceArea < 0.0 )
            {
                std::reverse( tResult.begin(), tResult.end() );
            }
            return tResult;
        }

        // Sutherland-Hodgman: aSubject clipped by the convex counter-clockwise aClip
        Polygon
        clip_polygon( const Polygon& aSubject, const Polygon& aClip )
        {
            Polygon tOutput = aSubject;
            for ( uint iEdge = 0; iEdge < aClip.size() and not tOutput.empty(); iEdge++ )
            {
                const std::array< real, 2 >& tA     = aClip[ iEdge ];
                const std::array< real, 2 >& tB     = aClip[ ( iEdge + 1 ) % aClip.size() ];
                Polygon                      tInput = tOutput;
                tOutput.clear();
                for ( uint iV = 0; iV < tInput.size(); iV++ )
                {
                    const std::array< real, 2 >& tP     = tInput[ iV ];
                    const std::array< real, 2 >& tQ     = tInput[ ( iV + 1 ) % tInput.size() ];
                    real                         tSideP = ( tB[ 0 ] - tA[ 0 ] ) * ( tP[ 1 ] - tA[ 1 ] ) - ( tB[ 1 ] - tA[ 1 ] ) * ( tP[ 0 ] - tA[ 0 ] );
                    real                         tSideQ = ( tB[ 0 ] - tA[ 0 ] ) * ( tQ[ 1 ] - tA[ 1 ] ) - ( tB[ 1 ] - tA[ 1 ] ) * ( tQ[ 0 ] - tA[ 0 ] );
                    if ( tSideP >= 0.0 )
                    {
                        tOutput.push_back( tP );
                    }
                    if ( ( tSideP >= 0.0 ) != ( tSideQ >= 0.0 ) )
                    {
                        real tT = tSideP / ( tSideP - tSideQ );
                        tOutput.push_back( { tP[ 0 ] + tT * ( tQ[ 0 ] - tP[ 0 ] ), tP[ 1 ] + tT * ( tQ[ 1 ] - tP[ 1 ] ) } );
                    }
                }
            }
            return tOutput;
        }

        Polygon
        rectangle_polygon( const std::array< real, 4 >& aRectangle )
        {
            return { { aRectangle[ 0 ], aRectangle[ 2 ] }, { aRectangle[ 1 ], aRectangle[ 2 ] }, { aRectangle[ 1 ], aRectangle[ 3 ] }, { aRectangle[ 0 ], aRectangle[ 3 ] } };
        }

        // the corner vertices of an IG cell (TET4, or the corners of a TET10)
        Vector< mtk::Vertex* >
        corner_vertices( const mtk::Cell& aCell )
        {
            Vector< mtk::Vertex* > tVertices = aCell.get_vertex_pointers();
            if ( tVertices.size() > 4 and aCell.get_geometry_type() == mtk::Geometry_Type::TET )
            {
                tVertices.resize( 4 );
            }
            return tVertices;
        }

        // ----------------------------------------------------------------------------------
        // section 2: reference checks of the three layers

        struct Connectivity_Counts
        {
            uint mCases = 0;    // when accumulated over a sweep: decomposed cases

            // layer 1: subphases
            uint mCutCells                  = 0;
            uint mSubphases                 = 0;
            uint mSubphasePartitionMismatch = 0;    // cut cells whose subphases differ from the reference
            uint mCoverage                  = 0;    // IG cells in no subphase, in two, or outside their cell's group
            uint mMixedPhaseSubphase        = 0;    // subphases whose IG cells differ in bulk phase
            uint mVolumeMismatch            = 0;    // background cells whose subphase volumes do not add up to the cell
            real mMaxVolumeError            = 0.0;

            // layer 2: neighbourhood
            uint mLinks                = 0;    // XTK links (directed)
            uint mReferenceLinks       = 0;    // reference links (directed)
            uint mMissingLinks         = 0;    // reference link that XTK lacks
            uint mExtraLinks           = 0;    // XTK link without reference contact (and not marginal)
            uint mExtraNoContact       = 0;    // ... of which the two subphases have no footprint overlap at all
            uint mExtraTransition      = 0;    // ... of which across a refinement transition
            uint mMarginalLinks        = 0;    // XTK links whose contact area is below the tolerance (degenerate facets)
            uint mSelfLinks            = 0;
            uint mDuplicateLinks       = 0;
            uint mAsymmetricLinks      = 0;
            uint mCrossPhaseLinks      = 0;
            uint mSameCellLinks        = 0;
            uint mNonAdjacentCells     = 0;    // XTK link between background cells that share no face
            uint mBadSideOrdinals      = 0;    // side ordinal whose vertices are not on the shared face
            uint mNonconformingFaces   = 0;    // same-level cut faces whose triangles overlap without matching
            uint mDuplicateVertexTriangles = 0;    // identical face triangles (by coordinates) built from different vertices
            uint mTransitionFaces      = 0;    // shared faces between cells of different refinement level

            // layer 3: enrichment
            uint mBases                = 0;
            uint mEnrichedBases        = 0;    // sum over bases of the XTK level count
            uint mReferenceComponents  = 0;    // sum over bases of the reference component count
            uint mCrossTalk            = 0;    // bases with a level spanning several reference components
            uint mSplitComponents      = 0;    // bases with a reference component split over several levels
            uint mLevelCountMismatch   = 0;    // bases whose level count differs from the component count
            uint mMissingBasis         = 0;    // subphase in support without an entry for the basis
            uint mEnrichedTotalMismatch = 0;   // mNumEnrichedBasisFunctions differs from the sum of level counts
            uint mZeroRowBases         = 0;    // enriched bases whose subphases have no positive integration volume (E3)
            uint mTinyRowBases         = 0;    // ... below 1e-12 of a background cell
            uint mMaxSupportSubphases  = 0;
            uint mMaxGraphDegree       = 0;    // widest row of the pruned subphase graph (50 columns allocated in Enrichment::generate_pruned_subphase_graph_in_basis_support)
            uint mFrontOverflow        = 0;    // pruned graphs whose flood-fill front would exceed 10 N + 1
            real mMaxFrontRatio        = 0.0;

            // bodies per phase, on this rank (only meaningful on one rank)
            uint mBodyMismatch = 0;    // phases whose XTK body count differs from the reference

            // ghost (E10)
            uint mGhostPairs         = 0;
            uint mGhostWithoutContact = 0;    // ghost facet between subphases that the reference does not link

            std::string mFirstFailure;

            void
            note( const std::string& aWhat )
            {
                if ( mFirstFailure.empty() )
                {
                    mFirstFailure = aWhat;
                }
            }

            uint
            failures() const
            {
                return mSubphasePartitionMismatch + mCoverage + mMixedPhaseSubphase + mVolumeMismatch + mMissingLinks + mExtraLinks + mSelfLinks
                     + mDuplicateLinks + mAsymmetricLinks + mCrossPhaseLinks + mSameCellLinks + mNonAdjacentCells + mBadSideOrdinals + mNonconformingFaces
                     + mCrossTalk + mSplitComponents + mLevelCountMismatch + mMissingBasis + mEnrichedTotalMismatch + mGhostWithoutContact;
            }

            void
            add( const Connectivity_Counts& aOther )
            {
                mCases += 1;
                mCutCells += aOther.mCutCells;
                mSubphases += aOther.mSubphases;
                mSubphasePartitionMismatch += aOther.mSubphasePartitionMismatch;
                mCoverage += aOther.mCoverage;
                mMixedPhaseSubphase += aOther.mMixedPhaseSubphase;
                mVolumeMismatch += aOther.mVolumeMismatch;
                mMaxVolumeError = std::max( mMaxVolumeError, aOther.mMaxVolumeError );
                mLinks += aOther.mLinks;
                mReferenceLinks += aOther.mReferenceLinks;
                mMissingLinks += aOther.mMissingLinks;
                mExtraLinks += aOther.mExtraLinks;
                mExtraNoContact += aOther.mExtraNoContact;
                mExtraTransition += aOther.mExtraTransition;
                mMarginalLinks += aOther.mMarginalLinks;
                mSelfLinks += aOther.mSelfLinks;
                mDuplicateLinks += aOther.mDuplicateLinks;
                mAsymmetricLinks += aOther.mAsymmetricLinks;
                mCrossPhaseLinks += aOther.mCrossPhaseLinks;
                mSameCellLinks += aOther.mSameCellLinks;
                mNonAdjacentCells += aOther.mNonAdjacentCells;
                mBadSideOrdinals += aOther.mBadSideOrdinals;
                mNonconformingFaces += aOther.mNonconformingFaces;
                mDuplicateVertexTriangles += aOther.mDuplicateVertexTriangles;
                mTransitionFaces += aOther.mTransitionFaces;
                mBases += aOther.mBases;
                mEnrichedBases += aOther.mEnrichedBases;
                mReferenceComponents += aOther.mReferenceComponents;
                mCrossTalk += aOther.mCrossTalk;
                mSplitComponents += aOther.mSplitComponents;
                mLevelCountMismatch += aOther.mLevelCountMismatch;
                mMissingBasis += aOther.mMissingBasis;
                mEnrichedTotalMismatch += aOther.mEnrichedTotalMismatch;
                mZeroRowBases += aOther.mZeroRowBases;
                mTinyRowBases += aOther.mTinyRowBases;
                mMaxSupportSubphases = std::max( mMaxSupportSubphases, aOther.mMaxSupportSubphases );
                mMaxGraphDegree      = std::max( mMaxGraphDegree, aOther.mMaxGraphDegree );
                mFrontOverflow += aOther.mFrontOverflow;
                mMaxFrontRatio = std::max( mMaxFrontRatio, aOther.mMaxFrontRatio );
                mBodyMismatch += aOther.mBodyMismatch;
                mGhostPairs += aOther.mGhostPairs;
                mGhostWithoutContact += aOther.mGhostWithoutContact;
            }

            // the same fields summed over the ranks (maxima reduced as maxima)
            void
            reduce()
            {
                Vector< uint* > tSums = { &mCutCells, &mSubphases, &mSubphasePartitionMismatch, &mCoverage, &mMixedPhaseSubphase, &mVolumeMismatch,
                    &mLinks, &mReferenceLinks, &mMissingLinks, &mExtraLinks, &mExtraNoContact, &mExtraTransition, &mMarginalLinks, &mSelfLinks,
                    &mDuplicateLinks, &mAsymmetricLinks, &mCrossPhaseLinks, &mSameCellLinks, &mNonAdjacentCells, &mBadSideOrdinals,
                    &mNonconformingFaces, &mDuplicateVertexTriangles, &mTransitionFaces, &mBases, &mEnrichedBases, &mReferenceComponents, &mCrossTalk, &mSplitComponents,
                    &mLevelCountMismatch, &mMissingBasis, &mEnrichedTotalMismatch, &mZeroRowBases, &mTinyRowBases, &mFrontOverflow, &mBodyMismatch,
                    &mGhostPairs, &mGhostWithoutContact };
                for ( uint* iValue : tSums )
                {
                    *iValue = sum_all( *iValue );
                }
                mMaxVolumeError      = max_all( mMaxVolumeError );
                mMaxSupportSubphases = max_all( mMaxSupportSubphases );
                mMaxGraphDegree      = max_all( mMaxGraphDegree );
                mMaxFrontRatio       = max_all( mMaxFrontRatio );
            }
        };

        // reference status of a pair of subphases
        enum class Contact
        {
            NONE,        // no overlap of the footprints
            MARGINAL,    // overlap area in ( 0, 1e-12 |face| ], or footprints made of degenerate triangles only
            LINK         // identical triangles on a conforming face, or overlap area above 1e-12 |face|
        };

        // a subphase's piece of a background face: a triangle of one of its IG cells, or the whole face (uncut cell)
        struct Footprint_Piece
        {
            moris_index                  mSubphase = 0;
            Polygon                      mPolygon;
            std::array< moris_index, 3 > mVertices  = { -1, -1, -1 };    // sorted IG vertex indices (triangles only)
            std::array< std::array< std::int64_t, 3 >, 3 > mCoordinates = {};    // sorted rounded vertex coordinates (triangles only)
            bool                         mWholeFace = false;
        };

        // data gathered for the multi-rank comparison
        struct Parallel_Rows
        {
            Vector< Vector< real > > mSubphases;           // subphase id, owned, key high, key low, background cell id, bulk phase
            Vector< Vector< real > > mEdges;               // owned subphase: key (high, low) of it and of the neighbour
            Vector< Vector< real > > mBases;               // owned basis: id, level count, hash of the levels (high, low)
            Vector< Vector< real > > mEnrichedIds;         // subphase id, basis id, enriched basis id, subphase owned
            Vector< Vector< real > > mClusterSubphases;    // subphase id of every cluster's interpolation cell on this rank
            Vector< Vector< real > > mOwnedSubphaseIds;    // owned subphase id
        };

        class Connectivity_Checker
        {
          public:
            mtk::Interpolation_Mesh* mBgMesh   = nullptr;
            Model*                   mModel    = nullptr;
            Cut_Integration_Mesh*    mCutMesh  = nullptr;
            moris_index              mMeshIndex = 0;
            uint                     mNumPhases = 0;

            Connectivity_Counts mCounts;

            Vector< Box >                                                  mBoxes;
            std::map< std::pair< moris_index, moris_index >, Shared_Face > mFaces;    // background cell pair (smaller first)
            std::map< std::pair< moris_index, moris_index >, Contact >     mContact;  // subphase pair (both orders)
            std::set< std::pair< moris_index, moris_index > >              mXtkLinks;
            Vector< Key >                                                  mSubphaseKeys;

            Connectivity_Checker( mtk::Interpolation_Mesh* aBgMesh, Model* aModel, uint aNumPhases, moris_index aMeshIndex )
                    : mBgMesh( aBgMesh )
                    , mModel( aModel )
                    , mCutMesh( aModel->get_cut_integration_mesh() )
                    , mMeshIndex( aMeshIndex )
                    , mNumPhases( aNumPhases )
            {
            }

            void
            check_all( const Vector< moris_index >& aMeshIndices, bool aEnrichment, bool aGhost, Parallel_Rows* aRows )
            {
                this->compute_boxes_and_faces();
                this->check_subphases();
                this->compute_subphase_keys();
                this->compute_reference_contacts();
                this->check_neighborhood();
                this->check_bodies();
                if ( aEnrichment )
                {
                    for ( uint iMesh = 0; iMesh < aMeshIndices.size(); iMesh++ )
                    {
                        mMeshIndex = aMeshIndices( iMesh );
                        this->check_enrichment( iMesh == 0 ? aRows : nullptr );
                    }
                }
                if ( aGhost )
                {
                    this->check_ghost();
                }
                if ( aRows != nullptr )
                {
                    this->collect_parallel_rows( *aRows, aEnrichment );
                }
            }

            // ------------------------------------------------------------------------------

            void
            compute_boxes_and_faces()
            {
                uint tNumCells = mBgMesh->get_num_elems();
                mBoxes.resize( tNumCells );
                for ( uint iCell = 0; iCell < tNumCells; iCell++ )
                {
                    mBoxes( iCell ) = cell_box( mBgMesh->get_mtk_cell( iCell ) );
                }
                for ( uint iFirst = 0; iFirst < tNumCells; iFirst++ )
                {
                    for ( uint iSecond = iFirst + 1; iSecond < tNumCells; iSecond++ )
                    {
                        Shared_Face tFace;
                        if ( shared_face( mBoxes( iFirst ), mBoxes( iSecond ), tFace ) )
                        {
                            mFaces[ { iFirst, iSecond } ] = tFace;
                        }
                    }
                }
            }

            // ------------------------------------------------------------------------------
            // layer 1

            void
            check_subphases()
            {
                uint tNumCells = mBgMesh->get_num_elems();
                for ( uint iCell = 0; iCell < tNumCells; iCell++ )
                {
                    Vector< moris_index > const & tSubphases = mCutMesh->get_parent_cell_subphases( iCell );
                    mCounts.mSubphases += tSubphases.size();

                    if ( not mCutMesh->parent_cell_has_children( iCell ) )
                    {
                        // an uncut cell is its own subphase with a single IG cell (itself)
                        bool tGood = tSubphases.size() == 1 and tSubphases( 0 ) == (moris_index)iCell
                                 and mCutMesh->get_subphase_ig_cells( iCell )->mIgCellGroup.size() == 1
                                 and mCutMesh->get_subphase_bulk_phase( iCell ) == mCutMesh->get_cell_bulk_phase( iCell );
                        if ( not tGood )
                        {
                            mCounts.mCoverage++;
                            mCounts.note( "uncut cell " + std::to_string( iCell ) + " is not a single subphase" );
                        }
                        continue;
                    }
                    mCounts.mCutCells++;

                    Vector< mtk::Cell* > const & tCells = mCutMesh->get_ig_cell_group( iCell )->mIgCellGroup;
                    std::map< moris_index, uint > tLocal;
                    for ( uint iIg = 0; iIg < tCells.size(); iIg++ )
                    {
                        tLocal[ tCells( iIg )->get_index() ] = iIg;
                    }

                    // XTK's partition of the cell's IG cells
                    Vector< moris_index > tXtk( tCells.size(), -1 );
                    real                  tVolume = 0.0;
                    for ( moris_index iSP : tSubphases )
                    {
                        moris_index tPhase = mCutMesh->get_subphase_bulk_phase( iSP );
                        bool        tMixed = false;
                        for ( mtk::Cell* iIgCell : mCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup )
                        {
                            std::map< moris_index, uint >::iterator tIt = tLocal.find( iIgCell->get_index() );
                            if ( tIt == tLocal.end() or tXtk( tIt->second ) != -1 )
                            {
                                mCounts.mCoverage++;
                                mCounts.note( "IG cell " + std::to_string( iIgCell->get_index() ) + " outside its group or in two subphases" );
                                continue;
                            }
                            tXtk( tIt->second ) = iSP;
                            tMixed              = tMixed or mCutMesh->get_cell_bulk_phase( iIgCell->get_index() ) != tPhase;
                            tVolume += iIgCell->compute_cell_measure();
                        }
                        mCounts.mMixedPhaseSubphase += tMixed;
                    }
                    for ( moris_index iLabel : tXtk )
                    {
                        if ( iLabel == -1 )
                        {
                            mCounts.mCoverage++;
                            mCounts.note( "IG cell of background cell " + std::to_string( iCell ) + " in no subphase" );
                        }
                    }

                    // reference: union-find over same-phase IG cells sharing a triangle
                    Union_Find tUnionFind( tCells.size() );
                    std::map< std::array< moris_index, 3 >, Vector< uint > > tCellsOfTriangle;
                    for ( uint iIg = 0; iIg < tCells.size(); iIg++ )
                    {
                        Vector< mtk::Vertex* > tCorners = corner_vertices( *tCells( iIg ) );
                        for ( uint iSkip = 0; iSkip < 4; iSkip++ )
                        {
                            std::array< moris_index, 3 > tTriangle;
                            uint                         tCount = 0;
                            for ( uint iV = 0; iV < 4; iV++ )
                            {
                                if ( iV != iSkip )
                                {
                                    tTriangle[ tCount++ ] = tCorners( iV )->get_index();
                                }
                            }
                            std::sort( tTriangle.begin(), tTriangle.end() );
                            tCellsOfTriangle[ tTriangle ].push_back( iIg );
                        }
                    }
                    for ( const std::pair< const std::array< moris_index, 3 >, Vector< uint > >& iTriangle : tCellsOfTriangle )
                    {
                        for ( uint iA = 0; iA < iTriangle.second.size(); iA++ )
                        {
                            for ( uint iB = iA + 1; iB < iTriangle.second.size(); iB++ )
                            {
                                uint tA = iTriangle.second( iA );
                                uint tB = iTriangle.second( iB );
                                if ( mCutMesh->get_cell_bulk_phase( tCells( tA )->get_index() ) == mCutMesh->get_cell_bulk_phase( tCells( tB )->get_index() ) )
                                {
                                    tUnionFind.unite( tA, tB );
                                }
                            }
                        }
                    }
                    if ( not compare_partitions( tUnionFind.labels(), tXtk ).equal() )
                    {
                        mCounts.mSubphasePartitionMismatch++;
                        mCounts.note( "subphases of background cell " + std::to_string( iCell ) + " differ from the reference" );
                    }

                    real tError = std::abs( tVolume - mBoxes( iCell ).volume() ) / mBoxes( iCell ).volume();
                    mCounts.mMaxVolumeError = std::max( mCounts.mMaxVolumeError, tError );
                    if ( tError > 1.0e-10 )
                    {
                        mCounts.mVolumeMismatch++;
                        mCounts.note( "subphase volumes of background cell " + std::to_string( iCell ) + " do not add up" );
                    }
                }
            }

            // ------------------------------------------------------------------------------
            // canonical subphase keys: background cell id, bulk phase, rounded IG cell centroids (independent of ranks)

            void
            compute_subphase_keys()
            {
                uint tNumSubphases = mCutMesh->get_num_subphases();
                mSubphaseKeys.resize( tNumSubphases );
                for ( uint iSP = 0; iSP < tNumSubphases; iSP++ )
                {
                    std::vector< std::array< std::int64_t, 3 > > tCentroids;
                    for ( mtk::Cell* iIgCell : mCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup )
                    {
                        Point                  tC        = { 0.0, 0.0, 0.0 };
                        Vector< mtk::Vertex* > tVertices = iIgCell->get_vertex_pointers();
                        for ( mtk::Vertex* iVertex : tVertices )
                        {
                            Matrix< DDRMat > tX = iVertex->get_coords();
                            for ( uint iDim = 0; iDim < 3; iDim++ )
                            {
                                tC[ iDim ] += tX( iDim ) / tVertices.size();
                            }
                        }
                        tCentroids.push_back( { std::llround( tC[ 0 ] * 1.0e9 ), std::llround( tC[ 1 ] * 1.0e9 ), std::llround( tC[ 2 ] * 1.0e9 ) } );
                    }
                    std::sort( tCentroids.begin(), tCentroids.end() );
                    Hash tHash;
                    // the parent cell by its centroid: background cell ids depend on the number of ranks
                    std::array< std::int64_t, 3 > tParent = rounded_centroid( mBoxes( mCutMesh->get_subphase_parent_cell( iSP )->get_index() ) );
                    tHash.add( tParent[ 0 ] );
                    tHash.add( tParent[ 1 ] );
                    tHash.add( tParent[ 2 ] );
                    tHash.add( mCutMesh->get_subphase_bulk_phase( iSP ) );
                    for ( const std::array< std::int64_t, 3 >& iC : tCentroids )
                    {
                        tHash.add( iC[ 0 ] );
                        tHash.add( iC[ 1 ] );
                        tHash.add( iC[ 2 ] );
                    }
                    mSubphaseKeys( iSP ) = key_of( tHash );
                }
            }

            // ------------------------------------------------------------------------------
            // layer 2 reference: footprints on shared faces

            Vector< Footprint_Piece >
            footprint( moris_index aCell, const Shared_Face& aFace )
            {
                Vector< Footprint_Piece > tPieces;
                uint                      tU = ( aFace.mAxis + 1 ) % 3;
                uint                      tV = ( aFace.mAxis + 2 ) % 3;
                if ( not mCutMesh->parent_cell_has_children( aCell ) )
                {
                    Footprint_Piece tPiece;
                    tPiece.mSubphase  = aCell;
                    tPiece.mPolygon   = rectangle_polygon( aFace.mRectangle );
                    tPiece.mWholeFace = true;
                    tPieces.push_back( tPiece );
                    return tPieces;
                }
                // vertices on a background face are interpolated between points on it, so they lie on it to round-off;
                // the tolerance must stay below the 1e-12 offsets of the near-face sweeps (MB-037 layers)
                real tTol = 1.0e-13 * mBoxes( aCell ).size( aFace.mAxis );
                for ( moris_index iSP : mCutMesh->get_parent_cell_subphases( aCell ) )
                {
                    for ( mtk::Cell* iIgCell : mCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup )
                    {
                        Vector< mtk::Vertex* > tCorners = corner_vertices( *iIgCell );
                        Vector< uint >         tOnFace;
                        for ( uint iV = 0; iV < tCorners.size(); iV++ )
                        {
                            if ( std::abs( tCorners( iV )->get_coords()( aFace.mAxis ) - aFace.mPlane ) <= tTol )
                            {
                                tOnFace.push_back( iV );
                            }
                        }
                        if ( tOnFace.size() != 3 )
                        {
                            continue;
                        }
                        Footprint_Piece tPiece;
                        tPiece.mSubphase = iSP;
                        for ( uint iK = 0; iK < 3; iK++ )
                        {
                            Matrix< DDRMat > tX = tCorners( tOnFace( iK ) )->get_coords();
                            tPiece.mPolygon.push_back( { tX( tU ), tX( tV ) } );
                            tPiece.mVertices[ iK ] = tCorners( tOnFace( iK ) )->get_index();
                            tPiece.mCoordinates[ iK ] = { std::llround( tX( 0 ) * 1.0e12 ), std::llround( tX( 1 ) * 1.0e12 ), std::llround( tX( 2 ) * 1.0e12 ) };
                        }
                        std::sort( tPiece.mVertices.begin(), tPiece.mVertices.end() );
                        std::sort( tPiece.mCoordinates.begin(), tPiece.mCoordinates.end() );
                        tPiece.mPolygon = counter_clockwise( tPiece.mPolygon );
                        tPieces.push_back( tPiece );
                    }
                }
                return tPieces;
            }

            void
            raise_contact( moris_index aFirst, moris_index aSecond, Contact aContact )
            {
                std::map< std::pair< moris_index, moris_index >, Contact >::iterator tIt = mContact.find( { aFirst, aSecond } );
                if ( tIt == mContact.end() or (int)tIt->second < (int)aContact )
                {
                    mContact[ { aFirst, aSecond } ] = aContact;
                    mContact[ { aSecond, aFirst } ] = aContact;
                }
            }

            void
            compute_reference_contacts()
            {
                for ( const std::pair< const std::pair< moris_index, moris_index >, Shared_Face >& iPair : mFaces )
                {
                    moris_index        tFirst  = iPair.first.first;
                    moris_index        tSecond = iPair.first.second;
                    const Shared_Face& tFace   = iPair.second;
                    bool               tSameLevel = mBgMesh->get_mtk_cell( tFirst ).get_level() == mBgMesh->get_mtk_cell( tSecond ).get_level();
                    mCounts.mTransitionFaces += not tSameLevel;

                    Vector< Footprint_Piece > tFirstPieces  = this->footprint( tFirst, tFace );
                    Vector< Footprint_Piece > tSecondPieces = this->footprint( tSecond, tFace );
                    real                      tFaceArea     = tFace.area();
                    Polygon                   tRegion       = rectangle_polygon( tFace.mRectangle );

                    // total contact area per subphase pair, and pairs that touch through identical triangles
                    std::map< std::pair< moris_index, moris_index >, real > tArea;
                    std::set< std::pair< moris_index, moris_index > >      tExact;
                    std::set< std::pair< moris_index, moris_index > >      tAnyPiece;
                    for ( const Footprint_Piece& iA : tFirstPieces )
                    {
                        for ( const Footprint_Piece& iB : tSecondPieces )
                        {
                            if ( mCutMesh->get_subphase_bulk_phase( iA.mSubphase ) != mCutMesh->get_subphase_bulk_phase( iB.mSubphase ) )
                            {
                                continue;
                            }
                            std::pair< moris_index, moris_index > tPair = { iA.mSubphase, iB.mSubphase };
                            tAnyPiece.insert( tPair );
                            bool tBothTriangles = not iA.mWholeFace and not iB.mWholeFace;
                            if ( tBothTriangles and tSameLevel )
                            {
                                // conforming face: contact is an identical triangle
                                if ( iA.mVertices == iB.mVertices )
                                {
                                    tExact.insert( tPair );
                                }
                                else if ( iA.mCoordinates == iB.mCoordinates )
                                {
                                    // the same triangle built from different (coincident) vertices: geometric contact, but
                                    // the IG mesh does not conform topologically here
                                    tExact.insert( tPair );
                                    mCounts.mDuplicateVertexTriangles++;
                                    this->dump_triangles( "duplicate-vertex", tFirst, tSecond, iA, iB );
                                }
                                else
                                {
                                    real tOverlap = polygon_area( clip_polygon( iA.mPolygon, iB.mPolygon ) );
                                    if ( tOverlap > 1.0e-10 * tFaceArea )
                                    {
                                        this->dump_triangles( "non-conforming", tFirst, tSecond, iA, iB );
                                        mCounts.mNonconformingFaces++;
                                        mCounts.note( "non-conforming triangles between cells " + std::to_string( tFirst ) + " and " + std::to_string( tSecond ) );
                                    }
                                }
                                continue;
                            }
                            Polygon tClipped = clip_polygon( clip_polygon( iA.mPolygon, tRegion ), iB.mPolygon );
                            tArea[ tPair ] += polygon_area( tClipped );
                        }
                    }
                    for ( const std::pair< moris_index, moris_index >& iSubphasePair : tAnyPiece )
                    {
                        Contact tContact = Contact::NONE;
                        if ( tExact.count( iSubphasePair ) )
                        {
                            tContact = Contact::LINK;
                        }
                        else if ( tArea.count( iSubphasePair ) )
                        {
                            real tA  = tArea[ iSubphasePair ];
                            tContact = tA > 1.0e-12 * tFaceArea ? Contact::LINK : tA > 0.0 ? Contact::MARGINAL : Contact::NONE;
                        }
                        // a subphase whose only pieces on the face are degenerate triangles: XTK sees facets there
                        bool tOneSideUncut = not mCutMesh->parent_cell_has_children( tFirst ) or not mCutMesh->parent_cell_has_children( tSecond );
                        if ( tContact == Contact::NONE and tSameLevel and tOneSideUncut
                                and this->has_only_degenerate_pieces( iSubphasePair, tFirstPieces, tSecondPieces, tFaceArea ) )
                        {
                            tContact = Contact::MARGINAL;
                        }
                        this->raise_contact( iSubphasePair.first, iSubphasePair.second, tContact );
                    }
                }
                for ( const std::pair< const std::pair< moris_index, moris_index >, Contact >& iContact : mContact )
                {
                    mCounts.mReferenceLinks += iContact.second == Contact::LINK;
                }
            }

            // MORIS_XTK_CONNECTIVITY_DUMP=1: print the first few face-triangle pairs that do not conform (indices, coordinates)
            uint mNumDumped = 0;

            void
            dump_triangles( const std::string& aWhat, moris_index aFirstCell, moris_index aSecondCell, const Footprint_Piece& aA, const Footprint_Piece& aB )
            {
                if ( std::getenv( "MORIS_XTK_CONNECTIVITY_DUMP" ) == nullptr or mNumDumped >= 6 )
                {
                    return;
                }
                mNumDumped++;
                std::ostringstream tStream;
                tStream << std::setprecision( 17 ) << aWhat << " pair, cells " << aFirstCell << " / " << aSecondCell << ", subphases " << aA.mSubphase << " / "
                        << aB.mSubphase << "\n";
                for ( const Footprint_Piece* iPiece : { &aA, &aB } )
                {
                    tStream << "   vertices";
                    for ( uint iV = 0; iV < 3; iV++ )
                    {
                        Matrix< DDRMat > tX = mCutMesh->get_mtk_vertex( iPiece->mVertices[ iV ] ).get_coords();
                        tStream << "  #" << iPiece->mVertices[ iV ] << " (" << tX( 0 ) << ", " << tX( 1 ) << ", " << tX( 2 ) << ")";
                    }
                    tStream << "\n";
                }
                std::cout << tStream.str();
            }

            // one side of the pair is an uncut whole face and the other side's pieces are all of (near) zero area
            bool
            has_only_degenerate_pieces(
                    const std::pair< moris_index, moris_index >& aPair,
                    const Vector< Footprint_Piece >&             aFirst,
                    const Vector< Footprint_Piece >&             aSecond,
                    real                                         aFaceArea )
            {
                bool tAnyTriangle = false;
                for ( const Vector< Footprint_Piece >* iSide : { &aFirst, &aSecond } )
                {
                    for ( const Footprint_Piece& iPiece : *iSide )
                    {
                        if ( ( iPiece.mSubphase == aPair.first or iPiece.mSubphase == aPair.second ) and not iPiece.mWholeFace )
                        {
                            tAnyTriangle = true;
                            if ( polygon_area( iPiece.mPolygon ) > 1.0e-12 * aFaceArea )
                            {
                                return false;
                            }
                        }
                    }
                }
                return tAnyTriangle;
            }

            Contact
            contact( moris_index aFirst, moris_index aSecond ) const
            {
                std::map< std::pair< moris_index, moris_index >, Contact >::const_iterator tIt = mContact.find( { aFirst, aSecond } );
                return tIt == mContact.end() ? Contact::NONE : tIt->second;
            }

            bool
            side_on_face( const mtk::Cell& aCell, moris_index aOrdinal, const Shared_Face& aFace )
            {
                if ( aOrdinal < 0 or aOrdinal >= 6 )
                {
                    return false;
                }
                real tTol = 1.0e-9 * mBoxes( aCell.get_index() ).size( aFace.mAxis );
                for ( mtk::Vertex const * iVertex : aCell.get_vertices_on_side_ordinal( aOrdinal ) )
                {
                    if ( std::abs( iVertex->get_coords()( aFace.mAxis ) - aFace.mPlane ) > tTol )
                    {
                        return false;
                    }
                }
                return true;
            }

            void
            check_neighborhood()
            {
                std::shared_ptr< Subphase_Neighborhood_Connectivity > tNeighborhood = mCutMesh->get_subphase_neighborhood();
                uint                                                  tNumSubphases = mCutMesh->get_num_subphases();
                for ( uint iSP = 0; iSP < tNumSubphases; iSP++ )
                {
                    Vector< moris_index > const & tNeighbors  = *tNeighborhood->mSubphaseToSubPhase( iSP );
                    Vector< moris_index > const & tMyOrds     = *tNeighborhood->mSubphaseToSubPhaseMySideOrds( iSP );
                    Vector< moris_index > const & tOtherOrds  = *tNeighborhood->mSubphaseToSubPhaseNeighborSideOrds( iSP );
                    Vector< moris_index > const & tTransition = *tNeighborhood->mTransitionNeighborCellLocation( iSP );
                    std::set< moris_index >       tSeen;
                    moris_index                   tMyCell = mCutMesh->get_subphase_parent_cell( iSP )->get_index();
                    for ( uint iN = 0; iN < tNeighbors.size(); iN++ )
                    {
                        moris_index tOther = tNeighbors( iN );
                        mCounts.mLinks++;
                        mXtkLinks.insert( { (moris_index)iSP, tOther } );
                        if ( tOther == (moris_index)iSP )
                        {
                            mCounts.mSelfLinks++;
                            mCounts.note( "self link of subphase " + std::to_string( iSP ) );
                            continue;
                        }
                        if ( not tSeen.insert( tOther ).second )
                        {
                            mCounts.mDuplicateLinks++;
                            mCounts.note( "duplicate link " + std::to_string( iSP ) + " -> " + std::to_string( tOther ) );
                        }
                        if ( mCutMesh->get_subphase_bulk_phase( iSP ) != mCutMesh->get_subphase_bulk_phase( tOther ) )
                        {
                            mCounts.mCrossPhaseLinks++;
                            mCounts.note( "cross-phase link " + std::to_string( iSP ) + " -> " + std::to_string( tOther ) );
                        }
                        moris_index tOtherCell = mCutMesh->get_subphase_parent_cell( tOther )->get_index();
                        if ( tOtherCell == tMyCell )
                        {
                            mCounts.mSameCellLinks++;
                            mCounts.note( "link inside background cell " + std::to_string( tMyCell ) );
                            continue;
                        }
                        std::map< std::pair< moris_index, moris_index >, Shared_Face >::iterator tFace =
                                mFaces.find( { std::min( tMyCell, tOtherCell ), std::max( tMyCell, tOtherCell ) } );
                        if ( tFace == mFaces.end() )
                        {
                            mCounts.mNonAdjacentCells++;
                            mCounts.note( "link between background cells " + std::to_string( tMyCell ) + " and " + std::to_string( tOtherCell ) + " that share no face" );
                            continue;
                        }
                        if ( not this->side_on_face( mBgMesh->get_mtk_cell( tMyCell ), tMyOrds( iN ), tFace->second )
                                or not this->side_on_face( mBgMesh->get_mtk_cell( tOtherCell ), tOtherOrds( iN ), tFace->second ) )
                        {
                            mCounts.mBadSideOrdinals++;
                            mCounts.note( "side ordinals of link " + std::to_string( iSP ) + " -> " + std::to_string( tOther ) + " not on the shared face" );
                        }
                        Contact tContact = this->contact( iSP, tOther );
                        if ( tContact == Contact::MARGINAL )
                        {
                            mCounts.mMarginalLinks++;
                        }
                        else if ( tContact == Contact::NONE )
                        {
                            mCounts.mExtraLinks++;
                            mCounts.mExtraNoContact += mContact.count( { (moris_index)iSP, tOther } ) == 0;
                            mCounts.mExtraTransition += tTransition( iN ) != MORIS_INDEX_MAX;
                            mCounts.note( "XTK links subphases " + std::to_string( iSP ) + " (cell " + std::to_string( tMyCell ) + ") and " + std::to_string( tOther )
                                          + " (cell " + std::to_string( tOtherCell ) + ") that do not touch"
                                          + ( tTransition( iN ) != MORIS_INDEX_MAX ? " [refinement transition]" : "" ) );
                        }
                    }
                }
                for ( const std::pair< moris_index, moris_index >& iLink : mXtkLinks )
                {
                    if ( iLink.first != iLink.second and mXtkLinks.count( { iLink.second, iLink.first } ) == 0 )
                    {
                        mCounts.mAsymmetricLinks++;
                        mCounts.note( "asymmetric link " + std::to_string( iLink.first ) + " -> " + std::to_string( iLink.second ) );
                    }
                }
                for ( const std::pair< const std::pair< moris_index, moris_index >, Contact >& iContact : mContact )
                {
                    if ( iContact.second == Contact::LINK and mXtkLinks.count( iContact.first ) == 0 )
                    {
                        mCounts.mMissingLinks++;
                        mCounts.note( "XTK lacks the link " + std::to_string( iContact.first.first ) + " -> " + std::to_string( iContact.first.second ) );
                    }
                }
            }

            // ------------------------------------------------------------------------------
            // bodies: connected components of all subphases, per phase, with the reference and with XTK's links

            void
            check_bodies()
            {
                uint       tNumSubphases = mCutMesh->get_num_subphases();
                Union_Find tReference( tNumSubphases );
                Union_Find tXtk( tNumSubphases );
                for ( const std::pair< const std::pair< moris_index, moris_index >, Contact >& iContact : mContact )
                {
                    if ( iContact.second != Contact::NONE )
                    {
                        tReference.unite( iContact.first.first, iContact.first.second );
                    }
                }
                for ( const std::pair< moris_index, moris_index >& iLink : mXtkLinks )
                {
                    tXtk.unite( iLink.first, iLink.second );
                }
                mReferenceBodies = Vector< uint >( mNumPhases, 0 );
                mXtkBodies       = Vector< uint >( mNumPhases, 0 );
                for ( uint iSP = 0; iSP < tNumSubphases; iSP++ )
                {
                    moris_index tPhase = mCutMesh->get_subphase_bulk_phase( iSP );
                    if ( tPhase < 0 or tPhase >= (moris_index)mNumPhases )
                    {
                        continue;
                    }
                    mReferenceBodies( tPhase ) += tReference.find( iSP ) == iSP;
                    mXtkBodies( tPhase ) += tXtk.find( iSP ) == iSP;
                }
                for ( uint iPhase = 0; iPhase < mNumPhases; iPhase++ )
                {
                    mCounts.mBodyMismatch += mReferenceBodies( iPhase ) != mXtkBodies( iPhase );
                }
            }

            Vector< uint > mReferenceBodies;
            Vector< uint > mXtkBodies;

            // ------------------------------------------------------------------------------
            // layer 3: enrichment levels per basis function

            void
            check_enrichment( Parallel_Rows* aRows )
            {
                Enrichment_Data const & tData        = mModel->get_basis_enrichment().get_enrichment_data()( mMeshIndex );
                uint                    tNumBases    = mBgMesh->get_num_basis_functions( mMeshIndex );
                moris_id                tRank        = par_rank();
                real                    tCellVolume  = mBoxes.size() > 0 ? mBoxes( 0 ).volume() : 1.0;
                std::shared_ptr< Subphase_Neighborhood_Connectivity > tNeighborhood = mCutMesh->get_subphase_neighborhood();

                // positive integration volume per subphase (FEM drops clusters with no positive cell volume, E3)
                uint           tNumSubphases = mCutMesh->get_num_subphases();
                Vector< real > tUsableVolume( tNumSubphases, 0.0 );
                for ( uint iSP = 0; iSP < tNumSubphases; iSP++ )
                {
                    for ( mtk::Cell* iIgCell : mCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup )
                    {
                        tUsableVolume( iSP ) += std::max( iIgCell->compute_cell_measure(), 0.0 );
                    }
                }

                uint tLevelTotal = 0;
                for ( uint iBasis = 0; iBasis < tNumBases; iBasis++ )
                {
                    mCounts.mBases++;
                    Matrix< IndexMat > tCellsInSupport;
                    mBgMesh->get_elements_in_support_of_basis( mMeshIndex, iBasis, tCellsInSupport );

                    Vector< moris_index >         tSubphases;
                    std::map< moris_index, uint > tLocal;
                    for ( uint iCell = 0; iCell < tCellsInSupport.numel(); iCell++ )
                    {
                        for ( moris_index iSP : mCutMesh->get_parent_cell_subphases( tCellsInSupport( iCell ) ) )
                        {
                            tLocal[ iSP ] = tSubphases.size();
                            tSubphases.push_back( iSP );
                        }
                    }
                    uint tNum                   = tSubphases.size();
                    mCounts.mMaxSupportSubphases = std::max( mCounts.mMaxSupportSubphases, tNum );

                    // reference components in the support
                    Union_Find tUnionFind( tNum );
                    for ( uint iA = 0; iA < tNum; iA++ )
                    {
                        for ( uint iB = iA + 1; iB < tNum; iB++ )
                        {
                            if ( this->contact( tSubphases( iA ), tSubphases( iB ) ) != Contact::NONE )
                            {
                                tUnionFind.unite( iA, iB );
                            }
                        }
                    }
                    Vector< moris_index > tReference  = tUnionFind.labels();
                    uint                  tComponents = count_distinct( tReference );
                    mCounts.mReferenceComponents += tComponents;

                    // XTK's levels
                    Vector< moris_index > tLevels( tNum, -1 );
                    for ( uint iLocal = 0; iLocal < tNum; iLocal++ )
                    {
                        Vector< moris_index > const & tBasisOfSP = tData.mSubphaseBGBasisIndices( tSubphases( iLocal ) );
                        Vector< moris_index > const & tLevelOfSP = tData.mSubphaseBGBasisEnrLev( tSubphases( iLocal ) );
                        for ( uint iEntry = 0; iEntry < tBasisOfSP.size(); iEntry++ )
                        {
                            if ( tBasisOfSP( iEntry ) == (moris_index)iBasis )
                            {
                                tLevels( iLocal ) = tLevelOfSP( iEntry );
                            }
                        }
                        if ( tLevels( iLocal ) == -1 )
                        {
                            mCounts.mMissingBasis++;
                            mCounts.note( "subphase " + std::to_string( tSubphases( iLocal ) ) + " has no level for basis " + std::to_string( iBasis ) );
                        }
                    }
                    uint tNumLevels = tData.mBasisEnrichmentIndices( iBasis ).numel();
                    tLevelTotal += tNumLevels;
                    mCounts.mEnrichedBases += tNumLevels;
                    if ( tNumLevels != tComponents )
                    {
                        mCounts.mLevelCountMismatch++;
                    }
                    Partition_Comparison tComparison = compare_partitions( tReference, tLevels );
                    if ( tComparison.mMerged > 0 )
                    {
                        mCounts.mCrossTalk++;
                        mCounts.note( "basis " + std::to_string( iBasis ) + ": an enrichment level spans " + "two reference components (cross-talk)" );
                    }
                    if ( tComparison.mSplit > 0 )
                    {
                        mCounts.mSplitComponents++;
                        mCounts.note( "basis " + std::to_string( iBasis ) + ": a connected component is split over two levels" );
                    }

                    // the pruned graph XTK flood-fills (Enrichment::generate_pruned_subphase_graph_in_basis_support):
                    // its widest row against the 50 allocated columns, and the replica's peak front against 10 N + 1
                    Graph tPruned( tNum );
                    for ( uint iLocal = 0; iLocal < tNum; iLocal++ )
                    {
                        for ( moris_index iNeighbor : *tNeighborhood->mSubphaseToSubPhase( tSubphases( iLocal ) ) )
                        {
                            std::map< moris_index, uint >::iterator tIt = tLocal.find( iNeighbor );
                            if ( tIt != tLocal.end() )
                            {
                                tPruned[ iLocal ].push_back( tIt->second );
                            }
                        }
                    }
                    uint tDegree           = max_degree( tPruned );
                    mCounts.mMaxGraphDegree = std::max( mCounts.mMaxGraphDegree, tDegree );
                    if ( tNum > 0 )
                    {
                        Replica_Result tReplica = flood_fill_replica( graph_matrix( tPruned, std::max( 50u, tDegree ) ), Vector< moris_index >( tNum, 1 ) );
                        mCounts.mFrontOverflow += tReplica.mOverflow;
                        mCounts.mMaxFrontRatio = std::max( mCounts.mMaxFrontRatio, real( tReplica.mPeakFront ) / tReplica.mFrontSize );
                    }

                    // E3: every enriched basis needs integration volume on some subphase of its level
                    std::map< moris_index, real > tLevelVolume;
                    for ( uint iLocal = 0; iLocal < tNum; iLocal++ )
                    {
                        tLevelVolume[ tLevels( iLocal ) ] += tUsableVolume( tSubphases( iLocal ) );
                    }
                    for ( const std::pair< const moris_index, real >& iLevel : tLevelVolume )
                    {
                        mCounts.mZeroRowBases += iLevel.second <= 0.0;
                        mCounts.mTinyRowBases += iLevel.second > 0.0 and iLevel.second < 1.0e-12 * tCellVolume;
                    }

                    if ( aRows != nullptr )
                    {
                        moris_id tBasisId    = mBgMesh->get_glb_entity_id_from_entity_loc_index( iBasis, mtk::EntityRank::BSPLINE, mMeshIndex );
                        uint     tBasisOwner = mBgMesh->get_entity_owner( iBasis, mtk::EntityRank::BSPLINE, mMeshIndex );
                        for ( uint iLocal = 0; iLocal < tNum; iLocal++ )
                        {
                            moris_index tSP = tSubphases( iLocal );
                            if ( tLevels( iLocal ) < 0 or tLevels( iLocal ) >= (moris_index)tNumLevels )
                            {
                                continue;
                            }
                            moris_index tEnriched = tData.mBasisEnrichmentIndices( iBasis )( tLevels( iLocal ) );
                            aRows->mEnrichedIds.push_back( { (real)mCutMesh->get_subphase_id( tSP ), (real)tBasisId, (real)tData.mEnrichedBasisIndexToId( tEnriched ),
                                    (real)( mCutMesh->get_subphase_parent_cell( tSP )->get_owner() == tRank ) } );
                        }
                        if ( tBasisOwner == (uint)tRank )
                        {
                            // the levels as sets of subphase keys, independent of the level numbering
                            std::map< moris_index, std::vector< Key > > tKeysOfLevel;
                            for ( uint iLocal = 0; iLocal < tNum; iLocal++ )
                            {
                                tKeysOfLevel[ tLevels( iLocal ) ].push_back( mSubphaseKeys( tSubphases( iLocal ) ) );
                            }
                            std::vector< Key > tLevelKeys;
                            for ( std::pair< const moris_index, std::vector< Key > >& iLevel : tKeysOfLevel )
                            {
                                std::sort( iLevel.second.begin(), iLevel.second.end() );
                                Hash tHash;
                                for ( const Key& iKey : iLevel.second )
                                {
                                    tHash.add( iKey.first );
                                    tHash.add( iKey.second );
                                }
                                tLevelKeys.push_back( key_of( tHash ) );
                            }
                            std::sort( tLevelKeys.begin(), tLevelKeys.end() );
                            Hash tHash;
                            for ( const Key& iKey : tLevelKeys )
                            {
                                tHash.add( iKey.first );
                                tHash.add( iKey.second );
                            }
                            // the basis by its support (rounded centroids of the support cells): background cell and basis
                            // ids depend on the number of ranks
                            std::vector< std::array< std::int64_t, 3 > > tSupportCentroids;
                            for ( uint iCell = 0; iCell < tCellsInSupport.numel(); iCell++ )
                            {
                                tSupportCentroids.push_back( rounded_centroid( mBoxes( tCellsInSupport( iCell ) ) ) );
                            }
                            std::sort( tSupportCentroids.begin(), tSupportCentroids.end() );
                            Hash tSupportHash;
                            for ( const std::array< std::int64_t, 3 >& iC : tSupportCentroids )
                            {
                                tSupportHash.add( iC[ 0 ] );
                                tSupportHash.add( iC[ 1 ] );
                                tSupportHash.add( iC[ 2 ] );
                            }
                            aRows->mBases.push_back( { tSupportHash.high(), tSupportHash.low(), (real)tNumLevels, tHash.high(), tHash.low() } );
                        }
                    }
                }
                if ( tLevelTotal != tData.mNumEnrichedBasisFunctions )
                {
                    mCounts.mEnrichedTotalMismatch++;
                    mCounts.note( "mNumEnrichedBasisFunctions " + std::to_string( tData.mNumEnrichedBasisFunctions ) + " != sum of level counts "
                                  + std::to_string( tLevelTotal ) );
                }
            }

            // ------------------------------------------------------------------------------
            // E10: every ghost facet joins two subphases that touch

            void
            check_ghost()
            {
                Enriched_Integration_Mesh& tEnrIgMesh = mModel->get_enriched_integ_mesh( 0 );
                for ( uint iSet = 0; iSet < tEnrIgMesh.get_num_sets(); iSet++ )
                {
                    mtk::Set* tSet = tEnrIgMesh.get_set_by_index( iSet );
                    if ( tSet->get_set_name().rfind( "ghost", 0 ) != 0 or tSet->get_set_type() != mtk::SetType::DOUBLE_SIDED_SIDESET )
                    {
                        continue;
                    }
                    for ( const mtk::Cluster* iCluster : tSet->get_clusters_on_set() )
                    {
                        Interpolation_Cell_Unzipped const * tLeader =
                                dynamic_cast< Interpolation_Cell_Unzipped const * >( &iCluster->get_interpolation_cell( mtk::Leader_Follower::LEADER ) );
                        Interpolation_Cell_Unzipped const * tFollower =
                                dynamic_cast< Interpolation_Cell_Unzipped const * >( &iCluster->get_interpolation_cell( mtk::Leader_Follower::FOLLOWER ) );
                        if ( tLeader == nullptr or tFollower == nullptr )
                        {
                            continue;
                        }
                        mCounts.mGhostPairs++;
                        if ( this->contact( tLeader->get_subphase_index(), tFollower->get_subphase_index() ) == Contact::NONE )
                        {
                            mCounts.mGhostWithoutContact++;
                            mCounts.note( "ghost facet between subphases " + std::to_string( tLeader->get_subphase_index() ) + " and "
                                          + std::to_string( tFollower->get_subphase_index() ) + " that do not touch" );
                        }
                    }
                }
            }

            // ------------------------------------------------------------------------------

            void
            collect_parallel_rows( Parallel_Rows& aRows, bool aEnrichment )
            {
                moris_id tRank         = par_rank();
                uint     tNumSubphases = mCutMesh->get_num_subphases();
                for ( uint iSP = 0; iSP < tNumSubphases; iSP++ )
                {
                    mtk::Cell const * tParent = mCutMesh->get_subphase_parent_cell( iSP );
                    bool              tOwned  = tParent->get_owner() == tRank;
                    aRows.mSubphases.push_back( { (real)mCutMesh->get_subphase_id( iSP ), (real)tOwned, (real)mSubphaseKeys( iSP ).first,
                            (real)mSubphaseKeys( iSP ).second, (real)tParent->get_id(), (real)mCutMesh->get_subphase_bulk_phase( iSP ) } );
                    if ( tOwned )
                    {
                        aRows.mOwnedSubphaseIds.push_back( { (real)mCutMesh->get_subphase_id( iSP ) } );
                    }
                }
                for ( const std::pair< moris_index, moris_index >& iLink : mXtkLinks )
                {
                    if ( mCutMesh->get_subphase_parent_cell( iLink.first )->get_owner() == tRank )
                    {
                        aRows.mEdges.push_back( { (real)mSubphaseKeys( iLink.first ).first, (real)mSubphaseKeys( iLink.first ).second,
                                (real)mSubphaseKeys( iLink.second ).first, (real)mSubphaseKeys( iLink.second ).second } );
                    }
                }
                if ( not aEnrichment )
                {
                    return;
                }
                Enriched_Integration_Mesh& tEnrIgMesh = mModel->get_enriched_integ_mesh( 0 );
                for ( uint iSet = 0; iSet < tEnrIgMesh.get_num_sets(); iSet++ )
                {
                    mtk::Set* tSet = tEnrIgMesh.get_set_by_index( iSet );
                    if ( tSet->get_set_type() != mtk::SetType::BULK )
                    {
                        continue;
                    }
                    for ( const mtk::Cluster* iCluster : tSet->get_clusters_on_set() )
                    {
                        Interpolation_Cell_Unzipped const * tCell =
                                dynamic_cast< Interpolation_Cell_Unzipped const * >( &iCluster->get_interpolation_cell( mtk::Leader_Follower::LEADER ) );
                        if ( tCell != nullptr )
                        {
                            aRows.mClusterSubphases.push_back( { (real)mCutMesh->get_subphase_id( tCell->get_subphase_index() ) } );
                        }
                    }
                }
            }
        };

        // ----------------------------------------------------------------------------------
        // multi-rank summary, evaluated on rank 0

        struct Parallel_Summary
        {
            uint mOwnedSubphases        = 0;
            uint mDuplicateOwned        = 0;    // owned subphase ids held as owned by two ranks
            uint mSubphaseCopyMismatch  = 0;    // ids whose copies on different ranks have different keys
            uint mEdges                 = 0;
            uint mOwnedBases            = 0;
            uint mDuplicateOwnedBases   = 0;
            uint mEnrichedIdMismatch    = 0;    // ( subphase id, basis id ) pairs whose enriched basis id differs between copies
            uint mClusterWithoutOwner   = 0;    // subphase ids used by a cluster that no rank owns
            Key  mSubphaseHash          = { 0, 0 };    // owned subphase keys
            Key  mSubphaseIdHash        = { 0, 0 };    // owned ( key, subphase id ): do the ids depend on the rank count?
            Key  mBgIdHash              = { 0, 0 };    // owned ( key, background cell id ): do the HMR ids depend on it?
            Key  mEdgeHash              = { 0, 0 };
            Key  mBasisHash             = { 0, 0 };
        };

        Key
        hash_rows( std::vector< std::vector< std::int64_t > >& aRows )
        {
            std::sort( aRows.begin(), aRows.end() );
            Hash tHash;
            for ( const std::vector< std::int64_t >& iRow : aRows )
            {
                for ( std::int64_t iValue : iRow )
                {
                    tHash.add( iValue );
                }
                tHash.add( -7 );
            }
            return key_of( tHash );
        }

        Parallel_Summary
        summarize_parallel( const Parallel_Rows& aRows )
        {
            Vector< Matrix< DDRMat > > tSubphases = gather_rows( aRows.mSubphases, 6 );
            Vector< Matrix< DDRMat > > tEdges     = gather_rows( aRows.mEdges, 4 );
            Vector< Matrix< DDRMat > > tBases     = gather_rows( aRows.mBases, 5 );
            Vector< Matrix< DDRMat > > tEnriched  = gather_rows( aRows.mEnrichedIds, 4 );
            Vector< Matrix< DDRMat > > tClusters  = gather_rows( aRows.mClusterSubphases, 1 );
            Vector< Matrix< DDRMat > > tOwnedIds  = gather_rows( aRows.mOwnedSubphaseIds, 1 );

            Parallel_Summary tSummary;
            if ( par_rank() != 0 )
            {
                return tSummary;
            }

            std::map< std::int64_t, Key > tKeyOfId;
            std::map< std::int64_t, uint > tOwnedCount;
            std::vector< std::vector< std::int64_t > > tOwnedKeys, tOwnedKeyIds, tOwnedKeyBgIds;
            for ( const Matrix< DDRMat >& iRows : tSubphases )
            {
                for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                {
                    std::int64_t tId  = (std::int64_t)iRows( iRow, 0 );
                    Key          tKey = { (std::uint64_t)iRows( iRow, 2 ), (std::uint64_t)iRows( iRow, 3 ) };
                    std::map< std::int64_t, Key >::iterator tIt = tKeyOfId.find( tId );
                    if ( tIt == tKeyOfId.end() )
                    {
                        tKeyOfId[ tId ] = tKey;
                    }
                    else if ( tIt->second != tKey )
                    {
                        tSummary.mSubphaseCopyMismatch++;
                    }
                    if ( iRows( iRow, 1 ) != 0.0 )
                    {
                        tOwnedCount[ tId ]++;
                        tOwnedKeys.push_back( { (std::int64_t)tKey.first, (std::int64_t)tKey.second } );
                        tOwnedKeyIds.push_back( { (std::int64_t)tKey.first, (std::int64_t)tKey.second, tId } );
                        tOwnedKeyBgIds.push_back( { (std::int64_t)tKey.first, (std::int64_t)tKey.second, (std::int64_t)iRows( iRow, 4 ) } );
                    }
                }
            }
            for ( const std::pair< const std::int64_t, uint >& iId : tOwnedCount )
            {
                tSummary.mDuplicateOwned += iId.second > 1;
            }
            tSummary.mOwnedSubphases = tOwnedKeys.size();
            tSummary.mSubphaseHash   = hash_rows( tOwnedKeys );
            tSummary.mSubphaseIdHash = hash_rows( tOwnedKeyIds );
            tSummary.mBgIdHash       = hash_rows( tOwnedKeyBgIds );

            std::vector< std::vector< std::int64_t > > tEdgeRows;
            for ( const Matrix< DDRMat >& iRows : tEdges )
            {
                for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                {
                    tEdgeRows.push_back( { (std::int64_t)iRows( iRow, 0 ), (std::int64_t)iRows( iRow, 1 ), (std::int64_t)iRows( iRow, 2 ), (std::int64_t)iRows( iRow, 3 ) } );
                }
            }
            tSummary.mEdges    = tEdgeRows.size();
            tSummary.mEdgeHash = hash_rows( tEdgeRows );

            std::vector< std::vector< std::int64_t > > tBasisRows;
            std::map< std::pair< std::int64_t, std::int64_t >, uint > tBasisCount;
            for ( const Matrix< DDRMat >& iRows : tBases )
            {
                for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                {
                    tBasisRows.push_back( { (std::int64_t)iRows( iRow, 0 ), (std::int64_t)iRows( iRow, 1 ), (std::int64_t)iRows( iRow, 2 ), (std::int64_t)iRows( iRow, 3 ),
                            (std::int64_t)iRows( iRow, 4 ) } );
                    tBasisCount[ { (std::int64_t)iRows( iRow, 0 ), (std::int64_t)iRows( iRow, 1 ) } ]++;
                }
            }
            for ( const std::pair< const std::pair< std::int64_t, std::int64_t >, uint >& iBasis : tBasisCount )
            {
                tSummary.mDuplicateOwnedBases += iBasis.second > 1;
            }
            tSummary.mOwnedBases = tBasisRows.size();
            tSummary.mBasisHash  = hash_rows( tBasisRows );

            std::map< std::pair< std::int64_t, std::int64_t >, std::int64_t > tEnrichedOf;
            std::set< std::pair< std::int64_t, std::int64_t > >               tBadPairs;
            for ( const Matrix< DDRMat >& iRows : tEnriched )
            {
                for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                {
                    std::pair< std::int64_t, std::int64_t > tPair = { (std::int64_t)iRows( iRow, 0 ), (std::int64_t)iRows( iRow, 1 ) };
                    std::int64_t tId = (std::int64_t)iRows( iRow, 2 );
                    std::map< std::pair< std::int64_t, std::int64_t >, std::int64_t >::iterator tIt = tEnrichedOf.find( tPair );
                    if ( tIt == tEnrichedOf.end() )
                    {
                        tEnrichedOf[ tPair ] = tId;
                    }
                    else if ( tIt->second != tId )
                    {
                        tBadPairs.insert( tPair );
                    }
                }
            }
            tSummary.mEnrichedIdMismatch = tBadPairs.size();

            std::set< std::int64_t > tOwned;
            for ( const Matrix< DDRMat >& iRows : tOwnedIds )
            {
                for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                {
                    tOwned.insert( (std::int64_t)iRows( iRow, 0 ) );
                }
            }
            for ( const Matrix< DDRMat >& iRows : tClusters )
            {
                for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                {
                    tSummary.mClusterWithoutOwner += tOwned.count( (std::int64_t)iRows( iRow, 0 ) ) == 0;
                }
            }
            return tSummary;
        }

        // ----------------------------------------------------------------------------------
        // one case: decompose, enrich, ghost, check

        struct Connectivity_Case
        {
            std::string      mName;
            Vector< Shape >  mShapes;
            Matrix< DDUMat > mBulkPhases;    // empty: the default 2^n phase table
        };

        struct Run_Options
        {
            bool                  mEnrichment           = true;
            bool                  mGhost                = true;
            bool                  mParallel             = false;    // gather the multi-rank rows (collective)
            bool                  mSortLevels           = false;    // sort_basis_enrichment_levels
            bool                  mTriangulateAllInPost = false;    // E6 (parameter-list Model)
            uint                  mIgOrder              = 1;        // E7: ig_element_order (parameter-list Model)
            Vector< moris_index > mMeshIndices          = { 0 };    // E8: enriched B-spline mesh indices
        };

        struct Case_Outcome
        {
            bool                mThrew      = false;    // on any rank
            bool                mGuard      = false;    // the MB-001 stop (expected)
            bool                mDecomposed = false;    // the level check passed on all ranks
            std::string         mMessage;
            Connectivity_Counts mCounts;
            Parallel_Summary    mParallel;
            Vector< uint >      mReferenceBodies;
            Vector< uint >      mXtkBodies;
            std::map< std::string, std::pair< uint, real > > mSets;    // rank 0: set name -> clusters, measure (summed)
        };

        const std::string gGuardMessage = "lies on a geometry interface, so its bulk phase cannot be determined";

        // per output set: clusters on this rank and their measure; gathered as name hashes like the aura sweep
        void
        collect_sets( Model& aModel, std::map< std::string, std::pair< uint, real > >& aSets )
        {
            Enriched_Integration_Mesh& tEnrIgMesh = aModel.get_enriched_integ_mesh( 0 );
            for ( uint iSet = 0; iSet < tEnrIgMesh.get_num_sets(); iSet++ )
            {
                mtk::Set*                tSet   = tEnrIgMesh.get_set_by_index( iSet );
                std::pair< uint, real >& tEntry = aSets[ tSet->get_set_name() ];
                for ( const mtk::Cluster* iCluster : tSet->get_clusters_on_set() )
                {
                    tEntry.first++;
                    tEntry.second += tSet->get_set_type() == mtk::SetType::BULK
                                           ? iCluster->compute_cluster_cell_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::LEADER )
                                           : iCluster->compute_cluster_cell_side_measure( mtk::Primary_Void::PRIMARY, mtk::Leader_Follower::LEADER );
                }
            }
        }

        std::map< std::string, std::pair< uint, real > >
        reduce_sets( const std::map< std::string, std::pair< uint, real > >& aLocal )
        {
            Vector< Vector< real > > tRows;
            Vector< std::string >    tNames;
            std::map< std::pair< real, real >, std::string > tNameOfHash;
            for ( const std::pair< const std::string, std::pair< uint, real > >& iSet : aLocal )
            {
                Hash tHash;
                for ( char iChar : iSet.first )
                {
                    tHash.add( iChar );
                }
                tRows.push_back( { tHash.high(), tHash.low(), (real)iSet.second.first, iSet.second.second } );
                tNameOfHash[ { tHash.high(), tHash.low() } ] = iSet.first;
            }
            Vector< Matrix< DDRMat > > tGathered = gather_rows( tRows, 4 );
            std::map< std::string, std::pair< uint, real > > tSets;
            for ( const Matrix< DDRMat >& iRows : tGathered )
            {
                for ( uint iRow = 0; iRow < iRows.n_rows(); iRow++ )
                {
                    std::pair< real, real > tKey = { iRows( iRow, 0 ), iRows( iRow, 1 ) };
                    std::string             tName;
                    if ( tNameOfHash.count( tKey ) )
                    {
                        tName = tNameOfHash[ tKey ];
                    }
                    else
                    {
                        std::ostringstream tStream;
                        tStream << "set#" << std::setprecision( 12 ) << tKey.first << "_" << tKey.second;
                        tName = tStream.str();
                    }
                    tSets[ tName ].first += (uint)iRows( iRow, 2 );
                    tSets[ tName ].second += iRows( iRow, 3 );
                }
            }
            return tSets;
        }

        Case_Outcome
        run_case( mtk::Interpolation_Mesh* aMesh, const Connectivity_Case& aCase, Regime aRegime, const Run_Options& aOptions )
        {
            Case_Outcome  tOutcome;
            Parallel_Rows tRows;
            std::map< std::string, std::pair< uint, real > > tLocalSets;
            bool          tThrew = false;
            try
            {
                MORIS_ERROR( aCase.mShapes.size() <= gMaxShapes, "run_case() - at most 6 geometries" );
                Vector< std::shared_ptr< gen::Geometry > > tGeometries;
                for ( uint iShape = 0; iShape < aCase.mShapes.size(); iShape++ )
                {
                    gShapes[ iShape ] = aCase.mShapes( iShape );
                    tGeometries.push_back( create_geometry( iShape, aRegime ) );
                }
                gen::Geometry_Engine_Parameters tParameters;
                tParameters.mGeometries = tGeometries;
                if ( aCase.mBulkPhases.numel() > 0 )
                {
                    tParameters.mBulkPhases = aCase.mBulkPhases;
                }
                gen::Geometry_Engine tGeometryEngine( aMesh, tParameters );

                std::ostringstream tMeshIndexList;
                Matrix< IndexMat > tMeshIndices( 1, aOptions.mMeshIndices.size() );
                for ( uint iMesh = 0; iMesh < aOptions.mMeshIndices.size(); iMesh++ )
                {
                    tMeshIndexList << ( iMesh ? "," : "" ) << aOptions.mMeshIndices( iMesh );
                    tMeshIndices( iMesh ) = aOptions.mMeshIndices( iMesh );
                }

                // the standalone Model, or (E6, E7) the parameter-list Model, which reads triangulate_all_in_post and
                // ig_element_order in perform_decomposition
                std::unique_ptr< Model > tModelPointer;
                if ( aOptions.mTriangulateAllInPost or aOptions.mIgOrder > 1 )
                {
                    Parameter_List tXtkParameters = prm::create_xtk_parameter_list();
                    tXtkParameters.set( "triangulate_all_in_post", aOptions.mTriangulateAllInPost );
                    tXtkParameters.set( "ig_element_order", aOptions.mIgOrder );
                    tXtkParameters.set( "enrich_mesh_indices", tMeshIndexList.str() );
                    tModelPointer = std::make_unique< Model >( tXtkParameters );
                    tModelPointer->set_geometry_engine( &tGeometryEngine );
                    tModelPointer->set_mtk_background_mesh( aMesh );
                    tOutcome.mDecomposed = tModelPointer->perform_decomposition();
                }
                else
                {
                    tModelPointer           = std::make_unique< Model >( 3, aMesh, &tGeometryEngine );
                    tModelPointer->mVerbose = false;
                    tOutcome.mDecomposed    = tModelPointer->decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );
                }
                Model& tModel = *tModelPointer;
                if ( tOutcome.mDecomposed )
                {
                    if ( aOptions.mEnrichment )
                    {
                        tModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, tMeshIndices, aOptions.mSortLevels );
                    }
                    if ( aOptions.mEnrichment and aOptions.mGhost )
                    {
                        tModel.construct_face_oriented_ghost_penalization_cells();
                    }
                    Connectivity_Checker tChecker( aMesh, &tModel, tGeometryEngine.get_num_bulk_phase(), 0 );
                    tChecker.check_all( aOptions.mMeshIndices, aOptions.mEnrichment, aOptions.mEnrichment and aOptions.mGhost, aOptions.mParallel ? &tRows : nullptr );
                    tOutcome.mCounts          = tChecker.mCounts;
                    tOutcome.mReferenceBodies = tChecker.mReferenceBodies;
                    tOutcome.mXtkBodies       = tChecker.mXtkBodies;
                    if ( aOptions.mParallel and aOptions.mEnrichment )
                    {
                        collect_sets( tModel, tLocalSets );
                    }
                }
            }
            catch ( const std::exception& aException )
            {
                tThrew            = true;
                tOutcome.mMessage = aException.what();
                tOutcome.mGuard   = tOutcome.mMessage.find( gGuardMessage ) != std::string::npos;
            }
            if ( aOptions.mParallel )
            {
                tOutcome.mThrew      = sum_all( (uint)tThrew ) > 0;
                tOutcome.mDecomposed = sum_all( (uint)tOutcome.mDecomposed ) == (uint)par_size();
                tOutcome.mCounts.reduce();
                tOutcome.mParallel = summarize_parallel( tRows );
                tOutcome.mSets     = reduce_sets( tLocalSets );
            }
            else
            {
                tOutcome.mThrew = tThrew;
            }
            return tOutcome;
        }

        // ----------------------------------------------------------------------------------
        // sweep bookkeeping and output

        struct Sweep_Totals
        {
            uint                mCases      = 0;
            uint                mThrew      = 0;
            uint                mGuard      = 0;
            uint                mRefused    = 0;    // decomposition refused by the level check
            uint                mFailing    = 0;    // cases with any reference failure
            uint                mBodyCases  = 0;    // cases where XTK's body count differs from the reference
            Connectivity_Counts mCounts;
            std::map< std::string, uint > mMessages;
            std::string         mFirstFailingCase;
            std::string         mFirstFailure;
            std::string         mFirstThrowCase;

            void
            add( const std::string& aCase, const Case_Outcome& aOutcome )
            {
                mCases++;
                if ( aOutcome.mThrew )
                {
                    ( aOutcome.mGuard ? mGuard : mThrew )++;
                    if ( not aOutcome.mGuard and mFirstThrowCase.empty() )
                    {
                        mFirstThrowCase = aCase;
                    }
                    mMessages[ aOutcome.mMessage.substr( 0, 160 ) ]++;
                    return;
                }
                if ( not aOutcome.mDecomposed )
                {
                    mRefused++;
                    return;
                }
                mCounts.add( aOutcome.mCounts );
                mBodyCases += aOutcome.mCounts.mBodyMismatch > 0;
                if ( aOutcome.mCounts.failures() > 0 )
                {
                    mFailing++;
                    if ( mFirstFailingCase.empty() )
                    {
                        mFirstFailingCase = aCase;
                        mFirstFailure     = aOutcome.mCounts.mFirstFailure;
                    }
                }
            }
        };

        void
        print_totals( const std::string& aTitle, const Sweep_Totals& aTotals )
        {
            if ( par_rank() != 0 )
            {
                return;
            }
            const Connectivity_Counts& tC = aTotals.mCounts;
            std::cout << "\n=== " << aTitle << " (" << par_size() << " rank(s)) ===\n"
                      << "  cases " << aTotals.mCases << " | threw " << aTotals.mThrew << " | MB-001 stops " << aTotals.mGuard << " | level check refused "
                      << aTotals.mRefused << " | checked " << tC.mCases << " | cases failing a reference check " << aTotals.mFailing << "\n"
                      << "  layer 1: cut cells " << tC.mCutCells << ", subphases " << tC.mSubphases << ", partition mismatch " << tC.mSubphasePartitionMismatch
                      << ", coverage " << tC.mCoverage << ", mixed-phase subphase " << tC.mMixedPhaseSubphase << ", volume mismatch " << tC.mVolumeMismatch
                      << " (max rel. error " << tC.mMaxVolumeError << ")\n"
                      << "  layer 2: XTK links " << tC.mLinks << ", reference links " << tC.mReferenceLinks << ", missing " << tC.mMissingLinks << ", extra "
                      << tC.mExtraLinks << " (no contact " << tC.mExtraNoContact << ", at a transition " << tC.mExtraTransition << "), marginal " << tC.mMarginalLinks
                      << ", self " << tC.mSelfLinks << ", duplicate " << tC.mDuplicateLinks << ", asymmetric " << tC.mAsymmetricLinks << ", cross-phase "
                      << tC.mCrossPhaseLinks << ", same cell " << tC.mSameCellLinks << ", non-adjacent cells " << tC.mNonAdjacentCells << ", bad ordinals "
                      << tC.mBadSideOrdinals << ", non-conforming faces " << tC.mNonconformingFaces << ", duplicate-vertex triangles " << tC.mDuplicateVertexTriangles << ", transition faces " << tC.mTransitionFaces << "\n"
                      << "  layer 3: bases " << tC.mBases << ", enriched bases " << tC.mEnrichedBases << ", reference components " << tC.mReferenceComponents
                      << ", cross-talk " << tC.mCrossTalk << ", split " << tC.mSplitComponents << ", level count mismatch " << tC.mLevelCountMismatch
                      << ", missing basis " << tC.mMissingBasis << ", total mismatch " << tC.mEnrichedTotalMismatch << "\n"
                      << "           zero-volume enriched bases (E3) " << tC.mZeroRowBases << ", below 1e-12 " << tC.mTinyRowBases << ", max subphases in a support "
                      << tC.mMaxSupportSubphases << ", max pruned-graph degree " << tC.mMaxGraphDegree << " (of 50), front overflow " << tC.mFrontOverflow
                      << ", max front / (10N+1) " << tC.mMaxFrontRatio << "\n"
                      << "  bodies: cases with XTK body count != reference " << aTotals.mBodyCases << " | ghost pairs " << tC.mGhostPairs << ", without contact "
                      << tC.mGhostWithoutContact << "\n";
            if ( not aTotals.mFirstFailingCase.empty() )
            {
                std::cout << "  first failing case: " << aTotals.mFirstFailingCase << ": " << aTotals.mFirstFailure << "\n";
            }
            if ( not aTotals.mFirstThrowCase.empty() )
            {
                std::cout << "  first unexpected throw: " << aTotals.mFirstThrowCase << "\n";
            }
            for ( const std::pair< const std::string, uint >& iMessage : aTotals.mMessages )
            {
                std::cout << "  throw x" << iMessage.second << ": " << iMessage.first << "\n";
            }
        }

        // CSV on rank 0 (MORIS_XTK_CONNECTIVITY_CSV), one row per case
        std::ofstream*
        case_csv()
        {
            static std::ofstream sStream;
            static bool          sOpened = false;
            if ( not sOpened )
            {
                sOpened = true;
                if ( const char* tPath = std::getenv( "MORIS_XTK_CONNECTIVITY_CSV" ) )
                {
                    if ( par_rank() == 0 )
                    {
                        sStream.open( tPath );
                        sStream << "num_ranks,sweep,regime,case,threw,guard,decomposed,failures,cut_cells,subphases,sp_mismatch,coverage,mixed_phase_sp,volume_mismatch,"
                                   "links,reference_links,missing_links,extra_links,extra_no_contact,extra_transition,marginal_links,self,duplicate,asymmetric,"
                                   "cross_phase,same_cell,non_adjacent,bad_ordinals,nonconforming,transition_faces,bases,enriched_bases,reference_components,"
                                   "cross_talk,split,level_count_mismatch,missing_basis,total_mismatch,zero_row,tiny_row,max_support_sp,max_degree,front_overflow,"
                                   "max_front_ratio,body_mismatch,ref_bodies,xtk_bodies,ghost_pairs,ghost_without_contact,owned_sp,dup_owned_sp,sp_copy_mismatch,"
                                   "sp_hash,sp_id_hash,bg_id_hash,edges,edge_hash,owned_bases,dup_owned_bases,basis_hash,enriched_id_mismatch,cluster_without_owner,"
                                   "first_failure,message\n";
                    }
                }
            }
            return sStream.is_open() ? &sStream : nullptr;
        }

        std::ofstream*
        sets_csv()
        {
            static std::ofstream sStream;
            static bool          sOpened = false;
            if ( not sOpened )
            {
                sOpened = true;
                if ( const char* tPath = std::getenv( "MORIS_XTK_CONNECTIVITY_SETS_CSV" ) )
                {
                    if ( par_rank() == 0 )
                    {
                        sStream.open( tPath );
                        sStream << "num_ranks,sweep,regime,case,set,clusters,measure\n";
                    }
                }
            }
            return sStream.is_open() ? &sStream : nullptr;
        }

        std::string
        csv_safe( std::string aText )
        {
            for ( char& iChar : aText )
            {
                if ( iChar == ',' or iChar == '\n' or iChar == '"' )
                {
                    iChar = ' ';
                }
            }
            return aText;
        }

        std::string
        key_string( const Key& aKey )
        {
            std::ostringstream tStream;
            tStream << std::hex << aKey.first << ":" << aKey.second;
            return tStream.str();
        }

        std::string
        bodies_string( const Vector< uint >& aBodies )
        {
            std::ostringstream tStream;
            for ( uint iPhase = 0; iPhase < aBodies.size(); iPhase++ )
            {
                tStream << ( iPhase ? "/" : "" ) << aBodies( iPhase );
            }
            return tStream.str();
        }

        void
        write_case( const std::string& aSweep, Regime aRegime, const std::string& aCase, const Case_Outcome& aOutcome )
        {
            std::ofstream* tStream = case_csv();
            if ( tStream != nullptr )
            {
                const Connectivity_Counts& tC = aOutcome.mCounts;
                const Parallel_Summary&    tP = aOutcome.mParallel;
                *tStream << par_size() << "," << aSweep << "," << regime_name( aRegime ) << "," << csv_safe( aCase ) << "," << aOutcome.mThrew << ","
                         << aOutcome.mGuard << "," << aOutcome.mDecomposed << "," << tC.failures() << "," << tC.mCutCells << "," << tC.mSubphases << ","
                         << tC.mSubphasePartitionMismatch << "," << tC.mCoverage << "," << tC.mMixedPhaseSubphase << "," << tC.mVolumeMismatch << ","
                         << tC.mLinks << "," << tC.mReferenceLinks << "," << tC.mMissingLinks << "," << tC.mExtraLinks << "," << tC.mExtraNoContact << ","
                         << tC.mExtraTransition << "," << tC.mMarginalLinks << "," << tC.mSelfLinks << "," << tC.mDuplicateLinks << "," << tC.mAsymmetricLinks
                         << "," << tC.mCrossPhaseLinks << "," << tC.mSameCellLinks << "," << tC.mNonAdjacentCells << "," << tC.mBadSideOrdinals << ","
                         << tC.mNonconformingFaces << "," << tC.mTransitionFaces << "," << tC.mBases << "," << tC.mEnrichedBases << "," << tC.mReferenceComponents
                         << "," << tC.mCrossTalk << "," << tC.mSplitComponents << "," << tC.mLevelCountMismatch << "," << tC.mMissingBasis << ","
                         << tC.mEnrichedTotalMismatch << "," << tC.mZeroRowBases << "," << tC.mTinyRowBases << "," << tC.mMaxSupportSubphases << ","
                         << tC.mMaxGraphDegree << "," << tC.mFrontOverflow << "," << tC.mMaxFrontRatio << "," << tC.mBodyMismatch << ","
                         << bodies_string( aOutcome.mReferenceBodies ) << "," << bodies_string( aOutcome.mXtkBodies ) << "," << tC.mGhostPairs << ","
                         << tC.mGhostWithoutContact << "," << tP.mOwnedSubphases << "," << tP.mDuplicateOwned << "," << tP.mSubphaseCopyMismatch << ","
                         << key_string( tP.mSubphaseHash ) << "," << key_string( tP.mSubphaseIdHash ) << "," << key_string( tP.mBgIdHash ) << "," << tP.mEdges << "," << key_string( tP.mEdgeHash )
                         << "," << tP.mOwnedBases << "," << tP.mDuplicateOwnedBases << "," << key_string( tP.mBasisHash ) << "," << tP.mEnrichedIdMismatch << ","
                         << tP.mClusterWithoutOwner << "," << csv_safe( tC.mFirstFailure.substr( 0, 200 ) ) << "," << csv_safe( aOutcome.mMessage.substr( 0, 160 ) )
                         << "\n"
                         << std::flush;
            }
            std::ofstream* tSets = sets_csv();
            if ( tSets != nullptr )
            {
                for ( const std::pair< const std::string, std::pair< uint, real > >& iSet : aOutcome.mSets )
                {
                    *tSets << par_size() << "," << aSweep << "," << regime_name( aRegime ) << "," << csv_safe( aCase ) << "," << iSet.first << ","
                           << iSet.second.first << "," << std::setprecision( 17 ) << iSet.second.second << std::setprecision( 6 ) << "\n";
                }
                *tSets << std::flush;
            }
        }

        // MORIS_XTK_CONNECTIVITY_CASE=<substring>: run only cases whose name contains it
        bool
        case_selected( const std::string& aName )
        {
            const char* tFilter = std::getenv( "MORIS_XTK_CONNECTIVITY_CASE" );
            return tFilter == nullptr or aName.find( tFilter ) != std::string::npos;
        }

        // MORIS_XTK_CONNECTIVITY_REGIMES=<stored,analytic,multilinear> (default: all)
        Vector< Regime >
        selected_regimes( const Vector< Regime >& aDefault )
        {
            const char* tList = std::getenv( "MORIS_XTK_CONNECTIVITY_REGIMES" );
            if ( tList == nullptr )
            {
                return aDefault;
            }
            Vector< Regime > tRegimes;
            std::string      tText( tList );
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
            {
                if ( tText.find( regime_name( iRegime ) ) != std::string::npos )
                {
                    tRegimes.push_back( iRegime );
                }
            }
            return tRegimes;
        }

        // runs a list of cases on a mesh, prints a line per case (if aVerbose) and the totals
        Sweep_Totals
        run_cases(
                const std::string&                  aSweep,
                mtk::Interpolation_Mesh*            aMesh,
                const Vector< Connectivity_Case >&  aCases,
                const Vector< Regime >&             aRegimes,
                const Run_Options&                  aOptions,
                bool                                aVerbose )
        {
            Sweep_Totals tTotals;
            if ( aVerbose and par_rank() == 0 )
            {
                std::cout << "\n--- " << aSweep << ", " << par_size() << " rank(s) ---\n";
            }
            for ( Regime iRegime : selected_regimes( aRegimes ) )
            {
                for ( const Connectivity_Case& iCase : aCases )
                {
                    if ( not case_selected( iCase.mName ) )
                    {
                        continue;
                    }
                    Case_Outcome tOutcome = run_case( aMesh, iCase, iRegime, aOptions );
                    tTotals.add( regime_name( iRegime ) + " | " + iCase.mName, tOutcome );
                    write_case( aSweep, iRegime, iCase.mName, tOutcome );
                    if ( aVerbose and par_rank() == 0 )
                    {
                        const Connectivity_Counts& tC = tOutcome.mCounts;
                        std::cout << std::setw( 12 ) << regime_name( iRegime ) << "  " << std::left << std::setw( 44 ) << iCase.mName << std::right;
                        if ( tOutcome.mThrew )
                        {
                            std::cout << ( tOutcome.mGuard ? "  MB-001 stop" : "  THREW: " + tOutcome.mMessage.substr( 0, 120 ) ) << "\n";
                            continue;
                        }
                        if ( not tOutcome.mDecomposed )
                        {
                            std::cout << "  refused by the level check\n";
                            continue;
                        }
                        std::cout << "  sp " << std::setw( 4 ) << tC.mSubphases << "  links " << std::setw( 4 ) << tC.mLinks << "/" << std::setw( 4 ) << tC.mReferenceLinks
                                  << "  levels " << std::setw( 4 ) << tC.mEnrichedBases << "/" << std::setw( 4 ) << tC.mReferenceComponents << "  bodies ref "
                                  << bodies_string( tOutcome.mReferenceBodies ) << " xtk " << bodies_string( tOutcome.mXtkBodies ) << "  fail " << tC.failures();
                        if ( tC.failures() > 0 )
                        {
                            std::cout << "  [" << tC.mFirstFailure.substr( 0, 140 ) << "]";
                        }
                        std::cout << "\n";
                    }
                }
            }
            return tTotals;
        }

        // where the processor boundaries are
        void
        print_ownership( mtk::Interpolation_Mesh* aMesh )
        {
            Point tMin = { MORIS_REAL_MAX, MORIS_REAL_MAX, MORIS_REAL_MAX };
            Point tMax = { -MORIS_REAL_MAX, -MORIS_REAL_MAX, -MORIS_REAL_MAX };
            uint  tOwned = 0;
            for ( uint iCell = 0; iCell < aMesh->get_num_elems(); iCell++ )
            {
                mtk::Cell& tCell = aMesh->get_mtk_cell( iCell );
                if ( tCell.get_owner() != (moris_id)par_rank() )
                {
                    continue;
                }
                tOwned++;
                Box tBox = cell_box( tCell );
                for ( uint iDim = 0; iDim < 3; iDim++ )
                {
                    tMin[ iDim ] = std::min( tMin[ iDim ], tBox.mMin[ iDim ] );
                    tMax[ iDim ] = std::max( tMax[ iDim ], tBox.mMax[ iDim ] );
                }
            }
            std::ostringstream tLine;
            tLine << "rank " << par_rank() << " of " << par_size() << " owns " << tOwned << " of " << aMesh->get_num_elems() << " held background cells in ["
                  << tMin[ 0 ] << "," << tMax[ 0 ] << "] x [" << tMin[ 1 ] << "," << tMax[ 1 ] << "] x [" << tMin[ 2 ] << "," << tMax[ 2 ] << "]\n";
            std::cout << tLine.str() << std::flush;
        }

        // ----------------------------------------------------------------------------------
        // case generators

        // single hex [0,1]^3: the node grid holds the 8 corner values
        void
        set_hex_corners( const std::array< real, 8 >& aCorners )
        {
            gNodeGrid.mCells   = { 1, 1, 1 };
            gNodeGrid.mOrigin  = { 0.0, 0.0, 0.0 };
            gNodeGrid.mSpacing = { 1.0, 1.0, 1.0 };
            gNodeGrid.mValues.clear();
            for ( real iValue : aCorners )
            {
                gNodeGrid.mValues.push_back( iValue );
            }
        }

        std::string
        values_string( const Vector< real >& aValues )
        {
            std::ostringstream tStream;
            for ( uint iV = 0; iV < aValues.size(); iV++ )
            {
                tStream << ( iV ? ";" : "" ) << aValues( iV );
            }
            return tStream.str();
        }

        // gap sweep: two same-phase half spaces separated by a void slab of width aWidth centred at aCentre (normal aNormal)
        Connectivity_Case
        gap_case( const std::string& aName, const Point& aNormal, const Point& aCentre, real aWidth )
        {
            Point tReversed = { -aNormal[ 0 ], -aNormal[ 1 ], -aNormal[ 2 ] };
            return { aName, { plane_through( aNormal, aCentre, -0.5 * aWidth ), plane_through( tReversed, aCentre, -0.5 * aWidth ) },
                phase_table( 2, &solid_if_any_negative ) };
        }

        // kissing at an edge: solid where x < a and y < b, or x > a and y > b (aOverlap > 0: the two quadrants overlap)
        Connectivity_Case
        edge_kissing_case( const std::string& aName, real aA, real aB, real aOverlap )
        {
            if ( aOverlap == 0.0 )
            {
                return { aName, { plane( 1.0, 0.0, 0.0, -aA ), plane( 0.0, 1.0, 0.0, -aB ) }, phase_table( 2, &solid_if_all_same_sign ) };
            }
            return { aName,
                { plane( 1.0, 0.0, 0.0, -aA - aOverlap ), plane( 0.0, 1.0, 0.0, -aB - aOverlap ), plane( 1.0, 0.0, 0.0, -aA + aOverlap ),
                        plane( 0.0, 1.0, 0.0, -aB + aOverlap ) },
                phase_table( 4, &solid_if_first_half_negative_or_second_half_positive ) };
        }

        // kissing at a vertex: solid in the octants ( -, -, - ) and ( +, +, + ) around aPoint
        Connectivity_Case
        vertex_kissing_case( const std::string& aName, const Point& aPoint, real aOverlap )
        {
            if ( aOverlap == 0.0 )
            {
                return { aName, { plane( 1.0, 0.0, 0.0, -aPoint[ 0 ] ), plane( 0.0, 1.0, 0.0, -aPoint[ 1 ] ), plane( 0.0, 0.0, 1.0, -aPoint[ 2 ] ) },
                    phase_table( 3, &solid_if_all_same_sign ) };
            }
            return { aName,
                { plane( 1.0, 0.0, 0.0, -aPoint[ 0 ] - aOverlap ), plane( 0.0, 1.0, 0.0, -aPoint[ 1 ] - aOverlap ), plane( 0.0, 0.0, 1.0, -aPoint[ 2 ] - aOverlap ),
                        plane( 1.0, 0.0, 0.0, -aPoint[ 0 ] + aOverlap ), plane( 0.0, 1.0, 0.0, -aPoint[ 1 ] + aOverlap ),
                        plane( 0.0, 0.0, 1.0, -aPoint[ 2 ] + aOverlap ) },
                phase_table( 6, &solid_if_first_half_negative_or_second_half_positive ) };
        }

        std::string
        number_string( real aValue )
        {
            std::ostringstream tStream;
            tStream << aValue;
            return tStream.str();
        }

        // gap cases on a box whose processor face (2 ranks) is x = aFace; positions inside a cell, on a background face
        // (x = aFace - h), on the processor face, and +-1e-11 off it
        Vector< Connectivity_Case >
        gap_cases( real aFace, real aCellSize, bool aSlanted )
        {
            Vector< std::pair< std::string, real > > tPositions = {
                { "in a cell", aFace - 0.37 * aCellSize },
                { "on a bg face", aFace - aCellSize },
                { "on the proc face", aFace },
                { "proc face +1e-11", aFace + 1.0e-11 },
                { "proc face -1e-11", aFace - 1.0e-11 }
            };
            Vector< Connectivity_Case > tCases;
            const Point                 tNormal = aSlanted ? Point{ 1.0, 0.3, 0.2 } : Point{ 1.0, 0.0, 0.0 };
            for ( const std::pair< std::string, real >& iPosition : tPositions )
            {
                for ( real tWidth : { 1.0e-1, 1.0e-2, 1.0e-3, 1.0e-6, 1.0e-9, 1.0e-11, 1.0e-12, 1.0e-13 } )
                {
                    tCases.push_back( gap_case( std::string( aSlanted ? "slanted gap " : "gap " ) + iPosition.first + ", w=" + number_string( tWidth ), tNormal,
                            { iPosition.second, 0.5, 0.5 }, tWidth ) );
                }
            }
            return tCases;
        }

        // kissing cases: edge contact along a line ( x = a, y = b ) and vertex contact at a point, inside a cell, on a
        // background face, on a background edge (diagonal across cells) and on the processor face, with an epsilon-overlap
        // control
        Vector< Connectivity_Case >
        kissing_cases( real aFace, real aCellSize )
        {
            real tIn = aFace - 0.63 * aCellSize;    // off the grid and off the subdivision vertices
            Vector< std::pair< std::string, std::array< real, 3 > > > tEdgeLines = {
                { "in a cell", { tIn, 0.37 * aCellSize, 0.0 } },
                { "x on a bg face", { aFace - aCellSize, 0.37 * aCellSize, 0.0 } },
                { "on a bg edge", { aFace - aCellSize, aCellSize, 0.0 } },
                { "x on the proc face", { aFace, 0.37 * aCellSize, 0.0 } },
                { "on a proc-face edge", { aFace, aCellSize, 0.0 } },
                { "proc face +1e-11", { aFace + 1.0e-11, 0.37 * aCellSize, 0.0 } }
            };
            Vector< Connectivity_Case > tCases;
            for ( const std::pair< std::string, std::array< real, 3 > >& iLine : tEdgeLines )
            {
                for ( real tOverlap : { 0.0, 1.0e-3, 1.0e-9 } )
                {
                    tCases.push_back( edge_kissing_case( "edge kiss " + iLine.first + ( tOverlap > 0.0 ? ", overlap " + number_string( tOverlap ) : "" ),
                            iLine.second[ 0 ], iLine.second[ 1 ], tOverlap ) );
                }
            }
            Vector< std::pair< std::string, Point > > tVertices = {
                { "in a cell", { tIn, 0.37 * aCellSize, 0.41 * aCellSize } },
                { "on a bg face", { aFace - aCellSize, 0.37 * aCellSize, 0.41 * aCellSize } },
                { "on a bg edge", { aFace - aCellSize, aCellSize, 0.41 * aCellSize } },
                { "on a bg vertex", { aFace - aCellSize, aCellSize, aCellSize } },
                { "on a proc-face vertex", { aFace, aCellSize, aCellSize } }
            };
            for ( const std::pair< std::string, Point >& iVertex : tVertices )
            {
                for ( real tOverlap : { 0.0, 1.0e-3 } )
                {
                    tCases.push_back( vertex_kissing_case( "vertex kiss " + iVertex.first + ( tOverlap > 0.0 ? ", overlap " + number_string( tOverlap ) : "" ),
                            iVertex.second, tOverlap ) );
                }
            }
            return tCases;
        }

        // processor-boundary cases on the 8 x 4 x 4 box [0,2] x [0,1] x [0,1]; processor face x = 1 on 2 ranks
        Vector< Connectivity_Case >
        parallel_cases()
        {
            Vector< Connectivity_Case > tCases;
            // P1: one region split by the processor face
            tCases.push_back( { "P1 sphere on the proc face", { sphere( { 1.0, 0.5, 0.5 }, 0.3 ) }, {} } );
            tCases.push_back( { "P1 sphere 1e-11 off the proc face", { sphere( { 1.0 + 1.0e-11, 0.5, 0.5 }, 0.3 ) }, {} } );
            tCases.push_back( { "P1 plane x=1 (proc face)", { plane( 1.0, 0.0, 0.0, -1.0 ) }, {} } );
            // P2: supports straddling ranks, several pieces per support
            tCases.push_back( { "P2 two spheres across the proc face", { sphere( { 0.85, 0.5, 0.5 }, 0.13 ), sphere( { 1.15, 0.5, 0.5 }, 0.13 ) },
                    phase_table( 2, &solid_if_any_negative ) } );
            tCases.push_back( { "P2 slanted plane", { plane_through( { 1.0, 0.3, 0.2 }, { 1.0, 0.5, 0.5 }, 0.0 ) }, {} } );
            // P3: two prongs on rank 0 (x < 1) joined only by a bridge owned by rank 1
            tCases.push_back( { "P3 U joined only on the other rank",
                    { box_union( { { 0.3, 1.2, 0.2, 0.45, 0.2, 0.8 }, { 0.3, 1.2, 0.55, 0.8, 0.2, 0.8 }, { 1.1, 1.4, 0.2, 0.8, 0.2, 0.8 } } ) }, {} } );
            tCases.push_back( { "P3 U with the bridge at x=1.3..1.7",
                    { box_union( { { 0.3, 1.45, 0.2, 0.45, 0.2, 0.8 }, { 0.3, 1.45, 0.55, 0.8, 0.2, 0.8 }, { 1.3, 1.7, 0.2, 0.8, 0.2, 0.8 } } ) }, {} } );
            tCases.push_back( { "P3 prongs only (no bridge)", { box_union( { { 0.3, 1.2, 0.2, 0.45, 0.2, 0.8 }, { 0.3, 1.2, 0.55, 0.8, 0.2, 0.8 } } ) }, {} } );
            // P4: a gap and a kissing contact on the processor face
            tCases.push_back( gap_case( "P4 gap w=1e-3 on the proc face", { 1.0, 0.0, 0.0 }, { 1.0, 0.5, 0.5 }, 1.0e-3 ) );
            tCases.push_back( gap_case( "P4 gap w=1e-11 on the proc face", { 1.0, 0.0, 0.0 }, { 1.0, 0.5, 0.5 }, 1.0e-11 ) );
            tCases.push_back( gap_case( "P4 gap w=1e-3 at proc face +1e-11", { 1.0, 0.0, 0.0 }, { 1.0 + 1.0e-11, 0.5, 0.5 }, 1.0e-3 ) );
            tCases.push_back( edge_kissing_case( "P4 edge kiss on the proc face", 1.0, 0.37, 0.0 ) );
            tCases.push_back( edge_kissing_case( "P4 edge kiss on a proc-face edge", 1.0, 0.5, 0.0 ) );
            tCases.push_back( vertex_kissing_case( "P4 vertex kiss on a proc-face vertex", { 1.0, 0.5, 0.5 }, 0.0 ) );
            tCases.push_back( vertex_kissing_case( "P4 vertex kiss, overlap 1e-3", { 1.0, 0.5, 0.5 }, 1.0e-3 ) );
            // E4: three quadrants of two crossing planes are one phase (no interface on the geometric surfaces inside it)
            tCases.push_back( { "E4 crossing planes, 3 quadrants one phase", { plane( 1.0, 0.0, 0.0, -1.0 ), plane( 0.0, 1.0, 0.0, -0.37 ) },
                    phase_table( 2, &solid_if_any_negative ) } );
            // E5: a region touching only the outer boundary, and a region pinched at the boundary
            tCases.push_back( { "E5 slab touching only the boundary x=0", { plane( 1.0, 0.0, 0.0, -0.1 ) }, {} } );
            tCases.push_back( { "E5 sphere pinched at the boundary y=0", { sphere( { 1.0, -0.3, 0.5 }, 0.3 ) }, {} } );
            tCases.push_back( { "E5 sphere just inside the boundary y=0", { sphere( { 1.0, -0.29, 0.5 }, 0.3 ) }, {} } );
            // E9: no cell cut, and every cell cut
            tCases.push_back( { "E9 no cell cut", { plane( 1.0, 0.0, 0.0, -3.0 ) }, {} } );
            tCases.push_back( { "E9 every cell cut", { node_grid() }, {} } );
            return tCases;
        }

        // node values ( -1 )^( i + j + k ) on the 9 x 5 x 5 nodes of the 8 x 4 x 4 box: every cell is cut (E9)
        void
        set_checkerboard_grid( const std::array< uint, 3 >& aCells, const Point& aSpacing )
        {
            gNodeGrid.mCells   = aCells;
            gNodeGrid.mOrigin  = { 0.0, 0.0, 0.0 };
            gNodeGrid.mSpacing = aSpacing;
            gNodeGrid.mValues.clear();
            for ( uint iK = 0; iK <= aCells[ 2 ]; iK++ )
            {
                for ( uint iJ = 0; iJ <= aCells[ 1 ]; iJ++ )
                {
                    for ( uint iI = 0; iI <= aCells[ 0 ]; iI++ )
                    {
                        // offset 0.25 keeps the cell-centre values away from 0
                        gNodeGrid.mValues.push_back( ( ( iI + iJ + iK ) % 2 ? -1.0 : 1.0 ) + 0.25 );
                    }
                }
            }
        }

        // random node values from { -1, 0, 1 } on a grid of aCells; returns false if some cell has all corners zero
        // (the MB-001 stop, which on several ranks could hang a collective)
        bool
        set_random_grid( const std::array< uint, 3 >& aCells, const Point& aSpacing, std::mt19937& aRandom, uint aPattern )
        {
            gNodeGrid.mCells   = aCells;
            gNodeGrid.mOrigin  = { 0.0, 0.0, 0.0 };
            gNodeGrid.mSpacing = aSpacing;
            gNodeGrid.mValues.clear();
            for ( uint iK = 0; iK <= aCells[ 2 ]; iK++ )
            {
                for ( uint iJ = 0; iJ <= aCells[ 1 ]; iJ++ )
                {
                    for ( uint iI = 0; iI <= aCells[ 0 ]; iI++ )
                    {
                        real tValue;
                        switch ( aPattern )
                        {
                            case 0:    // checkerboard
                                tValue = ( iI + iJ + iK ) % 2 ? -1.0 : 1.0;
                                break;
                            case 1:    // checkerboard in x-y, constant in z (saddle faces on every vertical face)
                                tValue = ( iI + iJ ) % 2 ? -1.0 : 1.0;
                                break;
                            case 2:    // checkerboard with zeros on one parity of z
                                tValue = iK % 2 ? 0.0 : ( ( iI + iJ ) % 2 ? -1.0 : 1.0 );
                                break;
                            default:
                                tValue = real( (int)( aRandom() % 3 ) - 1 );
                        }
                        gNodeGrid.mValues.push_back( tValue );
                    }
                }
            }
            for ( uint iK = 0; iK < aCells[ 2 ]; iK++ )
            {
                for ( uint iJ = 0; iJ < aCells[ 1 ]; iJ++ )
                {
                    for ( uint iI = 0; iI < aCells[ 0 ]; iI++ )
                    {
                        bool tAllZero = true;
                        for ( uint iCorner = 0; iCorner < 8; iCorner++ )
                        {
                            uint tNode = ( iI + ( iCorner & 1 ) ) + ( aCells[ 0 ] + 1 ) * ( ( iJ + ( ( iCorner >> 1 ) & 1 ) ) + ( aCells[ 1 ] + 1 ) * ( iK + ( ( iCorner >> 2 ) & 1 ) ) );
                            tAllZero   = tAllZero and gNodeGrid.mValues( tNode ) == 0.0;
                        }
                        if ( tAllZero )
                        {
                            return false;
                        }
                    }
                }
            }
            return true;
        }
    }    // namespace

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
