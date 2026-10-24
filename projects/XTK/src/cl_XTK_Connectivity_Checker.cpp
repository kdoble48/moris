/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Connectivity_Checker.cpp
 *
 */

#include "cl_XTK_Connectivity_Checker.hpp"

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <numeric>
#include <sstream>

#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_XTK_Enrichment.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"
#include "cl_XTK_Interpolation_Cell_Unzipped.hpp"
#include "cl_MTK_Vertex.hpp"
#include "cl_MTK_Set.hpp"
#include "cl_MTK_Cluster.hpp"

namespace moris::xtk::connectivity_check
{
    // ----------------------------------------------------------------------------------

    Connectivity_Checker::Connectivity_Checker( mtk::Interpolation_Mesh* aBgMesh, Model* aModel, uint aNumPhases, moris_index aMeshIndex )
            : mBgMesh( aBgMesh )
            , mModel( aModel )
            , mCutMesh( aModel->get_cut_integration_mesh() )
            , mMeshIndex( aMeshIndex )
            , mNumPhases( aNumPhases )
    {
    }

    // ----------------------------------------------------------------------------------

    uint
    count_distinct( const Vector< moris_index >& aLabels )
    {
        std::set< moris_index > tDistinct( aLabels.begin(), aLabels.end() );
        return tDistinct.size();
    }

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

    Key
    key_of( const Hash& aHash )
    {
        return { aHash.mValue >> 32, aHash.mValue & 0xffffffffULL };
    }

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

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

    // ----------------------------------------------------------------------------------

    Polygon
    rectangle_polygon( const std::array< real, 4 >& aRectangle )
    {
        return { { aRectangle[ 0 ], aRectangle[ 2 ] }, { aRectangle[ 1 ], aRectangle[ 2 ] }, { aRectangle[ 1 ], aRectangle[ 3 ] }, { aRectangle[ 0 ], aRectangle[ 3 ] } };
    }

    // ----------------------------------------------------------------------------------

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

    void
    Connectivity_Checker::check_all( const Vector< moris_index >& aMeshIndices, bool aEnrichment, bool aGhost, Parallel_Rows* aRows )
    {
        // the reference treats background cells as axis-aligned boxes and reads three coordinates
        MORIS_ERROR( mBgMesh->get_mesh_type() == mtk::MeshType::HMR and mBgMesh->get_spatial_dim() == 3,
                "Connectivity_Checker::check_all() - only 3D HMR background meshes are supported." );

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

    // ----------------------------------------------------------------------------------

    Connectivity_Counts
    Connectivity_Checker::check_and_report( const Vector< moris_index >& aMeshIndices, bool aEnrichment, bool aGhost, Report_Mode aMode )
    {
        bool tDetailed        = aMode == Report_Mode::DETAILED;
        mCounts.mMaxExamples  = tDetailed ? 10 : 1;
        mMaxDumps             = tDetailed ? std::max( mMaxDumps, 6u ) : mMaxDumps;

        this->check_all( aMeshIndices, aEnrichment, aGhost, nullptr );

        Connectivity_Counts tTotal = mCounts;
        tTotal.reduce();

        if ( tDetailed )
        {
            MORIS_LOG_INFO( "XTK connectivity check (detailed):\n%s", tTotal.detailed_report().c_str() );
            for ( const std::string& iExample : mCounts.mExamples )
            {
                MORIS_LOG_INFO_ALL_PROCS( "XTK connectivity check, rank %d: %s", par_rank(), iExample.c_str() );
            }
        }

        // the one line that is always written: INFO when it passes, WARNING (kept at higher severity levels) when it fails
        if ( tTotal.failures() == 0 )
        {
            MORIS_LOG_INFO( "XTK connectivity: %s", tTotal.sparse_report().c_str() );
        }
        else
        {
            MORIS_LOG_WARNING( "XTK connectivity: %s", tTotal.sparse_report().c_str() );
        }
        return tTotal;
    }

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::compute_boxes_and_faces()
    {
        uint tNumCells = mBgMesh->get_num_elems();
        mBoxes.resize( tNumCells );
        for ( uint iCell = 0; iCell < tNumCells; iCell++ )
        {
            mBoxes( iCell ) = cell_box( mBgMesh->get_mtk_cell( iCell ) );
        }
        for ( uint iFirst = 0; iFirst < tNumCells; iFirst++ )
        {
            // face neighbours, including the finer or coarser cells across a refinement transition
            Matrix< IndexMat > tNeighbors = mBgMesh->get_elements_connected_to_element_and_face_ord_loc_inds( iFirst );
            for ( uint iNeighbor = 0; iNeighbor < tNeighbors.n_cols(); iNeighbor++ )
            {
                moris_index tSecond = tNeighbors( 0, iNeighbor );
                Shared_Face tFace;
                if ( tSecond > (moris_index)iFirst and shared_face( mBoxes( iFirst ), mBoxes( tSecond ), tFace ) )
                {
                    mFaces[ { (moris_index)iFirst, tSecond } ] = tFace;
                }
            }
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::check_subphases()
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::compute_subphase_keys()
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

    // ----------------------------------------------------------------------------------

    Vector< Footprint_Piece >
    Connectivity_Checker::footprint( moris_index aCell, const Shared_Face& aFace )
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::raise_contact( moris_index aFirst, moris_index aSecond, Contact aContact )
    {
        std::map< std::pair< moris_index, moris_index >, Contact >::iterator tIt = mContact.find( { aFirst, aSecond } );
        if ( tIt == mContact.end() or (int)tIt->second < (int)aContact )
        {
            mContact[ { aFirst, aSecond } ] = aContact;
            mContact[ { aSecond, aFirst } ] = aContact;
        }
    }

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::compute_reference_contacts()
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::dump_triangles( const std::string& aWhat, moris_index aFirstCell, moris_index aSecondCell, const Footprint_Piece& aA, const Footprint_Piece& aB )
    {
        if ( mNumDumped >= mMaxDumps )
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
        MORIS_LOG_INFO_ALL_PROCS( "%s", tStream.str().c_str() );
    }

    // ----------------------------------------------------------------------------------

    bool
    Connectivity_Checker::has_only_degenerate_pieces(
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

    // ----------------------------------------------------------------------------------

    Contact
    Connectivity_Checker::contact( moris_index aFirst, moris_index aSecond ) const
    {
        std::map< std::pair< moris_index, moris_index >, Contact >::const_iterator tIt = mContact.find( { aFirst, aSecond } );
        return tIt == mContact.end() ? Contact::NONE : tIt->second;
    }

    // ----------------------------------------------------------------------------------

    bool
    Connectivity_Checker::side_on_face( const mtk::Cell& aCell, moris_index aOrdinal, const Shared_Face& aFace )
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::check_neighborhood()
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::check_bodies()
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::check_enrichment( Parallel_Rows* aRows )
    {
        Enrichment_Data const & tData        = mModel->get_basis_enrichment().get_enrichment_data()( mMeshIndex );
        uint                    tNumBases    = mBgMesh->get_num_basis_functions( mMeshIndex );
        moris_id                tRank        = par_rank();
        real                    tCellVolume  = mBoxes.size() > 0 ? mBoxes( 0 ).volume() : 1.0;
        std::shared_ptr< Subphase_Neighborhood_Connectivity > tNeighborhood = mCutMesh->get_subphase_neighborhood();

        // positive integration volume per subphase (FEM drops clusters with no positive cell volume, which leaves a zero row)
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::check_ghost()
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

    // ----------------------------------------------------------------------------------

    void
    Connectivity_Checker::collect_parallel_rows( Parallel_Rows& aRows, bool aEnrichment )
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

    // ----------------------------------------------------------------------------------

    std::string
    Connectivity_Counts::sparse_report() const
    {
        std::ostringstream tStream;
        if ( this->failures() == 0 )
        {
            tStream << "PASS cut_cells=" << mCutCells << " subphases=" << mSubphases << " links=" << mLinks << "/" << mReferenceLinks
                    << " enriched_bases=" << mEnrichedBases << "/" << mReferenceComponents << " ghost_facets=" << mGhostPairs;
            return tStream.str();
        }
        const std::pair< const char*, uint > tFailures[] = {
            { "subphase_partition", mSubphasePartitionMismatch }, { "ig_cell_coverage", mCoverage }, { "mixed_phase_subphases", mMixedPhaseSubphase },
            { "subphase_volumes", mVolumeMismatch }, { "missing_links", mMissingLinks }, { "extra_links", mExtraLinks }, { "self_links", mSelfLinks },
            { "duplicate_links", mDuplicateLinks }, { "asymmetric_links", mAsymmetricLinks }, { "cross_phase_links", mCrossPhaseLinks },
            { "same_cell_links", mSameCellLinks }, { "non_adjacent_links", mNonAdjacentCells }, { "bad_side_ordinals", mBadSideOrdinals },
            { "nonconforming_faces", mNonconformingFaces }, { "cross_talk", mCrossTalk }, { "split_components", mSplitComponents },
            { "level_count_mismatch", mLevelCountMismatch }, { "missing_bases", mMissingBasis }, { "enriched_total_mismatch", mEnrichedTotalMismatch },
            { "ghost_without_contact", mGhostWithoutContact }
        };
        tStream << "FAIL failures=" << this->failures();
        for ( const std::pair< const char*, uint >& iFailure : tFailures )
        {
            if ( iFailure.second > 0 )
            {
                tStream << " " << iFailure.first << "=" << iFailure.second;
            }
        }
        tStream << " (details: Report_Mode::DETAILED)";
        return tStream.str();
    }

    // ----------------------------------------------------------------------------------

    std::string
    Connectivity_Counts::detailed_report() const
    {
        std::ostringstream tStream;
        tStream << std::setprecision( 3 );
        tStream << "  subphases:     " << mSubphases << " in " << mCutCells << " cut cells; partition mismatches " << mSubphasePartitionMismatch
                << ", IG cell coverage " << mCoverage << ", mixed-phase " << mMixedPhaseSubphase << ", volume mismatches " << mVolumeMismatch
                << " (max error " << mMaxVolumeError << ")\n";
        tStream << "  neighbourhood: links " << mLinks << " (reference " << mReferenceLinks << "); missing " << mMissingLinks << ", extra " << mExtraLinks
                << " (no contact " << mExtraNoContact << ", across a transition " << mExtraTransition << "), marginal " << mMarginalLinks << "\n"
                << "                 self " << mSelfLinks << ", duplicate " << mDuplicateLinks << ", asymmetric " << mAsymmetricLinks << ", cross-phase "
                << mCrossPhaseLinks << ", same-cell " << mSameCellLinks << ", non-adjacent cells " << mNonAdjacentCells << ", bad side ordinals "
                << mBadSideOrdinals << "\n"
                << "                 transition faces " << mTransitionFaces << ", nonconforming faces " << mNonconformingFaces
                << ", duplicate-vertex triangles " << mDuplicateVertexTriangles << "\n";
        tStream << "  enrichment:    bases " << mBases << ", enriched " << mEnrichedBases << " (reference components " << mReferenceComponents << "); cross-talk "
                << mCrossTalk << ", split components " << mSplitComponents << ", level count mismatches " << mLevelCountMismatch << ", missing "
                << mMissingBasis << ", enriched total mismatch " << mEnrichedTotalMismatch << "\n"
                << "                 zero-volume bases " << mZeroRowBases << ", below 1e-12 of a cell " << mTinyRowBases << ", max subphases per support "
                << mMaxSupportSubphases << ", max graph degree " << mMaxGraphDegree << ", flood-fill front overflows " << mFrontOverflow
                << " (max front ratio " << mMaxFrontRatio << ")\n";
        tStream << "  bodies:        phases with a body-count mismatch " << mBodyMismatch << " (one rank only)\n";
        tStream << "  ghost:         facets " << mGhostPairs << ", without contact " << mGhostWithoutContact << "\n";
        tStream << "  result:        " << this->sparse_report();
        return tStream.str();
    }

}    // namespace moris::xtk::connectivity_check
