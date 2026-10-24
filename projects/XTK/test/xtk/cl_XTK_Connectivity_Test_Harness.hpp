/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_XTK_Connectivity_Test_Harness.hpp
 *
 * The reference machinery of the XTK connectivity tests (UT_XTK_Connectivity_Pathologies.cpp,
 * UT_XTK_Flood_Fill_Pathologies.cpp): union-find and partition comparison, a replica of the generic flood fill,
 * test graphs, geometries, the reference checker of the three connectivity layers, the case runner and the sweep
 * bookkeeping. Header-only: free functions and mutable globals are inline.
 *
 */

#pragma once

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
#include "cl_XTK_Pathology_Test_Tools.hpp"

namespace moris::xtk
{
    namespace connectivity_test
    {
        using namespace pathology_test;
        using namespace connectivity_check;

        // ----------------------------------------------------------------------------------
        // union-find and partition comparison

        // ----------------------------------------------------------------------------------
        // section 1: the generic flood fill (fn_mesh_flood_fill.hpp) on arbitrary graphs

        // connected components of the same-phase subgraph
        inline Vector< moris_index >
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

        // the real flood_fill, called as Enrichment::assign_subphase_bin_enrichment_levels_in_basis_support does
        inline Vector< moris_index >
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

        inline void
        add_edge( Graph& aGraph, uint aFirst, uint aSecond )
        {
            aGraph[ aFirst ].push_back( aSecond );
            aGraph[ aSecond ].push_back( aFirst );
        }

        inline Graph
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
        inline Graph
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

        inline bool
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
        inline void
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

        inline Vector< moris_index >
        random_phases( uint aNumNodes, uint aNumPhases, std::mt19937& aRandom )
        {
            Vector< moris_index > tPhases( aNumNodes );
            for ( uint iNode = 0; iNode < aNumNodes; iNode++ )
            {
                tPhases( iNode ) = aRandom() % aNumPhases;
            }
            return tPhases;
        }

        // ----------------------------------------------------------------------------------
        // geometries: up to six level-set fields, each one of a few shapes

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
        inline std::array< Shape, 6 >   gShapes    = {};
        inline Node_Grid                gNodeGrid;

        inline real
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

        inline real
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
        inline real
        shape_field( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return shape_value( gShapes[ N ], aPoint );
        }

        const gen::Field_Function gShapeFields[ 6 ] = {
            &shape_field< 0 >, &shape_field< 1 >, &shape_field< 2 >, &shape_field< 3 >, &shape_field< 4 >, &shape_field< 5 >
        };

        inline std::shared_ptr< gen::Geometry >
        create_geometry( uint aSlot, Regime aRegime )
        {
            gen::Level_Set_Parameters tParameters;
            tParameters.mDiscretizationIndex         = aRegime == Regime::STORED ? -1 : -2;
            tParameters.mUseMultilinearInterpolation = aRegime == Regime::MULTILINEAR;

            std::shared_ptr< gen::Field > tField = std::make_shared< gen::User_Defined_Field >( gShapeFields[ aSlot ] );
            return std::make_shared< gen::Level_Set_Geometry >( tField, tParameters );
        }

        inline Shape
        plane( real aA, real aB, real aC, real aD )
        {
            Shape tShape;
            tShape.mType  = Shape_Type::PLANE;
            tShape.mPlane = { aA, aB, aC, aD };
            return tShape;
        }

        // unit normal through a point, shifted along the normal
        inline Shape
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

        inline Shape
        sphere( const Point& aCentre, real aRadius )
        {
            Shape tShape;
            tShape.mType   = Shape_Type::SPHERE;
            tShape.mCentre = aCentre;
            tShape.mRadius = aRadius;
            return tShape;
        }

        inline Shape
        box_union( const Vector< std::array< real, 6 > >& aBoxes )
        {
            Shape tShape;
            tShape.mType  = Shape_Type::BOX_UNION;
            tShape.mBoxes = aBoxes;
            return tShape;
        }

        inline Shape
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

        inline bool
        is_negative( uint aEntry, uint aGeometry, uint aNumGeometries )
        {
            return not( ( aEntry >> ( aNumGeometries - 1 - aGeometry ) ) & 1 );
        }

        // solid wherever some geometry is negative (union of the negative regions)
        inline bool
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
        inline bool
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
        inline bool
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
        inline bool
        solid_except_shell( uint aEntry, uint aNumGeometries )
        {
            return not( not is_negative( aEntry, 0, aNumGeometries ) and is_negative( aEntry, 1, aNumGeometries ) );
        }

        inline Matrix< DDUMat >
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
        // geometry helpers for the reference

        // ----------------------------------------------------------------------------------
        // section 2: reference checks of the three layers

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

        inline Key
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

        inline Parallel_Summary
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
        inline void
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

        inline std::map< std::string, std::pair< uint, real > >
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

        inline Case_Outcome
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
                    tChecker.mMaxDumps = std::getenv( "MORIS_XTK_CONNECTIVITY_DUMP" ) != nullptr ? 6 : 0;    // MORIS_XTK_CONNECTIVITY_DUMP=1: log the first non-conforming face-triangle pairs
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

        inline void
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
        inline std::ofstream*
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

        inline std::ofstream*
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

        inline std::string
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

        inline std::string
        key_string( const Key& aKey )
        {
            std::ostringstream tStream;
            tStream << std::hex << aKey.first << ":" << aKey.second;
            return tStream.str();
        }

        inline std::string
        bodies_string( const Vector< uint >& aBodies )
        {
            std::ostringstream tStream;
            for ( uint iPhase = 0; iPhase < aBodies.size(); iPhase++ )
            {
                tStream << ( iPhase ? "/" : "" ) << aBodies( iPhase );
            }
            return tStream.str();
        }

        inline void
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
        inline bool
        case_selected( const std::string& aName )
        {
            const char* tFilter = std::getenv( "MORIS_XTK_CONNECTIVITY_CASE" );
            return tFilter == nullptr or aName.find( tFilter ) != std::string::npos;
        }

        // MORIS_XTK_CONNECTIVITY_REGIMES=<stored,analytic,multilinear> (default: all)
        inline Vector< Regime >
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
        inline Sweep_Totals
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
        inline void
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
        inline void
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

        inline std::string
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
        inline Connectivity_Case
        gap_case( const std::string& aName, const Point& aNormal, const Point& aCentre, real aWidth )
        {
            Point tReversed = { -aNormal[ 0 ], -aNormal[ 1 ], -aNormal[ 2 ] };
            return { aName, { plane_through( aNormal, aCentre, -0.5 * aWidth ), plane_through( tReversed, aCentre, -0.5 * aWidth ) },
                phase_table( 2, &solid_if_any_negative ) };
        }

        // kissing at an edge: solid where x < a and y < b, or x > a and y > b (aOverlap > 0: the two quadrants overlap)
        inline Connectivity_Case
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
        inline Connectivity_Case
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

        inline std::string
        number_string( real aValue )
        {
            std::ostringstream tStream;
            tStream << aValue;
            return tStream.str();
        }

        // gap cases on a box whose processor face (2 ranks) is x = aFace; positions inside a cell, on a background face
        // (x = aFace - h), on the processor face, and +-1e-11 off it
        inline Vector< Connectivity_Case >
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
        inline Vector< Connectivity_Case >
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
        inline Vector< Connectivity_Case >
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
        inline void
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
        inline bool
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
    }    // namespace connectivity_test
}    // namespace moris::xtk
