/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_From_Interface.cpp
 *
 */

#include "cl_SDF_From_Interface.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "fn_cross.hpp"
#include "op_minus.hpp"
#include "op_plus.hpp"
#include "op_times.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    namespace
    {
        //-------------------------------------------------------------------------------
        // Allocation-free kernels of the brute-force reference below (same arithmetic as
        // point_to_segment_distance / point_to_triangle_distance, but on raw doubles; they return the
        // SQUARED distance so that the nearest-facet search compares without a sqrt per facet).

        real
        clamp_01( real aValue )
        {
            return aValue < 0.0 ? 0.0 : ( aValue > 1.0 ? 1.0 : aValue );
        }

        real
        segment_distance_squared( const real* aP, const real* aV0, const real* aV1, uint aDim )
        {
            real tEdge[ 3 ] = { 0.0, 0.0, 0.0 };
            real tD[ 3 ]    = { 0.0, 0.0, 0.0 };
            real tEE        = 0.0;
            real tDE        = 0.0;
            for ( uint iDim = 0; iDim < aDim; ++iDim )
            {
                tEdge[ iDim ] = aV1[ iDim ] - aV0[ iDim ];
                tD[ iDim ]    = aP[ iDim ] - aV0[ iDim ];
                tEE += tEdge[ iDim ] * tEdge[ iDim ];
                tDE += tD[ iDim ] * tEdge[ iDim ];
            }

            // degenerate (zero-length) segment: distance to V0
            real tT = tEE < 1.0e-14 ? 0.0 : clamp_01( tDE / tEE );

            real tDistSq = 0.0;
            for ( uint iDim = 0; iDim < aDim; ++iDim )
            {
                real tDiff = aP[ iDim ] - ( aV0[ iDim ] + tT * tEdge[ iDim ] );
                tDistSq += tDiff * tDiff;
            }
            return tDistSq;
        }

        real
        triangle_distance_squared( const real* aP, const real* aV0, const real* aV1, const real* aV2 )
        {
            real tE0[ 3 ], tE1[ 3 ], tD[ 3 ];
            for ( uint iDim = 0; iDim < 3; ++iDim )
            {
                tE0[ iDim ] = aV1[ iDim ] - aV0[ iDim ];
                tE1[ iDim ] = aV2[ iDim ] - aV0[ iDim ];
                tD[ iDim ]  = aV0[ iDim ] - aP[ iDim ];
            }
            real tA  = tE0[ 0 ] * tE0[ 0 ] + tE0[ 1 ] * tE0[ 1 ] + tE0[ 2 ] * tE0[ 2 ];
            real tB  = tE0[ 0 ] * tE1[ 0 ] + tE0[ 1 ] * tE1[ 1 ] + tE0[ 2 ] * tE1[ 2 ];
            real tC  = tE1[ 0 ] * tE1[ 0 ] + tE1[ 1 ] * tE1[ 1 ] + tE1[ 2 ] * tE1[ 2 ];
            real tD0 = tE0[ 0 ] * tD[ 0 ] + tE0[ 1 ] * tD[ 1 ] + tE0[ 2 ] * tD[ 2 ];
            real tD1 = tE1[ 0 ] * tD[ 0 ] + tE1[ 1 ] * tD[ 1 ] + tE1[ 2 ] * tD[ 2 ];

            real tDet = tA * tC - tB * tB;
            real tS   = tB * tD1 - tC * tD0;
            real tT   = tB * tD0 - tA * tD1;

            // Voronoi regions, identical to SDF_From_Interface::point_to_triangle_distance
            if ( tS + tT <= tDet )
            {
                if ( tS < 0.0 )
                {
                    if ( tT < 0.0 )
                    {
                        tS = clamp_01( -tD0 / tA );
                        tT = 0.0;
                    }
                    else
                    {
                        tS = 0.0;
                        tT = clamp_01( -tD1 / tC );
                    }
                }
                else if ( tT < 0.0 )
                {
                    tS = clamp_01( -tD0 / tA );
                    tT = 0.0;
                }
                else
                {
                    real tInvDet = 1.0 / tDet;
                    tS *= tInvDet;
                    tT *= tInvDet;
                }
            }
            else
            {
                if ( tS < 0.0 )
                {
                    real tTmp0 = tB + tD0;
                    real tTmp1 = tC + tD1;
                    if ( tTmp1 > tTmp0 )
                    {
                        tS = clamp_01( ( tTmp1 - tTmp0 ) / ( tA - 2.0 * tB + tC ) );
                        tT = 1.0 - tS;
                    }
                    else
                    {
                        tS = 0.0;
                        tT = clamp_01( -tD1 / tC );
                    }
                }
                else if ( tT < 0.0 )
                {
                    real tTmp0 = tB + tD1;
                    real tTmp1 = tA + tD0;
                    if ( tTmp1 > tTmp0 )
                    {
                        tT = clamp_01( ( tTmp1 - tTmp0 ) / ( tA - 2.0 * tB + tC ) );
                        tS = 1.0 - tT;
                    }
                    else
                    {
                        tS = clamp_01( -tD0 / tA );
                        tT = 0.0;
                    }
                }
                else
                {
                    tS = clamp_01( ( ( tC + tD1 ) - ( tB + tD0 ) ) / ( tA - 2.0 * tB + tC ) );
                    tT = 1.0 - tS;
                }
            }

            real tDistSq = 0.0;
            for ( uint iDim = 0; iDim < 3; ++iDim )
            {
                real tDiff = aP[ iDim ] - ( aV0[ iDim ] + tS * tE0[ iDim ] + tT * tE1[ iDim ] );
                tDistSq += tDiff * tDiff;
            }
            return tDistSq;
        }

        //-------------------------------------------------------------------------------
        /**
         * Uniform bucket grid over the facets' bounding box for an EXACT nearest-facet query.
         * Each facet is listed in every cell its bounding box overlaps. A query visits the cells in rings of growing
         * Chebyshev index distance r around the point's (clamped) cell and stops once the best distance found is no
         * larger than a lower bound on the distance to any cell outside ring r (per dimension, the gap from the point to
         * the first unvisited cell layer). The result is the same minimum as the brute-force loop over all facets.
         */
        class Facet_Grid
        {
          public:
            Facet_Grid(
                    const Matrix< DDRMat >&   aFacetNodeCoords,
                    const Matrix< IndexMat >& aFacetConn,
                    uint                      aDim )
                    : mDim( aDim )
                    , mNodesPerFacet( aFacetConn.n_cols() )
                    , mNumFacets( aFacetConn.n_rows() )
            {
                // the kernels use the first two facet nodes (segment) unless the facet is a 3D triangle
                mUseTriangle = mDim == 3 && mNodesPerFacet == 3;
                uint tNodesUsed = mUseTriangle ? 3 : 2;

                // facet vertex coordinates, contiguous
                mCoords.assign( mNumFacets * tNodesUsed * mDim, 0.0 );
                mStride = tNodesUsed * mDim;
                for ( uint iFacet = 0; iFacet < mNumFacets; ++iFacet )
                {
                    for ( uint iNode = 0; iNode < tNodesUsed; ++iNode )
                    {
                        moris_index tIdx = aFacetConn( iFacet, iNode );
                        for ( uint iDim = 0; iDim < mDim; ++iDim )
                        {
                            mCoords[ iFacet * mStride + iNode * mDim + iDim ] = aFacetNodeCoords( tIdx, iDim );
                        }
                    }
                }

                // bounding box of all facets
                for ( uint iDim = 0; iDim < 3; ++iDim )
                {
                    mLo[ iDim ] = 0.0;
                    mN[ iDim ]  = 1;
                    mH[ iDim ]  = 1.0;
                }
                std::vector< real > tFLo( mNumFacets * mDim ), tFHi( mNumFacets * mDim );
                real tHi[ 3 ] = { 0.0, 0.0, 0.0 };
                for ( uint iDim = 0; iDim < mDim; ++iDim )
                {
                    mLo[ iDim ] = std::numeric_limits< real >::max();
                    tHi[ iDim ] = -std::numeric_limits< real >::max();
                }
                for ( uint iFacet = 0; iFacet < mNumFacets; ++iFacet )
                {
                    for ( uint iDim = 0; iDim < mDim; ++iDim )
                    {
                        real tMin = std::numeric_limits< real >::max();
                        real tMax = -std::numeric_limits< real >::max();
                        for ( uint iNode = 0; iNode < tNodesUsed; ++iNode )
                        {
                            real tX = mCoords[ iFacet * mStride + iNode * mDim + iDim ];
                            tMin    = std::min( tMin, tX );
                            tMax    = std::max( tMax, tX );
                        }
                        tFLo[ iFacet * mDim + iDim ] = tMin;
                        tFHi[ iFacet * mDim + iDim ] = tMax;
                        mLo[ iDim ]                  = std::min( mLo[ iDim ], tMin );
                        tHi[ iDim ]                  = std::max( tHi[ iDim ], tMax );
                    }
                }

                // about one facet per cell: cell size = (box measure / facets)^(1/dim), at least 1 cell per dimension
                real tMeasure = 1.0;
                real tMaxExt  = 0.0;
                for ( uint iDim = 0; iDim < mDim; ++iDim )
                {
                    tMaxExt = std::max( tMaxExt, tHi[ iDim ] - mLo[ iDim ] );
                }
                real tFloor = std::max( tMaxExt * 1.0e-6, 1.0e-12 );
                for ( uint iDim = 0; iDim < mDim; ++iDim )
                {
                    tMeasure *= std::max( tHi[ iDim ] - mLo[ iDim ], tFloor );
                }
                real tCell = std::pow( tMeasure / std::max( (real)mNumFacets, 1.0 ), 1.0 / mDim );
                uint tNumCells = 1;
                for ( uint iDim = 0; iDim < mDim; ++iDim )
                {
                    real tExt = std::max( tHi[ iDim ] - mLo[ iDim ], tFloor );
                    mN[ iDim ] = (uint)std::max( 1.0, std::min( std::ceil( tExt / tCell ), 4096.0 ) );
                    mH[ iDim ] = tExt / mN[ iDim ];
                    tNumCells *= mN[ iDim ];
                }

                // cell -> facets (CSR)
                mCellStart.assign( tNumCells + 1, 0 );
                std::vector< uint > tLoC( mNumFacets * 3 ), tHiC( mNumFacets * 3 );
                for ( uint iPass = 0; iPass < 2; ++iPass )
                {
                    std::vector< uint > tFill;
                    if ( iPass == 1 )
                    {
                        for ( uint iCell = 0; iCell < tNumCells; ++iCell )
                        {
                            mCellStart[ iCell + 1 ] += mCellStart[ iCell ];
                        }
                        mCellFacets.assign( mCellStart[ tNumCells ], 0 );
                        tFill.assign( mCellStart.begin(), mCellStart.end() - 1 );
                    }
                    for ( uint iFacet = 0; iFacet < mNumFacets; ++iFacet )
                    {
                        uint tC0[ 3 ] = { 0, 0, 0 };
                        uint tC1[ 3 ] = { 0, 0, 0 };
                        for ( uint iDim = 0; iDim < mDim; ++iDim )
                        {
                            tC0[ iDim ] = this->cell_index( tFLo[ iFacet * mDim + iDim ], iDim );
                            tC1[ iDim ] = this->cell_index( tFHi[ iFacet * mDim + iDim ], iDim );
                        }
                        for ( uint k = tC0[ 2 ]; k <= tC1[ 2 ]; ++k )
                        {
                            for ( uint j = tC0[ 1 ]; j <= tC1[ 1 ]; ++j )
                            {
                                for ( uint i = tC0[ 0 ]; i <= tC1[ 0 ]; ++i )
                                {
                                    uint tCellId = i + mN[ 0 ] * ( j + mN[ 1 ] * k );
                                    if ( iPass == 0 )
                                    {
                                        mCellStart[ tCellId + 1 ]++;
                                    }
                                    else
                                    {
                                        mCellFacets[ tFill[ tCellId ]++ ] = iFacet;
                                    }
                                }
                            }
                        }
                    }
                }
                mStamp.assign( mNumFacets, 0 );
            }

            // squared distance from aP to the nearest facet (exact)
            real
            nearest_distance_squared( const real* aP )
            {
                ++mQuery;
                if ( mQuery == 0 )    // stamp counter wrapped: reset
                {
                    std::fill( mStamp.begin(), mStamp.end(), 0 );
                    mQuery = 1;
                }

                int tC[ 3 ] = { 0, 0, 0 };
                for ( uint iDim = 0; iDim < mDim; ++iDim )
                {
                    tC[ iDim ] = (int)this->cell_index( aP[ iDim ], iDim );
                }

                real tBest = std::numeric_limits< real >::max();
                int  tMaxRing = 0;
                for ( uint iDim = 0; iDim < mDim; ++iDim )
                {
                    tMaxRing = std::max( tMaxRing, std::max( tC[ iDim ], (int)mN[ iDim ] - 1 - tC[ iDim ] ) );
                }

                for ( int tR = 0; tR <= tMaxRing; ++tR )
                {
                    int tLo[ 3 ] = { 0, 0, 0 };
                    int tHi[ 3 ] = { 0, 0, 0 };
                    for ( uint iDim = 0; iDim < mDim; ++iDim )
                    {
                        tLo[ iDim ] = std::max( tC[ iDim ] - tR, 0 );
                        tHi[ iDim ] = std::min( tC[ iDim ] + tR, (int)mN[ iDim ] - 1 );
                    }
                    for ( int k = tLo[ 2 ]; k <= tHi[ 2 ]; ++k )
                    {
                        for ( int j = tLo[ 1 ]; j <= tHi[ 1 ]; ++j )
                        {
                            for ( int i = tLo[ 0 ]; i <= tHi[ 0 ]; ++i )
                            {
                                // only the shell of ring tR (inner cells were visited before)
                                int tCheb = std::abs( i - tC[ 0 ] );
                                if ( mDim > 1 ) tCheb = std::max( tCheb, std::abs( j - tC[ 1 ] ) );
                                if ( mDim > 2 ) tCheb = std::max( tCheb, std::abs( k - tC[ 2 ] ) );
                                if ( tCheb != tR )
                                {
                                    continue;
                                }
                                uint tCellId = i + mN[ 0 ] * ( j + mN[ 1 ] * k );
                                for ( uint iF = mCellStart[ tCellId ]; iF < mCellStart[ tCellId + 1 ]; ++iF )
                                {
                                    uint tFacet = mCellFacets[ iF ];
                                    if ( mStamp[ tFacet ] == mQuery )
                                    {
                                        continue;
                                    }
                                    mStamp[ tFacet ] = mQuery;
                                    const real* tV = &mCoords[ tFacet * mStride ];
                                    real tDistSq = mUseTriangle
                                                         ? triangle_distance_squared( aP, tV, tV + 3, tV + 6 )
                                                         : segment_distance_squared( aP, tV, tV + mDim, mDim );
                                    tBest = std::min( tBest, tDistSq );
                                }
                            }
                        }
                    }

                    // lower bound on the distance to any cell outside ring tR
                    real tBound = std::numeric_limits< real >::max();
                    for ( uint iDim = 0; iDim < mDim; ++iDim )
                    {
                        if ( tC[ iDim ] + tR + 1 <= (int)mN[ iDim ] - 1 )
                        {
                            real tGap = mLo[ iDim ] + ( tC[ iDim ] + tR + 1 ) * mH[ iDim ] - aP[ iDim ];
                            tBound    = std::min( tBound, std::max( tGap, 0.0 ) );
                        }
                        if ( tC[ iDim ] - tR - 1 >= 0 )
                        {
                            real tGap = aP[ iDim ] - ( mLo[ iDim ] + ( tC[ iDim ] - tR ) * mH[ iDim ] );
                            tBound    = std::min( tBound, std::max( tGap, 0.0 ) );
                        }
                    }
                    if ( tBound == std::numeric_limits< real >::max() || tBest <= tBound * tBound )
                    {
                        break;
                    }
                }
                return tBest;
            }

          private:
            uint
            cell_index( real aX, uint aDim ) const
            {
                real tF = ( aX - mLo[ aDim ] ) / mH[ aDim ];
                if ( !( tF > 0.0 ) )
                {
                    return 0;
                }
                uint tI = (uint)tF;
                return tI >= mN[ aDim ] ? mN[ aDim ] - 1 : tI;
            }

            uint                mDim;
            uint                mNodesPerFacet;
            uint                mNumFacets;
            bool                mUseTriangle = false;
            uint                mStride      = 0;
            std::vector< real > mCoords;
            real                mLo[ 3 ];
            real                mH[ 3 ];
            uint                mN[ 3 ];
            std::vector< uint > mCellStart;
            std::vector< uint > mCellFacets;
            std::vector< uint > mStamp;
            uint                mQuery = 0;
        };
    }    // namespace

    //-------------------------------------------------------------------------------

    void
    SDF_From_Interface::compute(
            uint                         aNumNodes,
            const Matrix< DDRMat >&      aNodeCoords,
            const Matrix< DDRMat >&      aFacetNodeCoords,
            const Matrix< IndexMat >&    aFacetConn,
            const Vector< moris_index >& aNodeBulkPhase,
            moris_index                  aMaterialPhase,
            Matrix< DDRMat >&            aSDF )
    {
        // Get spatial dimension from node coordinates
        uint tSpatialDim = aNodeCoords.n_cols();

        // Resize output
        aSDF.set_size( aNumNodes, 1 );

        // no facets: the brute-force reference returns the largest real for every node
        if ( aFacetConn.n_rows() == 0 || tSpatialDim == 0 || tSpatialDim > 3 )
        {
            for ( uint iNode = 0; iNode < aNumNodes; ++iNode )
            {
                real tSign    = ( aNodeBulkPhase( iNode ) == aMaterialPhase ) ? -1.0 : 1.0;
                aSDF( iNode ) = tSign * std::numeric_limits< real >::max();
            }
            return;
        }

        // nearest-facet search on a bucket grid; exact, the same minimum as compute_distance_to_interface()
        // (the per-node brute-force loop over all facets: 80M point-facet tests with heap temporaries per X18 L-bracket
        // redistance, 2.8 s of every evaluation, runs/studies/lunch_and_learn_xiga/perf/x18v8/PROFILE.md)
        Facet_Grid tGrid( aFacetNodeCoords, aFacetConn, tSpatialDim );

        real tPoint[ 3 ] = { 0.0, 0.0, 0.0 };
        for ( uint iNode = 0; iNode < aNumNodes; ++iNode )
        {
            for ( uint iDim = 0; iDim < tSpatialDim; ++iDim )
            {
                tPoint[ iDim ] = aNodeCoords( iNode, iDim );
            }

            // unsigned distance to the nearest facet
            real tDistance = std::sqrt( tGrid.nearest_distance_squared( tPoint ) );

            // Determine sign: material phase = negative, void = positive
            real tSign = ( aNodeBulkPhase( iNode ) == aMaterialPhase ) ? -1.0 : 1.0;

            // Store signed distance
            aSDF( iNode ) = tSign * tDistance;
        }
    }

    //-------------------------------------------------------------------------------

    real
    SDF_From_Interface::compute_distance_to_interface(
            const Matrix< DDRMat >&   aPoint,
            const Matrix< DDRMat >&   aFacetNodeCoords,
            const Matrix< IndexMat >& aFacetConn,
            uint                      aSpatialDim )
    {
        real tMinDistance = std::numeric_limits< real >::max();

        uint tNumFacets      = aFacetConn.n_rows();
        uint tNodesPerFacet  = aFacetConn.n_cols();

        // Loop over all facets
        for ( uint iFacet = 0; iFacet < tNumFacets; ++iFacet )
        {
            real tDistance = 0.0;

            if ( aSpatialDim == 3 && tNodesPerFacet == 3 )
            {
                // 3D: Triangle facet
                moris_index tIdx0 = aFacetConn( iFacet, 0 );
                moris_index tIdx1 = aFacetConn( iFacet, 1 );
                moris_index tIdx2 = aFacetConn( iFacet, 2 );

                Matrix< DDRMat > tV0( 3, 1 );
                Matrix< DDRMat > tV1( 3, 1 );
                Matrix< DDRMat > tV2( 3, 1 );

                for ( uint iDim = 0; iDim < 3; ++iDim )
                {
                    tV0( iDim ) = aFacetNodeCoords( tIdx0, iDim );
                    tV1( iDim ) = aFacetNodeCoords( tIdx1, iDim );
                    tV2( iDim ) = aFacetNodeCoords( tIdx2, iDim );
                }

                tDistance = point_to_triangle_distance( aPoint, tV0, tV1, tV2 );
            }
            else if ( aSpatialDim == 2 && tNodesPerFacet == 2 )
            {
                // 2D: Line segment facet
                moris_index tIdx0 = aFacetConn( iFacet, 0 );
                moris_index tIdx1 = aFacetConn( iFacet, 1 );

                Matrix< DDRMat > tV0( 2, 1 );
                Matrix< DDRMat > tV1( 2, 1 );

                for ( uint iDim = 0; iDim < 2; ++iDim )
                {
                    tV0( iDim ) = aFacetNodeCoords( tIdx0, iDim );
                    tV1( iDim ) = aFacetNodeCoords( tIdx1, iDim );
                }

                tDistance = point_to_segment_distance( aPoint, tV0, tV1 );
            }
            else
            {
                // Fallback: treat as segment between first two nodes
                moris_index tIdx0 = aFacetConn( iFacet, 0 );
                moris_index tIdx1 = aFacetConn( iFacet, 1 );

                Matrix< DDRMat > tV0( aSpatialDim, 1 );
                Matrix< DDRMat > tV1( aSpatialDim, 1 );

                for ( uint iDim = 0; iDim < aSpatialDim; ++iDim )
                {
                    tV0( iDim ) = aFacetNodeCoords( tIdx0, iDim );
                    tV1( iDim ) = aFacetNodeCoords( tIdx1, iDim );
                }

                tDistance = point_to_segment_distance( aPoint, tV0, tV1 );
            }

            // Track minimum distance
            if ( tDistance < tMinDistance )
            {
                tMinDistance = tDistance;
            }
        }

        return tMinDistance;
    }

    //-------------------------------------------------------------------------------

    real
    SDF_From_Interface::point_to_triangle_distance(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aV0,
            const Matrix< DDRMat >& aV1,
            const Matrix< DDRMat >& aV2 )
    {
        // Edge vectors
        Matrix< DDRMat > tE0 = aV1 - aV0;    // Edge 0: V0 -> V1
        Matrix< DDRMat > tE1 = aV2 - aV0;    // Edge 1: V0 -> V2
        Matrix< DDRMat > tD  = aV0 - aPoint; // Vector from point to V0

        // Precompute dot products
        real tA = dot( tE0, tE0 );
        real tB = dot( tE0, tE1 );
        real tC = dot( tE1, tE1 );
        real tD0 = dot( tE0, tD );
        real tD1 = dot( tE1, tD );

        real tDet = tA * tC - tB * tB;
        real tS   = tB * tD1 - tC * tD0;
        real tT   = tB * tD0 - tA * tD1;

        // Determine which Voronoi region the point is in
        if ( tS + tT <= tDet )
        {
            if ( tS < 0.0 )
            {
                if ( tT < 0.0 )
                {
                    // Region 4: closest to V0
                    tS = clamp01( -tD0 / tA );
                    tT = 0.0;
                }
                else
                {
                    // Region 3: closest to edge V0-V2
                    tS = 0.0;
                    tT = clamp01( -tD1 / tC );
                }
            }
            else if ( tT < 0.0 )
            {
                // Region 5: closest to edge V0-V1
                tS = clamp01( -tD0 / tA );
                tT = 0.0;
            }
            else
            {
                // Region 0: inside triangle
                real tInvDet = 1.0 / tDet;
                tS *= tInvDet;
                tT *= tInvDet;
            }
        }
        else
        {
            if ( tS < 0.0 )
            {
                // Region 2: closest to edge V0-V2 or V1-V2
                real tTmp0 = tB + tD0;
                real tTmp1 = tC + tD1;
                if ( tTmp1 > tTmp0 )
                {
                    real tNumer = tTmp1 - tTmp0;
                    real tDenom = tA - 2.0 * tB + tC;
                    tS = clamp01( tNumer / tDenom );
                    tT = 1.0 - tS;
                }
                else
                {
                    tS = 0.0;
                    tT = clamp01( -tD1 / tC );
                }
            }
            else if ( tT < 0.0 )
            {
                // Region 6: closest to edge V0-V1 or V1-V2
                real tTmp0 = tB + tD1;
                real tTmp1 = tA + tD0;
                if ( tTmp1 > tTmp0 )
                {
                    real tNumer = tTmp1 - tTmp0;
                    real tDenom = tA - 2.0 * tB + tC;
                    tT = clamp01( tNumer / tDenom );
                    tS = 1.0 - tT;
                }
                else
                {
                    tS = clamp01( -tD0 / tA );
                    tT = 0.0;
                }
            }
            else
            {
                // Region 1: closest to edge V1-V2
                real tNumer = ( tC + tD1 ) - ( tB + tD0 );
                real tDenom = tA - 2.0 * tB + tC;
                tS          = clamp01( tNumer / tDenom );
                tT          = 1.0 - tS;
            }
        }

        // Compute closest point: P_closest = V0 + s*E0 + t*E1
        Matrix< DDRMat > tClosest = aV0 + tS * tE0 + tT * tE1;

        // Return distance
        Matrix< DDRMat > tDiff = aPoint - tClosest;
        return norm( tDiff );
    }

    //-------------------------------------------------------------------------------

    real
    SDF_From_Interface::point_to_segment_distance(
            const Matrix< DDRMat >& aPoint,
            const Matrix< DDRMat >& aV0,
            const Matrix< DDRMat >& aV1 )
    {
        Matrix< DDRMat > tEdge = aV1 - aV0;
        Matrix< DDRMat > tD    = aPoint - aV0;

        real tEdgeLengthSq = dot( tEdge, tEdge );

        // Handle degenerate case (zero-length segment)
        if ( tEdgeLengthSq < 1.0e-14 )
        {
            return norm( tD );
        }

        // Project point onto line, clamp to segment
        real tT = clamp01( dot( tD, tEdge ) / tEdgeLengthSq );

        // Closest point on segment
        Matrix< DDRMat > tClosest = aV0 + tT * tEdge;

        // Return distance
        Matrix< DDRMat > tDiff = aPoint - tClosest;
        return norm( tDiff );
    }

    //-------------------------------------------------------------------------------

    real
    SDF_From_Interface::clamp01( real aValue )
    {
        if ( aValue < 0.0 )
        {
            return 0.0;
        }
        if ( aValue > 1.0 )
        {
            return 1.0;
        }
        return aValue;
    }

    //-------------------------------------------------------------------------------

}    // namespace moris::sdf
