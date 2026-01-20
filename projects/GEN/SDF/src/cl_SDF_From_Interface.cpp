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

#include <cmath>
#include <limits>

namespace moris::sdf
{
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

        // For each mesh node, compute distance to interface
        for ( uint iNode = 0; iNode < aNumNodes; ++iNode )
        {
            // Extract node coordinates
            Matrix< DDRMat > tPoint( tSpatialDim, 1 );
            for ( uint iDim = 0; iDim < tSpatialDim; ++iDim )
            {
                tPoint( iDim ) = aNodeCoords( iNode, iDim );
            }

            // Compute unsigned distance to nearest facet
            real tDistance = compute_distance_to_interface(
                    tPoint,
                    aFacetNodeCoords,
                    aFacetConn,
                    tSpatialDim );

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
