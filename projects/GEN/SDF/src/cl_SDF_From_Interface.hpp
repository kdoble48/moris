/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_SDF_From_Interface.hpp
 *
 */

#pragma once

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::sdf
{
    //-------------------------------------------------------------------------------

    /**
     * @brief Computes a signed distance field from interface facet geometry.
     *
     * This class provides static methods to compute exact signed distances
     * from a triangulated interface (in 3D) or segmented interface (in 2D).
     * It is used for SDF reinitialization during level-set topology optimization
     * to eliminate floating artifacts and maintain the signed distance property.
     *
     * Usage:
     * @code
     * Matrix<DDRMat> tSDF;
     * SDF_From_Interface::compute(
     *     tMesh,
     *     tFacetNodeCoords,
     *     tFacetConnectivity,
     *     tNodeBulkPhases,
     *     tMaterialPhaseIndex,
     *     tSDF);
     * @endcode
     */
    class SDF_From_Interface
    {
      public:
        //-------------------------------------------------------------------------------

        /**
         * @brief Compute SDF at all mesh nodes from interface facet geometry.
         *
         * @param[in]  aNumNodes       Number of nodes in the mesh
         * @param[in]  aNodeCoords     Node coordinates (NumNodes x SpatialDim)
         * @param[in]  aFacetNodeCoords Coordinates of facet vertices (NumFacetNodes x SpatialDim)
         * @param[in]  aFacetConn      Facet connectivity (NumFacets x NodesPerFacet)
         * @param[in]  aNodeBulkPhase  Bulk phase index for each mesh node
         * @param[in]  aMaterialPhase  Phase index for "inside" (negative SDF)
         * @param[out] aSDF            Computed SDF values at mesh nodes (NumNodes x 1)
         */
        static void compute(
                uint                       aNumNodes,
                const Matrix< DDRMat >&    aNodeCoords,
                const Matrix< DDRMat >&    aFacetNodeCoords,
                const Matrix< IndexMat >&  aFacetConn,
                const Vector< moris_index >& aNodeBulkPhase,
                moris_index                aMaterialPhase,
                Matrix< DDRMat >&          aSDF );

        //-------------------------------------------------------------------------------

        /**
         * @brief Compute signed distance from a point to the nearest interface facet.
         *
         * @param[in] aPoint           Point coordinates (SpatialDim x 1)
         * @param[in] aFacetNodeCoords Coordinates of facet vertices
         * @param[in] aFacetConn       Facet connectivity
         * @param[in] aSpatialDim      Spatial dimension (2 or 3)
         * @return Unsigned distance to nearest facet
         */
        static real compute_distance_to_interface(
                const Matrix< DDRMat >&   aPoint,
                const Matrix< DDRMat >&   aFacetNodeCoords,
                const Matrix< IndexMat >& aFacetConn,
                uint                      aSpatialDim );

        //-------------------------------------------------------------------------------

        /**
         * @brief Compute distance from a point to a triangle (3D).
         *
         * Uses barycentric coordinates to determine the closest point
         * on the triangle (vertex, edge, or face interior).
         *
         * @param[in] aPoint      Point coordinates (3 x 1)
         * @param[in] aV0, aV1, aV2 Triangle vertex coordinates (3 x 1 each)
         * @return Distance from point to triangle
         */
        static real point_to_triangle_distance(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aV0,
                const Matrix< DDRMat >& aV1,
                const Matrix< DDRMat >& aV2 );

        //-------------------------------------------------------------------------------

        /**
         * @brief Compute distance from a point to a line segment (2D/3D).
         *
         * @param[in] aPoint    Point coordinates
         * @param[in] aV0, aV1  Segment endpoint coordinates
         * @return Distance from point to segment
         */
        static real point_to_segment_distance(
                const Matrix< DDRMat >& aPoint,
                const Matrix< DDRMat >& aV0,
                const Matrix< DDRMat >& aV1 );

        //-------------------------------------------------------------------------------

      private:
        /**
         * @brief Clamp a value to [0, 1] range.
         */
        static real clamp01( real aValue );

        //-------------------------------------------------------------------------------
    };

    //-------------------------------------------------------------------------------

}    // namespace moris::sdf
