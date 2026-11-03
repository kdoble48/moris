/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Performer.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"

// This can be removed once we have the refinement interface
namespace moris::gen
{
    enum class Geometric_Region : signed char;
}

namespace moris::wrk
{
    class Performer
    {
        public:

            /**
             * Constructor
             */
            Performer() = default;

            /**
             * Destructor
             */
            virtual ~Performer() = default;

            /**
             * Return the number of fields that can be used for refinement
             *
             * @return Number of fields for refinement
             */
            virtual uint get_num_refinement_fields() = 0;

            /**
             * Returns fields so that HMR can perform refinement based on the data from this performer
             *
             * @param aFieldIndex Index of the field
             * @param aNodeIndex Index of the node
             * @param aCoordinates Coordinates of the node
             */
            virtual gen::Geometric_Region get_geometric_region(
                    uint                  aFieldIndex,
                    uint                  aNodeIndex,
                    const Matrix<DDRMat>& aCoordinates) = 0;

            virtual const Vector< uint >& get_num_refinements(uint aFieldIndex ) = 0;

            virtual const Vector< uint >& get_refinement_mesh_indices(uint aFieldIndex ) = 0;

            /**
             * Gets the index of an HMR user-defined refinement function for the given field index
             *
             * @param aFieldIndex Index of the field
             * @param aRefinementIndex The current refinement step being performed
             * @return User-defined function index, or -1 to use default refinement
             */
            virtual sint get_refinement_function_index(
                    uint aFieldIndex,
                    uint aRefinementIndex) = 0;

            /**
             * Sub-cell feature refinement (MB-084): lattice sub-intervals per cell direction for this field, 0 = off.
             * With it, number_of_refinements is the most refinement passes; passes stop once no cell hides a feature.
             */
            virtual uint get_feature_refinement_samples( uint aFieldIndex )
            {
                return 0;
            }

            /**
             * Whether the field can change region inside a cell whose nodes share one region (exact analytic fields).
             * Such cells can be cut for XTK although HMR's node-based criterion does not flag them (D-205).
             */
            virtual bool can_change_region_inside_cell( uint aFieldIndex )
            {
                return false;
            }

            /**
             * Sub-cell feature detection in one Lagrange cell (gen::Cell_Feature_Flag bits: crossing, hidden feature,
             * cut as seen by XTK's samples, feature below the depth tolerance).
             *
             * @param aFieldIndex Field index
             * @param aCorners Corner coordinates in exodus order (2^d x d)
             * @param aLagrangeOrder Order of the Lagrange cell
             * @param aSamples Lattice sub-intervals per direction (2: XTK's own samples only)
             */
            virtual uint detect_cell_features(
                    uint                  aFieldIndex,
                    const Matrix<DDRMat>& aCorners,
                    uint                  aLagrangeOrder,
                    uint                  aSamples )
            {
                return 0;
            }

            /**
             * Geometry-aware query points of a field used by the refinement (feature_query_refine): the number of
             * points whose cut the refinement checks (0: off), and the coordinates (1 x d) and radius of one point
             * (radius 0: one cell of the coarsest Lagrange level).
             */
            virtual uint get_num_feature_query_refine_points( uint aFieldIndex )
            {
                return 0;
            }

            virtual void get_feature_query_point(
                    uint              aFieldIndex,
                    uint              aPointIndex,
                    Matrix< DDRMat >& aCoordinates,
                    real&             aRadius )
            {
            }

            /**
             * Side of the field (+1 / -1, 0 on the interface) expected at a query point, and the side that XTK's cut of
             * the Lagrange cell aCorners (containing the point) would represent there.
             */
            virtual void get_feature_query_regions(
                    uint                    aFieldIndex,
                    uint                    aPointIndex,
                    const Matrix< DDRMat >& aCorners,
                    sint&                   aExpected,
                    sint&                   aRepresented )
            {
                aExpected    = 0;
                aRepresented = 0;
            }
    };
}
