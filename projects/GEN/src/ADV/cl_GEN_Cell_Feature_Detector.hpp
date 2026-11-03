/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Cell_Feature_Detector.hpp
 *
 * Sub-cell feature detection for level-set fields that are evaluated exactly at any point (analytic and user-defined
 * fields, MB-084). XTK cuts an edge only where its end values change sign and selects cells from their corners and
 * regular-subdivision vertices, so a feature that crosses an edge twice, or lies between those samples, is not
 * represented. This detector samples the field on an (m+1)^d lattice of one QUAD4/HEX8 cell (corners in exodus
 * order) and reports, compared with what XTK's cut can represent:
 *   - whether the cell contains a zero crossing at all (lattice samples of different region);
 *   - background edges whose samples cross the interface more often than the cut does (even crossing counts);
 *   - crossings that XTK's samples (corners + regular-subdivision vertices) do not see;
 *   - cells whose connected sign components on the lattice differ from those of the piecewise-linear template that
 *     the cut builds (4 triangles around the cell centre in 2D; 24 tets around the face and cell centres in 3D).
 * Each hidden feature is weighed by its geometric depth |phi| / |grad phi| (central difference) relative to the cell
 * size; features below the depth tolerance are counted, not flagged (near-tangent chords, sub-resolution slivers).
 *
 * Resolution limit: without a Lipschitz bound only features that contain a lattice sample are found (a feature must
 * contain a ball of radius about h sqrt(d) / (2 m)). With a user-supplied bound L >= |grad phi| the detector bisects
 * every lattice segment / sub-cell that the bound cannot certify (|phi(a)| + |phi(b)| > L |b - a|, or
 * |phi(centre)| > L * half-diagonal), down to a fixed depth; what remains uncertified there is counted.
 */

#pragma once

#include <functional>
#include <string>

#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"

namespace moris::gen
{
    /**
     * Settings of one detection (see detect_cell_features()).
     */
    struct Cell_Feature_Settings
    {
        uint mSamples        = 4;      // sub-intervals per cell direction (rounded up to an even number >= 2)
        uint mLagrangeOrder  = 1;      // 1: XTK samples corners and subdivision vertices; 2: all 3^d nodes
        real mLipschitz      = 0.0;    // bound on |grad phi| for certified gaps; 0 = no certificate
        real mDepthTolerance = 0.0;    // hidden features shallower than this times the cell size are not significant
        bool mCheckTopology  = true;   // compare the lattice sign components with the template's
        uint mMaxEdgeBisections = 10;  // bisection depth of an uncertified lattice segment (Lipschitz only)
        uint mMaxCellBisections = 0;   // bisection depth of an uncertified lattice sub-cell (0: 3 in 2D, 2 in 3D)
    };

    /**
     * Result of one detection.
     */
    struct Cell_Feature_Report
    {
        bool mContainsCrossing  = false;    // the lattice (and the certificate) find a zero crossing in the cell
        bool mXTKSeesCut        = false;    // XTK's own samples (corners + subdivision vertices) differ in region
        uint mHiddenEdgeMask    = 0;        // bit e: edge e (exodus order) has significant crossings the cut misses
        bool mUnseenCrossing    = false;    // significant crossing that XTK's samples do not see
        bool mTopologyMismatch  = false;    // significant difference in sign components between lattice and template
        uint mNumInsignificant  = 0;        // hidden features found below the depth tolerance
        uint mNumUncertified    = 0;        // lattice segments / sub-cells the Lipschitz bound could not certify
        real mMaxHiddenDepth    = 0.0;      // largest depth of a significant hidden feature, relative to the cell size
        real mCellMeasure       = 0.0;      // area / volume of the cell
        real mCellSize          = 0.0;      // measure^(1/d)
        real mSampledPositive   = 0.0;      // trapezoidal estimate of the positive-region measure from the lattice
        real mTemplatePositive  = 0.0;      // the same for the piecewise-linear template of XTK's samples

        uint num_hidden_edges() const;

        /**
         * The cell holds a significant feature that XTK's cut of this cell does not represent.
         */
        bool hides_feature() const
        {
            return mHiddenEdgeMask != 0 or mUnseenCrossing or mTopologyMismatch;
        }
    };

    // phi( x ) minus the isocontour threshold, at global coordinates x (1 x d)
    using Point_Field_Function = std::function< real( const Matrix< DDRMat >& ) >;

    /**
     * Detect sub-cell features of a level-set field in one QUAD4 / HEX8 cell.
     *
     * @param aPhi Field minus the isocontour threshold, evaluated exactly at a point
     * @param aCorners Cell corner coordinates (2^d x d, exodus order; further rows are ignored)
     * @param aInterfaceTolerance |phi| at or below this is on the interface (isocontour_tolerance)
     * @param aSettings Lattice, template, certificate and significance settings
     * @return Detection report
     */
    Cell_Feature_Report detect_cell_features(
            const Point_Field_Function&  aPhi,
            const Matrix< DDRMat >&      aCorners,
            real                         aInterfaceTolerance,
            const Cell_Feature_Settings& aSettings );

    /**
     * Corner pairs (exodus numbering) of the edges of a QUAD4 (4 edges) or HEX8 (12 edges), in the order of the
     * bits of Cell_Feature_Report::mHiddenEdgeMask.
     */
    const Vector< std::pair< uint, uint > >& cell_feature_edges( uint aNumDimensions );

    // bit flags of a detection, as passed through the WRK refinement interface
    enum Cell_Feature_Flag : uint
    {
        CELL_FEATURE_CROSSING   = 1,    // mContainsCrossing
        CELL_FEATURE_HIDDEN     = 2,    // hides_feature()
        CELL_FEATURE_XTK_CUT    = 4,    // mXTKSeesCut
        CELL_FEATURE_INSIGNIFICANT = 8  // mNumInsignificant > 0
    };

    /**
     * Encode a report as Cell_Feature_Flag bits.
     */
    uint cell_feature_flags( const Cell_Feature_Report& aReport );

    /**
     * Geometry-aware query point (level-set keys feature_query_points / feature_query_function): a location at which
     * the user declares the side of the level set the geometry has there, so the cut can be checked at that point and
     * the Lagrange mesh refined until the cut represents it (a wedge tip or a neck that the depth tolerance misses).
     */
    struct Feature_Query_Point
    {
        Matrix< DDRMat > mCoordinates;     // 1 x d
        sint             mPhase  = 0;      // +1 positive side (phi > isocontour_threshold), -1 negative side, 0 = sign of the exact field
        real             mRadius = 0.0;    // refinement neighbourhood; 0 = one cell of the coarsest Lagrange level
    };

    /**
     * Deck function of feature_query_function: fills one row per query point, x, y[, z], phase, radius.
     */
    typedef void ( *Feature_Query_Function )( Matrix< DDRMat >& aQueryPoints );

    /**
     * Parse feature_query_points: points separated by ';', each "x,y[,z][:phase[:radius]]" (blanks ignored).
     */
    Vector< Feature_Query_Point > parse_feature_query_points( const std::string& aString );

    /**
     * Query points from the rows of a matrix (feature_query_function): x, y[, z], phase, radius (d + 2 columns).
     */
    Vector< Feature_Query_Point > feature_query_points_from_rows( const Matrix< DDRMat >& aRows );

    /**
     * Side of the level set that XTK's cut of one QUAD4 / HEX8 cell represents at a point of the cell: the sub-simplex
     * of the regular subdivision (4 triangles around the cell centre in 2D; 24 tets around the face and cell centres in
     * 3D) that contains the point, split by the plane through its edge intersections (exact roots of the field when
     * aExactRoots, as use_analytic_intersection; otherwise the linear crossing of the end values). In the 3D case of
     * two vertices on each side the cut prism is split into three tets, as XTK does up to the choice of the diagonal.
     *
     * @param aPhi Field minus the isocontour threshold, evaluated exactly at a point
     * @param aCorners Cell corner coordinates (2^d x d, exodus order)
     * @param aInterfaceTolerance |phi| at or below this is on the interface
     * @param aPoint Point (1 x d) inside the cell
     * @param aExactRoots Edge intersections at the exact root (true) or the linear crossing (false)
     * @return +1 positive side, -1 negative side, 0 every template vertex on the interface
     */
    sint represented_region_at_point(
            const Point_Field_Function& aPhi,
            const Matrix< DDRMat >&     aCorners,
            real                        aInterfaceTolerance,
            const Matrix< DDRMat >&     aPoint,
            bool                        aExactRoots );

    /**
     * Whether a point lies in a simplex (TRI / TET; only the first d + 1 rows of aVertexCoordinates, the corners, count).
     */
    bool simplex_contains_point( const Matrix< DDRMat >& aVertexCoordinates, const Matrix< DDRMat >& aPoint );

    /**
     * Whether a point lies in the axis-aligned box of a cell's corners (tolerance relative to the cell size).
     */
    bool cell_contains_point( const Matrix< DDRMat >& aCorners, const Matrix< DDRMat >& aPoint );

    /**
     * Distance of a point to the axis-aligned box of a cell's corners (0 inside).
     */
    real cell_distance_to_point( const Matrix< DDRMat >& aCorners, const Matrix< DDRMat >& aPoint );

}    // namespace moris::gen
