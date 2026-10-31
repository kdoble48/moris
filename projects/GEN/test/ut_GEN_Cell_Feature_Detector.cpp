/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Cell_Feature_Detector.cpp
 *
 * MB-084: sub-cell feature detection of an exactly evaluated level set in one QUAD4 / HEX8 cell. Positive = material.
 * Cases: a slot between node lines at many offsets (translated, off-grid), a slot on a node line, an enclosed ring,
 * hairline bridges on and between node lines, a same-sign edge crossed twice, a near-tangent chord, a second void the
 * cut misses (topology), the Lipschitz certificate, and 3D slab / ball cases.
 *
 */

#include <catch.hpp>

#include <cmath>

#include "cl_GEN_Cell_Feature_Detector.hpp"

namespace moris::gen
{
    namespace
    {
        // corners (exodus order) of the box [x0, x0 + h] x [y0, y0 + h] (x [z0, z0 + h])
        Matrix< DDRMat >
        box_corners( real aX0, real aY0, real aH, int aZ0Flag = -1, real aZ0 = 0.0 )
        {
            uint             tNumDims = aZ0Flag < 0 ? 2 : 3;
            Matrix< DDRMat > tCorners( 1u << tNumDims, tNumDims );
            for ( uint iCorner = 0; iCorner < ( 1u << tNumDims ); iCorner++ )
            {
                uint tInFace = iCorner % 4;
                tCorners( iCorner, 0 ) = aX0 + ( ( tInFace == 1 or tInFace == 2 ) ? aH : 0.0 );
                tCorners( iCorner, 1 ) = aY0 + ( ( tInFace == 2 or tInFace == 3 ) ? aH : 0.0 );
                if ( tNumDims == 3 )
                {
                    tCorners( iCorner, 2 ) = aZ0 + ( iCorner >= 4 ? aH : 0.0 );
                }
            }
            return tCorners;
        }

        Cell_Feature_Settings
        settings( uint aSamples, real aDepthTolerance = 0.02, real aLipschitz = 0.0 )
        {
            Cell_Feature_Settings tSettings;
            tSettings.mSamples        = aSamples;
            tSettings.mDepthTolerance = aDepthTolerance;
            tSettings.mLipschitz      = aLipschitz;
            return tSettings;
        }

        const real gTol = 1.0e-12;

        // edge bits (exodus order): 0 bottom (y = y0), 1 right (x = x0 + h), 2 top, 3 left (x = x0)
        const uint gVerticalEdges = ( 1u << 1 ) | ( 1u << 3 );
    }    // namespace

    TEST_CASE( "MB-084 detector: slot between node lines, every offset", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // horizontal void slot of width 0.3 (positive outside) inside cell [0,1]^2, at offsets 0.17 ... 0.83:
        // both vertical edges cross it twice with the same sign at their ends
        real tWidth = 0.3;
        for ( real tCentre = 0.17; tCentre < 0.84; tCentre += 0.0233 )
        {
            Point_Field_Function tPhi = [ & ]( const Matrix< DDRMat >& aX ) { return std::abs( aX( 1 ) - tCentre ) - 0.5 * tWidth; };

            Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 4 ) );

            CHECK( tReport.mContainsCrossing );
            CHECK( tReport.mHiddenEdgeMask == gVerticalEdges );
            CHECK( tReport.hides_feature() );
            // XTK sees the cell only when the slot holds the cell centre (y = 0.5)
            CHECK( tReport.mXTKSeesCut == ( std::abs( tCentre - 0.5 ) < 0.5 * tWidth ) );
            CHECK( tReport.mUnseenCrossing == not tReport.mXTKSeesCut );
            // depth of the deepest sample in the slot: at most half the slot width (the lattice need not hit the slot
            // axis), and above the significance tolerance
            CHECK( tReport.mMaxHiddenDepth <= 0.5 * tWidth + 1.0e-9 );
            CHECK( tReport.mMaxHiddenDepth >= 0.02 );
        }
    }

    TEST_CASE( "MB-084 detector: translated and scaled cell", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // the same slot in the cell [10,12]x[-4,-2] (h = 2), slot 0.6 wide at y = -3.46 (off-grid)
        Point_Field_Function tPhi = []( const Matrix< DDRMat >& aX ) { return std::abs( aX( 1 ) + 3.46 ) - 0.3; };

        Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 10.0, -4.0, 2.0 ), gTol, settings( 4 ) );

        CHECK( tReport.mCellMeasure == Approx( 4.0 ) );
        CHECK( tReport.mHiddenEdgeMask == gVerticalEdges );
        CHECK( not tReport.mXTKSeesCut );
        CHECK( tReport.mUnseenCrossing );
        // deepest sample y = -3.5: |phi| = 0.26 of the 0.3 half-width, 0.13 cell sizes
        CHECK( tReport.mMaxHiddenDepth == Approx( 0.13 ).epsilon( 1e-6 ) );
    }

    TEST_CASE( "MB-084 detector: slot on a node line is represented", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // slot centred on the bottom edge y = 0 of the cell: every vertical edge changes sign once
        Point_Field_Function tPhi = []( const Matrix< DDRMat >& aX ) { return std::abs( aX( 1 ) - 0.0 ) - 0.15; };

        Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8 ) );

        CHECK( tReport.mContainsCrossing );
        CHECK( tReport.mXTKSeesCut );
        CHECK( tReport.mHiddenEdgeMask == 0 );
        CHECK( not tReport.hides_feature() );
    }

    TEST_CASE( "MB-084 detector: enclosed ring", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // ring wall (material) of radii 0.35 ... 0.45 around the cell centre: corners and centre are void
        Point_Field_Function tPhi = []( const Matrix< DDRMat >& aX ) {
            real tR = std::hypot( aX( 0 ) - 0.5, aX( 1 ) - 0.5 );
            return std::min( 0.45 - tR, tR - 0.35 );
        };

        Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8 ) );

        CHECK( tReport.mContainsCrossing );
        CHECK( not tReport.mXTKSeesCut );
        CHECK( tReport.mUnseenCrossing );
        CHECK( tReport.mHiddenEdgeMask == 0 );    // the ring touches no edge
        CHECK( tReport.hides_feature() );
    }

    TEST_CASE( "MB-084 detector: same-sign edge crossed twice", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // void disc of radius 0.15 centred on the bottom edge at x = 0.6: the edge crosses it at 0.45 and 0.75, both
        // ends and every XTK sample are material
        Point_Field_Function tPhi = []( const Matrix< DDRMat >& aX ) { return std::hypot( aX( 0 ) - 0.6, aX( 1 ) ) - 0.15; };

        Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8 ) );

        CHECK( tReport.mHiddenEdgeMask == 1u );
        CHECK( not tReport.mXTKSeesCut );
        CHECK( tReport.mUnseenCrossing );
    }

    TEST_CASE( "MB-084 detector: hairline bridges are below the resolution", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        real tS = 1.0e-4;

        SECTION( "between node lines, on a lattice line" )
        {
            // vertical material bridge of width s at x = 0.5 in a void cell: both horizontal edges cross it twice
            Point_Field_Function tPhi = [ & ]( const Matrix< DDRMat >& aX ) { return 0.5 * tS - std::abs( aX( 0 ) - 0.5 ); };

            Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8 ) );
            CHECK( tReport.mContainsCrossing );
            CHECK( not tReport.hides_feature() );    // depth s / 2 is below the tolerance 0.02 h
            CHECK( tReport.mNumInsignificant > 0 );

            // without the depth filter the bridge is a hidden feature on both horizontal edges
            Cell_Feature_Report tUnfiltered = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8, 0.0 ) );
            CHECK( tUnfiltered.mHiddenEdgeMask == ( 1u | ( 1u << 2 ) ) );
            CHECK( tUnfiltered.hides_feature() );
        }

        SECTION( "on a node line" )
        {
            // the same bridge on the left edge x = 0: the cut sees it at the corners; no hidden crossing, no change of
            // topology (the node-centred hourglass of MB-045 is a template error, not a detection miss)
            Point_Field_Function tPhi = [ & ]( const Matrix< DDRMat >& aX ) { return 0.5 * tS - std::abs( aX( 0 ) ); };

            Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8, 0.0 ) );
            CHECK( tReport.mXTKSeesCut );
            CHECK( tReport.mHiddenEdgeMask == 0 );
            CHECK( not tReport.mTopologyMismatch );
        }
    }

    TEST_CASE( "MB-084 detector: near-tangent chord is not significant", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // void disc of radius 5 grazing the bottom edge with a sagitta of 1e-4 (chord 0.063 around x = 0.5)
        Point_Field_Function tPhi = []( const Matrix< DDRMat >& aX ) { return std::hypot( aX( 0 ) - 0.5, aX( 1 ) + 5.0 - 1.0e-4 ) - 5.0; };

        Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8 ) );

        CHECK( tReport.mContainsCrossing );
        CHECK( not tReport.hides_feature() );
        CHECK( tReport.mNumInsignificant > 0 );
    }

    TEST_CASE( "MB-084 detector: a second void the cut misses changes the topology", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // void disc at the cell centre (seen by the centre vertex) plus a void disc near the corner (0, 1) that touches no
        // edge and holds no XTK sample
        Point_Field_Function tPhi = []( const Matrix< DDRMat >& aX ) {
            return std::min( std::hypot( aX( 0 ) - 0.5, aX( 1 ) - 0.5 ) - 0.15, std::hypot( aX( 0 ) - 0.22, aX( 1 ) - 0.78 ) - 0.1 );
        };

        Cell_Feature_Report tReport = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8 ) );

        CHECK( tReport.mXTKSeesCut );
        CHECK( tReport.mHiddenEdgeMask == 0 );
        CHECK( tReport.mTopologyMismatch );
        CHECK( tReport.hides_feature() );

        // a straight interface: no hidden feature, template and samples agree on the positive measure
        Point_Field_Function tLine         = []( const Matrix< DDRMat >& aX ) { return aX( 1 ) - 0.3 - 0.2 * aX( 0 ); };
        Cell_Feature_Report  tLineReport   = detect_cell_features( tLine, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 8 ) );
        CHECK( not tLineReport.hides_feature() );
        CHECK( tLineReport.mSampledPositive == Approx( tLineReport.mTemplatePositive ).epsilon( 1e-12 ) );
    }

    TEST_CASE( "MB-084 detector: Lipschitz certificate finds a slot between lattice lines", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // void slot of width 0.08 at y = 0.62 between the lattice lines y = 0.5 and 0.75 of m = 4 (1-Lipschitz field)
        Point_Field_Function tPhi = []( const Matrix< DDRMat >& aX ) { return std::abs( aX( 1 ) - 0.62 ) - 0.04; };

        Cell_Feature_Report tPlain = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 4 ) );
        CHECK( not tPlain.mContainsCrossing );    // resolution limit of plain sampling

        Cell_Feature_Report tCertified = detect_cell_features( tPhi, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 4, 0.02, 1.0 ) );
        CHECK( tCertified.mContainsCrossing );
        CHECK( tCertified.mHiddenEdgeMask == gVerticalEdges );
        CHECK( tCertified.hides_feature() );

        // far from any interface the bound certifies the cell without bisection
        Point_Field_Function tFar       = []( const Matrix< DDRMat >& aX ) { return 3.0 + aX( 0 ); };
        Cell_Feature_Report  tFarReport = detect_cell_features( tFar, box_corners( 0.0, 0.0, 1.0 ), gTol, settings( 4, 0.02, 1.0 ) );
        CHECK( not tFarReport.mContainsCrossing );
        CHECK( tFarReport.mNumUncertified == 0 );
    }

    TEST_CASE( "MB-084 detector: 3D slab and enclosed ball", "[GEN][MB-084][MB-084_DETECTOR]" )
    {
        // void slab of width 0.3 at z = 0.5: the four z-edges (8 ... 11) cross it twice; the cell centre sees it
        Point_Field_Function tSlab   = []( const Matrix< DDRMat >& aX ) { return std::abs( aX( 2 ) - 0.5 ) - 0.15; };
        Cell_Feature_Report  tReport = detect_cell_features( tSlab, box_corners( 0.0, 0.0, 1.0, 1, 0.0 ), gTol, settings( 4 ) );
        CHECK( tReport.mHiddenEdgeMask == ( ( 1u << 8 ) | ( 1u << 9 ) | ( 1u << 10 ) | ( 1u << 11 ) ) );
        CHECK( tReport.mXTKSeesCut );
        CHECK( tReport.mCellMeasure == Approx( 1.0 ) );

        // void ball of radius 0.15 at (0.3, 0.3, 0.3): no edge, face centre or cell centre in it
        Point_Field_Function tBall = []( const Matrix< DDRMat >& aX ) {
            return std::sqrt( ( aX( 0 ) - 0.3 ) * ( aX( 0 ) - 0.3 ) + ( aX( 1 ) - 0.3 ) * ( aX( 1 ) - 0.3 ) + ( aX( 2 ) - 0.3 ) * ( aX( 2 ) - 0.3 ) ) - 0.15;
        };
        Cell_Feature_Report tBallReport = detect_cell_features( tBall, box_corners( 0.0, 0.0, 1.0, 1, 0.0 ), gTol, settings( 4 ) );
        CHECK( tBallReport.mContainsCrossing );
        CHECK( not tBallReport.mXTKSeesCut );
        CHECK( tBallReport.mUnseenCrossing );
        CHECK( tBallReport.hides_feature() );

        // a plane through the cell: represented
        Point_Field_Function tPlane       = []( const Matrix< DDRMat >& aX ) { return aX( 0 ) + 0.5 * aX( 1 ) + 0.25 * aX( 2 ) - 0.8; };
        Cell_Feature_Report  tPlaneReport = detect_cell_features( tPlane, box_corners( 0.0, 0.0, 1.0, 1, 0.0 ), gTol, settings( 4 ) );
        CHECK( tPlaneReport.mXTKSeesCut );
        CHECK( not tPlaneReport.hides_feature() );
        CHECK( tPlaneReport.mSampledPositive == Approx( tPlaneReport.mTemplatePositive ).epsilon( 1e-12 ) );
    }
}    // namespace moris::gen
