/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_SDF_From_Interface.cpp
 *
 * Unit tests for SDF_From_Interface class - distance computation from interface facets
 *
 */

#include <catch.hpp>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "fn_norm.hpp"
#include "cl_Vector.hpp"
#include "cl_SDF_From_Interface.hpp"

using namespace moris;
using namespace sdf;

TEST_CASE(
        "gen::sdf::SDF_From_Interface - Point to Segment",
        "[geomeng],[sdf],[SDF_From_Interface]" )
{
    real tEpsilon = 1e-10;

    //-------------------------------------------------------------------------------

    SECTION( "Point closest to segment midpoint" )
    {
        // Segment from (0,0) to (2,0)
        Matrix< DDRMat > tV0( 2, 1 );
        Matrix< DDRMat > tV1( 2, 1 );
        tV0( 0 ) = 0.0; tV0( 1 ) = 0.0;
        tV1( 0 ) = 2.0; tV1( 1 ) = 0.0;

        // Point at (1, 1) - directly above midpoint
        Matrix< DDRMat > tPoint( 2, 1 );
        tPoint( 0 ) = 1.0;
        tPoint( 1 ) = 1.0;

        real tDist = SDF_From_Interface::point_to_segment_distance( tPoint, tV0, tV1 );

        // Expected distance: 1.0 (straight up from midpoint)
        REQUIRE( std::abs( tDist - 1.0 ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point closest to segment endpoint V0" )
    {
        // Segment from (0,0) to (2,0)
        Matrix< DDRMat > tV0( 2, 1 );
        Matrix< DDRMat > tV1( 2, 1 );
        tV0( 0 ) = 0.0; tV0( 1 ) = 0.0;
        tV1( 0 ) = 2.0; tV1( 1 ) = 0.0;

        // Point at (-1, 0) - past V0
        Matrix< DDRMat > tPoint( 2, 1 );
        tPoint( 0 ) = -1.0;
        tPoint( 1 ) = 0.0;

        real tDist = SDF_From_Interface::point_to_segment_distance( tPoint, tV0, tV1 );

        // Expected distance: 1.0 (distance to V0)
        REQUIRE( std::abs( tDist - 1.0 ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point closest to segment endpoint V1" )
    {
        // Segment from (0,0) to (2,0)
        Matrix< DDRMat > tV0( 2, 1 );
        Matrix< DDRMat > tV1( 2, 1 );
        tV0( 0 ) = 0.0; tV0( 1 ) = 0.0;
        tV1( 0 ) = 2.0; tV1( 1 ) = 0.0;

        // Point at (3, 0) - past V1
        Matrix< DDRMat > tPoint( 2, 1 );
        tPoint( 0 ) = 3.0;
        tPoint( 1 ) = 0.0;

        real tDist = SDF_From_Interface::point_to_segment_distance( tPoint, tV0, tV1 );

        // Expected distance: 1.0 (distance to V1)
        REQUIRE( std::abs( tDist - 1.0 ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point on segment" )
    {
        // Segment from (0,0) to (2,0)
        Matrix< DDRMat > tV0( 2, 1 );
        Matrix< DDRMat > tV1( 2, 1 );
        tV0( 0 ) = 0.0; tV0( 1 ) = 0.0;
        tV1( 0 ) = 2.0; tV1( 1 ) = 0.0;

        // Point at (1, 0) - on segment
        Matrix< DDRMat > tPoint( 2, 1 );
        tPoint( 0 ) = 1.0;
        tPoint( 1 ) = 0.0;

        real tDist = SDF_From_Interface::point_to_segment_distance( tPoint, tV0, tV1 );

        // Expected distance: 0.0
        REQUIRE( tDist < tEpsilon );
    }
}

//-------------------------------------------------------------------------------

TEST_CASE(
        "gen::sdf::SDF_From_Interface - Point to Triangle",
        "[geomeng],[sdf],[SDF_From_Interface]" )
{
    real tEpsilon = 1e-10;

    // Simple right triangle in xy-plane at z=0
    // V0 = (0,0,0), V1 = (1,0,0), V2 = (0,1,0)
    Matrix< DDRMat > tV0( 3, 1 );
    Matrix< DDRMat > tV1( 3, 1 );
    Matrix< DDRMat > tV2( 3, 1 );

    tV0( 0 ) = 0.0; tV0( 1 ) = 0.0; tV0( 2 ) = 0.0;
    tV1( 0 ) = 1.0; tV1( 1 ) = 0.0; tV1( 2 ) = 0.0;
    tV2( 0 ) = 0.0; tV2( 1 ) = 1.0; tV2( 2 ) = 0.0;

    //-------------------------------------------------------------------------------

    SECTION( "Point directly above triangle interior" )
    {
        // Point at (0.25, 0.25, 1.0) - above centroid region
        Matrix< DDRMat > tPoint( 3, 1 );
        tPoint( 0 ) = 0.25;
        tPoint( 1 ) = 0.25;
        tPoint( 2 ) = 1.0;

        real tDist = SDF_From_Interface::point_to_triangle_distance( tPoint, tV0, tV1, tV2 );

        // Expected distance: 1.0 (straight down to plane)
        REQUIRE( std::abs( tDist - 1.0 ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point directly below triangle interior" )
    {
        // Point at (0.25, 0.25, -2.0)
        Matrix< DDRMat > tPoint( 3, 1 );
        tPoint( 0 ) = 0.25;
        tPoint( 1 ) = 0.25;
        tPoint( 2 ) = -2.0;

        real tDist = SDF_From_Interface::point_to_triangle_distance( tPoint, tV0, tV1, tV2 );

        // Expected distance: 2.0
        REQUIRE( std::abs( tDist - 2.0 ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point closest to vertex V0" )
    {
        // Point at (-1, -1, 0) - closest to origin
        Matrix< DDRMat > tPoint( 3, 1 );
        tPoint( 0 ) = -1.0;
        tPoint( 1 ) = -1.0;
        tPoint( 2 ) = 0.0;

        real tDist = SDF_From_Interface::point_to_triangle_distance( tPoint, tV0, tV1, tV2 );

        // Expected distance: sqrt(2) ≈ 1.414
        real tExpected = std::sqrt( 2.0 );
        REQUIRE( std::abs( tDist - tExpected ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point closest to edge V0-V1" )
    {
        // Point at (0.5, -1, 0) - below edge V0-V1
        Matrix< DDRMat > tPoint( 3, 1 );
        tPoint( 0 ) = 0.5;
        tPoint( 1 ) = -1.0;
        tPoint( 2 ) = 0.0;

        real tDist = SDF_From_Interface::point_to_triangle_distance( tPoint, tV0, tV1, tV2 );

        // Expected distance: 1.0 (perpendicular to edge)
        REQUIRE( std::abs( tDist - 1.0 ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point on triangle surface" )
    {
        // Point at (0.25, 0.25, 0) - on triangle
        Matrix< DDRMat > tPoint( 3, 1 );
        tPoint( 0 ) = 0.25;
        tPoint( 1 ) = 0.25;
        tPoint( 2 ) = 0.0;

        real tDist = SDF_From_Interface::point_to_triangle_distance( tPoint, tV0, tV1, tV2 );

        // Expected distance: 0.0
        REQUIRE( tDist < tEpsilon );
    }
}

//-------------------------------------------------------------------------------

TEST_CASE(
        "gen::sdf::SDF_From_Interface - Full SDF Computation 2D",
        "[geomeng],[sdf],[SDF_From_Interface]" )
{
    real tEpsilon = 1e-10;

    // Simple 2D test: unit square with vertical interface at x=0.5
    // Interface is a vertical line segment from (0.5, 0) to (0.5, 1)

    // Four corner nodes
    uint tNumNodes = 4;
    Matrix< DDRMat > tNodeCoords( 4, 2 );
    tNodeCoords( 0, 0 ) = 0.0; tNodeCoords( 0, 1 ) = 0.0;  // Node 0: (0,0)
    tNodeCoords( 1, 0 ) = 1.0; tNodeCoords( 1, 1 ) = 0.0;  // Node 1: (1,0)
    tNodeCoords( 2, 0 ) = 1.0; tNodeCoords( 2, 1 ) = 1.0;  // Node 2: (1,1)
    tNodeCoords( 3, 0 ) = 0.0; tNodeCoords( 3, 1 ) = 1.0;  // Node 3: (0,1)

    // Interface: single line segment
    Matrix< DDRMat > tFacetNodeCoords( 2, 2 );
    tFacetNodeCoords( 0, 0 ) = 0.5; tFacetNodeCoords( 0, 1 ) = 0.0;  // Facet node 0: (0.5, 0)
    tFacetNodeCoords( 1, 0 ) = 0.5; tFacetNodeCoords( 1, 1 ) = 1.0;  // Facet node 1: (0.5, 1)

    // Connectivity: one segment
    Matrix< IndexMat > tFacetConn( 1, 2 );
    tFacetConn( 0, 0 ) = 0;
    tFacetConn( 0, 1 ) = 1;

    // Bulk phases: left side (x < 0.5) is material (phase 0), right is void (phase 1)
    Vector< moris_index > tNodeBulkPhase( 4 );
    tNodeBulkPhase( 0 ) = 0;  // Node at x=0 is material
    tNodeBulkPhase( 1 ) = 1;  // Node at x=1 is void
    tNodeBulkPhase( 2 ) = 1;  // Node at x=1 is void
    tNodeBulkPhase( 3 ) = 0;  // Node at x=0 is material

    moris_index tMaterialPhase = 0;

    // Compute SDF
    Matrix< DDRMat > tSDF;
    SDF_From_Interface::compute(
            tNumNodes,
            tNodeCoords,
            tFacetNodeCoords,
            tFacetConn,
            tNodeBulkPhase,
            tMaterialPhase,
            tSDF );

    //-------------------------------------------------------------------------------

    SECTION( "Verify SDF values at corners" )
    {
        // Node 0 at (0,0): material, distance = 0.5, SDF = -0.5
        REQUIRE( std::abs( tSDF( 0 ) - ( -0.5 ) ) < tEpsilon );

        // Node 1 at (1,0): void, distance = 0.5, SDF = +0.5
        REQUIRE( std::abs( tSDF( 1 ) - 0.5 ) < tEpsilon );

        // Node 2 at (1,1): void, distance = 0.5, SDF = +0.5
        REQUIRE( std::abs( tSDF( 2 ) - 0.5 ) < tEpsilon );

        // Node 3 at (0,1): material, distance = 0.5, SDF = -0.5
        REQUIRE( std::abs( tSDF( 3 ) - ( -0.5 ) ) < tEpsilon );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Verify sign convention" )
    {
        // Material nodes have negative SDF
        REQUIRE( tSDF( 0 ) < 0.0 );
        REQUIRE( tSDF( 3 ) < 0.0 );

        // Void nodes have positive SDF
        REQUIRE( tSDF( 1 ) > 0.0 );
        REQUIRE( tSDF( 2 ) > 0.0 );
    }
}

//-------------------------------------------------------------------------------

TEST_CASE(
        "gen::sdf::SDF_From_Interface - Circle Interface",
        "[geomeng],[sdf],[SDF_From_Interface]" )
{
    // Approximate a circle of radius R=0.5 centered at (0.5, 0.5, 0)
    // using a regular polygon with N segments

    real tRadius = 0.5;
    real tCenterX = 0.5;
    real tCenterY = 0.5;
    uint tNumSegments = 32;

    // Create facet vertices on circle
    Matrix< DDRMat > tFacetNodeCoords( tNumSegments, 2 );
    for ( uint i = 0; i < tNumSegments; ++i )
    {
        real tTheta = 2.0 * M_PI * (real)i / (real)tNumSegments;
        tFacetNodeCoords( i, 0 ) = tCenterX + tRadius * std::cos( tTheta );
        tFacetNodeCoords( i, 1 ) = tCenterY + tRadius * std::sin( tTheta );
    }

    // Connectivity: closed polygon
    Matrix< IndexMat > tFacetConn( tNumSegments, 2 );
    for ( uint i = 0; i < tNumSegments; ++i )
    {
        tFacetConn( i, 0 ) = i;
        tFacetConn( i, 1 ) = ( i + 1 ) % tNumSegments;
    }

    // Test points: center, on circle, and outside
    uint tNumNodes = 4;
    Matrix< DDRMat > tNodeCoords( 4, 2 );
    tNodeCoords( 0, 0 ) = 0.5; tNodeCoords( 0, 1 ) = 0.5;   // Center
    tNodeCoords( 1, 0 ) = 0.5; tNodeCoords( 1, 1 ) = 0.0;   // On boundary (south)
    tNodeCoords( 2, 0 ) = 0.5; tNodeCoords( 2, 1 ) = 1.5;   // Outside (north)
    tNodeCoords( 3, 0 ) = 0.0; tNodeCoords( 3, 1 ) = 0.0;   // Outside (corner)

    Vector< moris_index > tNodeBulkPhase( 4 );
    tNodeBulkPhase( 0 ) = 0;  // Inside
    tNodeBulkPhase( 1 ) = 0;  // On boundary (treat as inside)
    tNodeBulkPhase( 2 ) = 1;  // Outside
    tNodeBulkPhase( 3 ) = 1;  // Outside

    Matrix< DDRMat > tSDF;
    SDF_From_Interface::compute( tNumNodes, tNodeCoords, tFacetNodeCoords, tFacetConn, tNodeBulkPhase, 0, tSDF );

    //-------------------------------------------------------------------------------

    SECTION( "Center point SDF approximately -R" )
    {
        // Center should have SDF ≈ -R = -0.5 (negative because inside)
        // Polygon approximation introduces small error
        real tExpected = -tRadius;
        real tTolerance = 0.01;  // 1% tolerance for polygon approximation
        REQUIRE( std::abs( tSDF( 0 ) - tExpected ) < tTolerance );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Point on boundary has SDF near zero" )
    {
        // Point at (0.5, 0) is at distance 0 from circle
        real tTolerance = 0.02;  // Polygon approximation
        REQUIRE( std::abs( tSDF( 1 ) ) < tTolerance );
    }

    //-------------------------------------------------------------------------------

    SECTION( "Outside point has positive SDF" )
    {
        // Point at (0.5, 1.5) is distance 0.5 from boundary (outside)
        real tExpected = 0.5;  // Distance from (0.5, 1.5) to (0.5, 1.0)
        real tTolerance = 0.02;
        REQUIRE( std::abs( tSDF( 2 ) - tExpected ) < tTolerance );
    }
}
