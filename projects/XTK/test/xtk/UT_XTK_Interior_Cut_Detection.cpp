/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Interior_Cut_Detection.cpp
 *
 * Regression tests for MB-031 / MB-016: the cut-cell selection must agree with the vertices at which the
 * decomposition evaluates an analytic level-set field. An analytic field (no discretization, no multilinear
 * interpolation) is evaluated exactly at the regular-subdivision vertices (cell centre in 2D; cell and face
 * centres in 3D) and at the intersection vertices of earlier geometries, so a sub-cell feature can change sign
 * there while every background corner keeps one sign.
 *
 */

#include "catch.hpp"

#include <cmath>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_User_Defined_Field.hpp"

namespace moris::xtk
{
    namespace
    {
        // radius of the 2D void used by the area-continuity test (set before each decomposition)
        real gMB031Radius = 0.2;

        // distance to a point minus a radius (negative inside the disc / ball)
        real
        distance_minus_radius(
                const Matrix< DDRMat >& aPoint,
                const Vector< real >&   aCentre,
                real                    aRadius )
        {
            real tSquaredDistance = 0.0;
            for ( uint iDim = 0; iDim < aPoint.numel(); iDim++ )
            {
                tSquaredDistance += ( aPoint( iDim ) - aCentre( iDim ) ) * ( aPoint( iDim ) - aCentre( iDim ) );
            }
            return std::sqrt( tSquaredDistance ) - aRadius;
        }

        // 2D: disc centred in background cell [0,1]x[0,1]
        real
        mb031_disc( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return distance_minus_radius( aPoint, { 0.5, 0.5 }, gMB031Radius );
        }

        // 3D: ball on the centre of the face x = 1 shared by the two cells, plus a ball on the corner (0,0,0)
        real
        mb016_balls_3d( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return std::min(
                    distance_minus_radius( aPoint, { 1.0, 0.5, 0.5 }, 0.3 ),
                    distance_minus_radius( aPoint, { 0.0, 0.0, 0.0 }, 0.3 ) );
        }

        // 2D, geometry 0: horizontal line y = 0.4 (cuts both cells and the shared edge x = 1 at (1, 0.4))
        real
        mb016_line( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return aPoint( 1 ) - 0.4;
        }

        // 2D, geometry 1: disc on the geometry-0 cut vertex (1, 0.4) plus a disc on the corner (0, 0)
        real
        mb016_discs_2d( const Matrix< DDRMat >& aPoint, const Vector< real >& aParameters )
        {
            return std::min(
                    distance_minus_radius( aPoint, { 1.0, 0.4 }, 0.1 ),
                    distance_minus_radius( aPoint, { 0.0, 0.0 }, 0.1 ) );
        }

        // finalized interpolation mesh of a linear HMR mesh
        hmr::Interpolation_Mesh_HMR*
        create_mesh( hmr::HMR& aHMR )
        {
            aHMR.finalize();
            return aHMR.create_interpolation_mesh( 0 );
        }

        hmr::Parameters
        mesh_parameters(
                const Vector< uint >& aNumElements,
                const Vector< real >& aDimensions )
        {
            hmr::Parameters tParameters;
            tParameters.set_number_of_elements_per_dimension( aNumElements );
            tParameters.set_domain_dimensions( aDimensions );
            tParameters.set_domain_offset( Vector< real >( aDimensions.size(), 0.0 ) );
            tParameters.set_bspline_truncation( true );
            tParameters.set_lagrange_orders( { 1 } );
            tParameters.set_lagrange_patterns( { 0 } );
            tParameters.set_bspline_orders( { 1 } );
            tParameters.set_bspline_patterns( { 0 } );
            tParameters.set_create_side_sets( true );
            return tParameters;
        }

        // index of the background cell whose vertex centroid lies closest to the given point
        moris_index
        find_cell( mtk::Interpolation_Mesh* aMesh, const Vector< real >& aPoint )
        {
            moris_index tBest     = 0;
            real        tBestDist = MORIS_REAL_MAX;
            for ( uint iCell = 0; iCell < aMesh->get_num_elems(); iCell++ )
            {
                Matrix< DDRMat > tCoords = aMesh->get_mtk_cell( iCell ).get_vertex_coords();
                real             tDist   = 0.0;
                for ( uint iDim = 0; iDim < tCoords.n_cols(); iDim++ )
                {
                    real tMean = 0.0;
                    for ( uint iV = 0; iV < tCoords.n_rows(); iV++ )
                    {
                        tMean += tCoords( iV, iDim ) / tCoords.n_rows();
                    }
                    tDist += ( tMean - aPoint( iDim ) ) * ( tMean - aPoint( iDim ) );
                }
                if ( tDist < tBestDist )
                {
                    tBestDist = tDist;
                    tBest     = iCell;
                }
            }
            return tBest;
        }

        // measure of all integration cells of a bulk phase
        real
        bulk_phase_measure( Cut_Integration_Mesh* aCutMesh, moris_index aBulkPhase )
        {
            real tMeasure = 0.0;
            for ( uint iSP = 0; iSP < aCutMesh->get_num_subphases(); iSP++ )
            {
                if ( aCutMesh->get_subphase_bulk_phase( iSP ) == aBulkPhase )
                {
                    for ( mtk::Cell* iCell : aCutMesh->get_subphase_ig_cells( iSP )->mIgCellGroup )
                    {
                        tMeasure += iCell->compute_cell_measure();
                    }
                }
            }
            return tMeasure;
        }

        // decompose the 2x2 mesh on [0,2]^2 with the disc of radius aRadius; returns the void (phase 0) area
        real
        void_area_2d( real aRadius, uint& aNumSubphasesInDiscCell )
        {
            gMB031Radius = aRadius;

            hmr::Parameters              tParameters = mesh_parameters( { 2, 2 }, { 2.0, 2.0 } );
            hmr::HMR                     tHMR( tParameters );
            hmr::Interpolation_Mesh_HMR* tInterpMesh = create_mesh( tHMR );

            auto tField = std::make_shared< gen::User_Defined_Field >( &mb031_disc );

            gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = { std::make_shared< gen::Level_Set_Geometry >( tField ) };
            gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            Model tXTKModel( 2, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;
            tXTKModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 } );

            Cut_Integration_Mesh* tCutMesh = tXTKModel.get_cut_integration_mesh();

            aNumSubphasesInDiscCell = tCutMesh->get_parent_cell_subphases( find_cell( tInterpMesh, { 0.5, 0.5 } ) ).size();

            real tVoidArea = bulk_phase_measure( tCutMesh, 0 );

            // the two phases tile the domain
            CHECK( tVoidArea + bulk_phase_measure( tCutMesh, 1 ) == Approx( 4.0 ).epsilon( 1e-10 ) );

            delete tInterpMesh;
            return tVoidArea;
        }
    }    // namespace

    TEST_CASE( "MB-031 2D void inside one background cell, all corners positive", "[XTK][MB-031][MB-031_2D]" )
    {
        if ( par_size() == 1 )
        {
            // disc of radius 0.2 centred in cell [0,1]^2: corners at distance 0.707 (phi > 0), cell centre phi = -0.2
            uint tNumSubphases = 0;
            real tVoidArea     = void_area_2d( 0.2, tNumSubphases );

            // the cell is decomposed (solid + void subphase)
            CHECK( tNumSubphases == 2 );

            // the void is the piecewise-linear interpolant on the 4-triangle subdivision: a diamond with
            // half-diagonals 0.2, area 2 * 0.2^2
            CHECK( tVoidArea == Approx( 0.08 ).epsilon( 1e-10 ) );
        }
    }

    TEST_CASE( "MB-031 2D void area is continuous when the corners change sign", "[XTK][MB-031][MB-031_2D_JUMP]" )
    {
        if ( par_size() == 1 )
        {
            // the four corners of cell [0,1]^2 lie on the circle of radius sqrt(0.5) around its centre
            real tCornerRadius = std::sqrt( 0.5 );
            uint tNumSubphases = 0;

            real tAreaBelow = void_area_2d( tCornerRadius - 1.0e-6, tNumSubphases );
            real tAreaAbove = void_area_2d( tCornerRadius + 1.0e-6, tNumSubphases );

            // before MB-031: the void is invisible below (area 0) and appears whole above (area ~1)
            CHECK( tAreaBelow > 0.99 );
            CHECK( std::abs( tAreaAbove - tAreaBelow ) < 1.0e-4 );
        }
    }

    TEST_CASE( "MB-016 3D void on a face centre shared by a cut and an uncut cell", "[XTK][MB-031][MB-016][MB-016_3D]" )
    {
        if ( par_size() == 1 )
        {
            // two hex cells [0,1]^3 and [1,2]x[0,1]^2. A ball of radius 0.3 on the shared face centre (1, 0.5, 0.5):
            // every corner and both cell centres are positive, the face centre is negative. A second ball on the
            // corner (0,0,0) makes the left cell intersected by its corners, the right one not.
            hmr::Parameters              tParameters = mesh_parameters( { 2, 1, 1 }, { 2.0, 1.0, 1.0 } );
            hmr::HMR                     tHMR( tParameters );
            hmr::Interpolation_Mesh_HMR* tInterpMesh = create_mesh( tHMR );

            auto tField = std::make_shared< gen::User_Defined_Field >( &mb016_balls_3d );

            gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = { std::make_shared< gen::Level_Set_Geometry >( tField ) };
            gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            Model tXTKModel( 3, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;
            tXTKModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 } );

            Cut_Integration_Mesh* tCutMesh = tXTKModel.get_cut_integration_mesh();

            // both cells see the face-centre void (left: solid + corner void + face void; right: solid + face void)
            CHECK( tCutMesh->get_parent_cell_subphases( find_cell( tInterpMesh, { 0.5, 0.5, 0.5 } ) ).size() == 3 );
            CHECK( tCutMesh->get_parent_cell_subphases( find_cell( tInterpMesh, { 1.5, 0.5, 0.5 } ) ).size() == 2 );

            // the two phases tile the domain
            CHECK( bulk_phase_measure( tCutMesh, 0 ) + bulk_phase_measure( tCutMesh, 1 ) == Approx( 2.0 ).epsilon( 1e-10 ) );

            delete tInterpMesh;
        }
    }

    TEST_CASE( "MB-016 2D second analytic geometry changes sign at a shared cut vertex", "[XTK][MB-031][MB-016][MB-016_2D]" )
    {
        if ( par_size() == 1 )
        {
            // two quad cells [0,1]^2 and [1,2]x[0,1]. Geometry 0 (y = 0.4) puts a vertex at (1, 0.4) on the shared edge.
            // Geometry 1 is negative there (disc of radius 0.1) and at the corner (0,0), positive at every other corner
            // and at both cell centres: only the left cell is intersected by its corners.
            hmr::Parameters              tParameters = mesh_parameters( { 2, 1 }, { 2.0, 1.0 } );
            hmr::HMR                     tHMR( tParameters );
            hmr::Interpolation_Mesh_HMR* tInterpMesh = create_mesh( tHMR );

            gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = {
                std::make_shared< gen::Level_Set_Geometry >( std::make_shared< gen::User_Defined_Field >( &mb016_line ) ),
                std::make_shared< gen::Level_Set_Geometry >( std::make_shared< gen::User_Defined_Field >( &mb016_discs_2d ) )
            };
            gen::Geometry_Engine tGeometryEngine( tInterpMesh, tGeometryEngineParameters );

            Model tXTKModel( 2, tInterpMesh, &tGeometryEngine );
            tXTKModel.mVerbose = false;
            tXTKModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 } );

            Cut_Integration_Mesh* tCutMesh = tXTKModel.get_cut_integration_mesh();

            // the right cell is cut by both geometries: (g0-, g1+), (g0+, g1+), (g0-, g1-), (g0+, g1-)
            CHECK( tCutMesh->get_parent_cell_subphases( find_cell( tInterpMesh, { 1.5, 0.5 } ) ).size() == 4 );

            // the four phases tile the domain
            real tTotal = 0.0;
            for ( moris_index iPhase = 0; iPhase < 4; iPhase++ )
            {
                tTotal += bulk_phase_measure( tCutMesh, iPhase );
            }
            CHECK( tTotal == Approx( 2.0 ).epsilon( 1e-10 ) );

            delete tInterpMesh;
        }
    }

}    // namespace moris::xtk
