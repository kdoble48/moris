/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Intersection_Analytic.cpp
 *
 * Analytic intersections (level-set key use_analytic_intersection): edge intersections found as the bracketed root of
 * the analytic field itself along the parent edge, instead of the zero of the linear interpolant of the two end values.
 *
 */

#include "catch.hpp"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include <limits>
#include "fn_eye.hpp"

#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Node_Manager.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Line.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Intersection_Node_Analytic.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"

#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    // circle (2D) and sphere (3D) of the precision tests, off-center so that cuts fall at generic edge positions
    static const real tCenterX = 0.1;
    static const real tCenterY = -0.05;
    static const real tCenterZ = 0.07;
    static const real tRadius  = 0.63;

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Statistics of the edge intersections an XTK regular subdivision queries on every cell of a mesh: the cell edges
     * (4 in 2D, 12 in 3D) and the spokes from each corner to the cell center (a derived node).
     */
    struct Cut_Statistics
    {
        real           mMaxDistance      = 0.0;    // max distance of an intersection from the analytic circle / sphere
        uint           mNumIntersections = 0;      // number of intersected edges
        Vector< real > mLocalCoordinates;          // all intersection local coordinates, in query order
    };

    //--------------------------------------------------------------------------------------------------------------

    static mtk::Interpolation_Mesh*
    create_cube_mesh( uint aNumElements )
    {
        Module_Parameter_Lists tParameters( Module_Type::HMR );

        tParameters.set( "number_of_elements_per_dimension", aNumElements, aNumElements, aNumElements );
        tParameters.set( "domain_dimensions", 2.0, 2.0, 2.0 );
        tParameters.set( "domain_offset", -1.0, -1.0, -1.0 );
        tParameters.set( "lagrange_output_meshes", "0" );
        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "bspline_pattern", "0" );
        tParameters.set( "lagrange_to_bspline", "0" );
        tParameters.set( "pattern_initial_refinement", 0 );
        tParameters.set( "refinement_buffer", 1 );
        tParameters.set( "staircase_buffer", 1 );
        tParameters.set( "severity_level", 2 );

        hmr::HMR tHMR( tParameters );
        tHMR.perform_initial_refinement();
        tHMR.finalize();

        return tHMR.create_interpolation_mesh( 0 );
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Circle (2D) or sphere (3D) parameter list.
     *
     * @param aNumDims 2 or 3
     * @param aAnalytic use_analytic_intersection (-1: leave the default)
     * @param aWithADVs Center and radius are ADVs (lower = value = upper)
     */
    static Submodule_Parameter_Lists
    create_round_parameters( uint aNumDims, sint aAnalytic, bool aWithADVs, uint aMaxIterations = 0 )
    {
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( aNumDims == 2 ? gen::Field_Type::CIRCLE : gen::Field_Type::SPHERE ) );
        if ( aWithADVs )
        {
            tFieldParameterLists.set( "center_x", tCenterX, tCenterX, tCenterX );
            tFieldParameterLists.set( "center_y", tCenterY, tCenterY, tCenterY );
            if ( aNumDims == 3 )
            {
                tFieldParameterLists.set( "center_z", tCenterZ, tCenterZ, tCenterZ );
            }
            tFieldParameterLists.set( "radius", tRadius, tRadius, tRadius );
        }
        else
        {
            tFieldParameterLists.set( "center_x", tCenterX );
            tFieldParameterLists.set( "center_y", tCenterY );
            if ( aNumDims == 3 )
            {
                tFieldParameterLists.set( "center_z", tCenterZ );
            }
            tFieldParameterLists.set( "radius", tRadius );
        }
        tFieldParameterLists.set( "discretization_mesh_index", -2 );
        if ( aMaxIterations > 0 )
        {
            tFieldParameterLists.set( "analytic_intersection_max_iterations", aMaxIterations );
        }
        if ( aAnalytic >= 0 )
        {
            tFieldParameterLists.set( "use_analytic_intersection", aAnalytic == 1 );
        }
        return tFieldParameterLists;
    }

    //--------------------------------------------------------------------------------------------------------------

    static real
    distance_to_round( const Matrix< DDRMat >& aCoordinates )
    {
        real tSquared = std::pow( aCoordinates( 0 ) - tCenterX, 2 ) + std::pow( aCoordinates( 1 ) - tCenterY, 2 );
        if ( aCoordinates.numel() == 3 )
        {
            tSquared += std::pow( aCoordinates( 2 ) - tCenterZ, 2 );
        }
        return std::abs( std::sqrt( tSquared ) - tRadius );
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Edges an XTK regular subdivision queries on one linear cell: { first, second } as cell-local indices, where index
     * aNumCorners is the cell center. Cell edges are the corner pairs that differ in exactly one parametric coordinate.
     */
    static Vector< Vector< uint > >
    get_cell_edges( const Vector< Matrix< DDRMat > >& aCornerXi )
    {
        uint                     tNumCorners = aCornerXi.size();
        Vector< Vector< uint > > tEdges;
        for ( uint iFirst = 0; iFirst < tNumCorners; iFirst++ )
        {
            for ( uint iSecond = iFirst + 1; iSecond < tNumCorners; iSecond++ )
            {
                uint tNumDifferent = 0;
                for ( uint iDim = 0; iDim < aCornerXi( iFirst ).numel(); iDim++ )
                {
                    tNumDifferent += aCornerXi( iFirst )( iDim ) != aCornerXi( iSecond )( iDim );
                }
                if ( tNumDifferent == 1 )
                {
                    tEdges.push_back( { iFirst, iSecond } );
                }
            }
        }
        for ( uint iCorner = 0; iCorner < tNumCorners; iCorner++ )
        {
            tEdges.push_back( { iCorner, tNumCorners } );
        }
        return tEdges;
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * One linear background cell: node indices, parametric corner coordinates (from the corner positions relative to
     * the cell center, so independent of the mesh's node order), geometry type.
     */
    struct Cell_Data
    {
        Matrix< IndexMat >         mSignedNodeIndices;
        Matrix< DDUMat >           mNodeIndices;
        Vector< Matrix< DDRMat > > mCornerXi;
        mtk::Geometry_Type         mGeometryType;
    };

    static Cell_Data
    get_cell_data( mtk::Interpolation_Mesh* aMesh, uint aCell )
    {
        Cell_Data tData;
        tData.mSignedNodeIndices = aMesh->get_nodes_connected_to_element_loc_inds( aCell );
        uint tNumCorners         = tData.mSignedNodeIndices.numel();
        uint tNumDims            = aMesh->get_spatial_dim();
        tData.mGeometryType      = tNumDims == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;
        tData.mNodeIndices.set_size( tNumCorners, 1 );

        Matrix< DDRMat > tCenter( 1, tNumDims, 0.0 );
        for ( uint iNode = 0; iNode < tNumCorners; iNode++ )
        {
            tData.mNodeIndices( iNode ) = tData.mSignedNodeIndices( iNode );
            tCenter                     = tCenter + aMesh->get_node_coordinate( tData.mNodeIndices( iNode ) ) / real( tNumCorners );
        }
        tData.mCornerXi.resize( tNumCorners );
        for ( uint iNode = 0; iNode < tNumCorners; iNode++ )
        {
            Matrix< DDRMat > tX = aMesh->get_node_coordinate( tData.mNodeIndices( iNode ) );
            tData.mCornerXi( iNode ).set_size( 1, tNumDims );
            for ( uint iDim = 0; iDim < tNumDims; iDim++ )
            {
                tData.mCornerXi( iNode )( iDim ) = tX( iDim ) > tCenter( iDim ) ? 1.0 : -1.0;
            }
        }
        return tData;
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Queues the regular-subdivision edges of every cell and collects the intersections. With aAdmitAndCheckFD, every
     * intersection is admitted and its analytic d(coordinates)/d(ADV) is compared to central finite differences of the
     * re-queued intersection coordinates; the maximum relative error is returned in aMaxFDError.
     */
    static Cut_Statistics
    cut_all_cells(
            mtk::Interpolation_Mesh* aMesh,
            Geometry_Engine_Test&    aGeometryEngine,
            bool                     aAdmitAndCheckFD = false,
            real*                    aMaxFDError      = nullptr,
            uint*                    aNumFDChecked    = nullptr )
    {
        Cut_Statistics tStatistics;
        Node_Manager&  tNodeManager = aGeometryEngine.get_node_manager();
        uint           tNumDims     = aMesh->get_spatial_dim();

        for ( uint iCell = 0; iCell < aMesh->get_num_elems(); iCell++ )
        {
            Cell_Data                tCell  = get_cell_data( aMesh, iCell );
            Vector< Vector< uint > > tEdges = get_cell_edges( tCell.mCornerXi );
            uint                     tNumCorners = tCell.mCornerXi.size();

            // cell center, a derived node of this cell
            uint tCenterNodeIndex = tNodeManager.get_total_number_of_nodes();
            aGeometryEngine.create_new_derived_nodes(
                    { tCell.mSignedNodeIndices }, { Matrix< DDRMat >( 1, tNumDims, 0.0 ) }, tCell.mGeometryType, mtk::Interpolation_Order::LINEAR );

            for ( const Vector< uint >& tEdge : tEdges )
            {
                uint             tFirstNode  = tCell.mNodeIndices( tEdge( 0 ) );
                uint             tSecondNode = tEdge( 1 ) < tNumCorners ? tCell.mNodeIndices( tEdge( 1 ) ) : tCenterNodeIndex;
                Matrix< DDRMat > tFirstXi    = tCell.mCornerXi( tEdge( 0 ) );
                Matrix< DDRMat > tSecondXi   = tEdge( 1 ) < tNumCorners ? tCell.mCornerXi( tEdge( 1 ) ) : Matrix< DDRMat >( 1, tNumDims, 0.0 );

                bool tIntersected = aGeometryEngine.queue_intersection(
                        tFirstNode, tSecondNode, tFirstXi, tSecondXi, tCell.mNodeIndices, tCell.mGeometryType, mtk::Interpolation_Order::LINEAR );
                if ( not tIntersected
                        or aGeometryEngine.queued_intersection_first_parent_on_interface()
                        or aGeometryEngine.queued_intersection_second_parent_on_interface() )
                {
                    continue;
                }

                Matrix< DDRMat > tCoordinates = aGeometryEngine.get_queued_intersection_global_coordinates();
                tStatistics.mMaxDistance      = std::max( tStatistics.mMaxDistance, distance_to_round( tCoordinates ) );
                tStatistics.mNumIntersections++;
                tStatistics.mLocalCoordinates.push_back( aGeometryEngine.get_queued_intersection_local_coordinate() );

                if ( not aAdmitAndCheckFD )
                {
                    continue;
                }

                // analytic d(coordinates)/d(ADV), summed per ADV ID
                uint tIntersectionIndex = tNodeManager.get_total_number_of_nodes();
                aGeometryEngine.admit_queued_intersection();
                Matrix< DDRMat > tI;
                eye( tNumDims, tNumDims, tI );
                Matrix< DDRMat > tAnalytic( 0, 0 );
                tNodeManager.append_dcoordinate_dadv_from_derived_node( tIntersectionIndex, tAnalytic, tI );
                Vector< sint > tADVIDs = tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( tIntersectionIndex );
                REQUIRE( tAnalytic.n_cols() == tADVIDs.size() );

                Vector< real >&  tADVs = aGeometryEngine.get_advs();
                uint             tNumADVs = tADVs.size();
                Matrix< DDRMat > tAnalyticPerADV( tNumDims, tNumADVs, 0.0 );
                for ( uint iColumn = 0; iColumn < tADVIDs.size(); iColumn++ )
                {
                    REQUIRE( tADVIDs( iColumn ) >= 0 );
                    REQUIRE( tADVIDs( iColumn ) < static_cast< sint >( tNumADVs ) );
                    for ( uint iDim = 0; iDim < tNumDims; iDim++ )
                    {
                        tAnalyticPerADV( iDim, tADVIDs( iColumn ) ) += tAnalytic( iDim, iColumn );
                    }
                }

                // central finite differences of the re-queued intersection, every ADV
                real           tEpsilon = 1e-6;
                Vector< real > tNominal = tADVs;
                for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
                {
                    Vector< real > tPerturbed = tNominal;
                    tPerturbed( iADV ) += tEpsilon;
                    aGeometryEngine.set_advs( tPerturbed );
                    aGeometryEngine.queue_intersection(
                            tFirstNode, tSecondNode, tFirstXi, tSecondXi, tCell.mNodeIndices, tCell.mGeometryType, mtk::Interpolation_Order::LINEAR );
                    Matrix< DDRMat > tPositive = aGeometryEngine.get_queued_intersection_global_coordinates();

                    tPerturbed( iADV ) -= 2.0 * tEpsilon;
                    aGeometryEngine.set_advs( tPerturbed );
                    aGeometryEngine.queue_intersection(
                            tFirstNode, tSecondNode, tFirstXi, tSecondXi, tCell.mNodeIndices, tCell.mGeometryType, mtk::Interpolation_Order::LINEAR );
                    Matrix< DDRMat > tNegative = aGeometryEngine.get_queued_intersection_global_coordinates();

                    aGeometryEngine.set_advs( tNominal );

                    for ( uint iDim = 0; iDim < tNumDims; iDim++ )
                    {
                        real tFD    = ( tPositive( iDim ) - tNegative( iDim ) ) / ( 2.0 * tEpsilon );
                        real tError = std::abs( tAnalyticPerADV( iDim, iADV ) - tFD ) / std::max( 1.0, std::abs( tFD ) );
                        *aMaxFDError = std::max( *aMaxFDError, tError );
                        INFO( "cell " << iCell << ", dimension " << iDim << ", ADV " << iADV
                                      << ": analytic " << tAnalyticPerADV( iDim, iADV ) << ", finite difference " << tFD );
                        CHECK( tError < 1e-6 );
                    }
                }
                ( *aNumFDChecked )++;
            }
        }
        return tStatistics;
    }

    //--------------------------------------------------------------------------------------------------------------

    static Cut_Statistics
    cut_round( uint aNumDims, uint aNumElements, sint aAnalytic, uint aMaxIterations = 0 )
    {
        mtk::Interpolation_Mesh* tMesh = aNumDims == 2 ? create_simple_mesh( aNumElements, aNumElements ) : create_cube_mesh( aNumElements );

        Submodule_Parameter_Lists  tFieldParameterLists = create_round_parameters( aNumDims, aAnalytic, false, aMaxIterations );
        Geometry_Engine_Parameters tGeometryEngineParameters;
        ADV_Manager                tADVManager;
        Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
        tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
        Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

        Cut_Statistics tStatistics = cut_all_cells( tMesh, tGeometryEngine );

        delete tMesh;
        return tStatistics;
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Analytic Intersection Precision", "[gen], [intersection], [analytic_intersection]" )
    {
        if ( par_size() == 1 )
        {
            Derived_Node::clear_override_linear_interpolation();

            SECTION( "2D circle: intersections on the circle at every h; the linear cut converges at O(h^2)" )
            {
                Cut_Statistics tLinear4    = cut_round( 2, 4, 0 );
                Cut_Statistics tLinear8    = cut_round( 2, 8, 0 );
                Cut_Statistics tLinear16   = cut_round( 2, 16, 0 );
                Cut_Statistics tAnalytic4  = cut_round( 2, 4, 1 );
                Cut_Statistics tAnalytic8  = cut_round( 2, 8, 1 );
                Cut_Statistics tAnalytic16 = cut_round( 2, 16, 1 );

                real tLinearRate = std::log2( tLinear8.mMaxDistance / tLinear16.mMaxDistance );
                std::printf(
                        "analytic intersection, 2D circle: max distance to the circle (h = 0.5 / 0.25 / 0.125): linear %.3e / %.3e / %.3e "
                        "(rate %.2f), analytic %.3e / %.3e / %.3e\n",
                        tLinear4.mMaxDistance,
                        tLinear8.mMaxDistance,
                        tLinear16.mMaxDistance,
                        tLinearRate,
                        tAnalytic4.mMaxDistance,
                        tAnalytic8.mMaxDistance,
                        tAnalytic16.mMaxDistance );

                // same edges are cut
                REQUIRE( tAnalytic16.mNumIntersections > 0 );
                CHECK( tAnalytic4.mNumIntersections == tLinear4.mNumIntersections );
                CHECK( tAnalytic8.mNumIntersections == tLinear8.mNumIntersections );
                CHECK( tAnalytic16.mNumIntersections == tLinear16.mNumIntersections );

                // linear: second order in h
                CHECK( tLinearRate > 1.5 );
                CHECK( tLinearRate < 2.5 );

                // analytic: on the circle to the root tolerance at every h, at least 1e6 times closer than linear
                for ( const Cut_Statistics* tPair : { &tAnalytic4, &tAnalytic8, &tAnalytic16 } )
                {
                    CHECK( tPair->mMaxDistance < 1e-12 );
                }
                CHECK( tAnalytic4.mMaxDistance < 1e-6 * tLinear4.mMaxDistance );
                CHECK( tAnalytic8.mMaxDistance < 1e-6 * tLinear8.mMaxDistance );
                CHECK( tAnalytic16.mMaxDistance < 1e-6 * tLinear16.mMaxDistance );
            }

            SECTION( "3D sphere: intersections on the sphere at every h; the linear cut converges at O(h^2)" )
            {
                Cut_Statistics tLinear4   = cut_round( 3, 4, 0 );
                Cut_Statistics tLinear8   = cut_round( 3, 8, 0 );
                Cut_Statistics tAnalytic4 = cut_round( 3, 4, 1 );
                Cut_Statistics tAnalytic8 = cut_round( 3, 8, 1 );

                real tLinearRate = std::log2( tLinear4.mMaxDistance / tLinear8.mMaxDistance );
                std::printf(
                        "analytic intersection, 3D sphere: max distance to the sphere (h = 0.5 / 0.25): linear %.3e / %.3e (rate %.2f), "
                        "analytic %.3e / %.3e\n",
                        tLinear4.mMaxDistance,
                        tLinear8.mMaxDistance,
                        tLinearRate,
                        tAnalytic4.mMaxDistance,
                        tAnalytic8.mMaxDistance );

                REQUIRE( tAnalytic8.mNumIntersections > 0 );
                CHECK( tAnalytic4.mNumIntersections == tLinear4.mNumIntersections );
                CHECK( tAnalytic8.mNumIntersections == tLinear8.mNumIntersections );

                CHECK( tLinearRate > 1.5 );
                CHECK( tLinearRate < 2.5 );

                CHECK( tAnalytic4.mMaxDistance < 1e-12 );
                CHECK( tAnalytic8.mMaxDistance < 1e-12 );
                CHECK( tAnalytic4.mMaxDistance < 1e-6 * tLinear4.mMaxDistance );
                CHECK( tAnalytic8.mMaxDistance < 1e-6 * tLinear8.mMaxDistance );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Analytic Intersection Sensitivities", "[gen], [pdv], [intersection], [analytic_intersection]" )
    {
        if ( par_size() == 1 )
        {
            Derived_Node::clear_override_linear_interpolation();

            // every cut of every cell, the circle / sphere center and radius as ADVs: analytic vs central FD
            for ( uint iNumDims = 2; iNumDims <= 3; iNumDims++ )
            {
                mtk::Interpolation_Mesh* tMesh = iNumDims == 2 ? create_simple_mesh( 4, 4 ) : create_cube_mesh( 3 );

                Submodule_Parameter_Lists  tFieldParameterLists = create_round_parameters( iNumDims, 1, true );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                ADV_Manager                tADVManager;
                Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                tGeometryEngineParameters.mADVManager = tADVManager;
                tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );
                REQUIRE( tGeometryEngine.get_advs().size() == iNumDims + 1 );

                real tMaxError  = 0.0;
                uint tNumChecked = 0;
                cut_all_cells( tMesh, tGeometryEngine, true, &tMaxError, &tNumChecked );

                std::printf( "analytic intersection, %uD: %u intersection nodes, max |analytic - FD| / max( 1, |FD| ) = %.2e\n", iNumDims, tNumChecked, tMaxError );
                CHECK( tNumChecked > 8 );
                CHECK( tMaxError < 1e-6 );

                delete tMesh;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * A chain of nested intersections on the single [-1, 1]^2 element: circle 0 cuts edges (0, 1) and (0, 3) (nodes 4
     * and 5), circle 1 cuts the segment between them (node 6), circle 2 cuts (6, corner 2) (node 7). Corner order
     * 0: (-1, -1), 1: (1, -1), 2: (1, 1), 3: (-1, 1). The ADVs are { c_x, c_y, r } per circle.
     *
     * @param aCircleADVs ADVs, 3 per circle
     * @param aAnalytic Per circle: 1 use_analytic_intersection, 2 the same with analytic_intersection_max_iterations = 1 (every root falls back to the linear crossing)
     * @param aCoordinates Global coordinates of nodes 4-7
     * @param aSensitivities Dense dX/dADV ( 2 x 9 ) of nodes 4-7 (only filled if aWithSensitivities)
     */
    static void
    build_nested_circle_intersections(
            mtk::Interpolation_Mesh*    aMesh,
            const Vector< real >&       aCircleADVs,
            const Vector< uint >&       aAnalytic,
            Vector< Matrix< DDRMat > >& aCoordinates,
            Vector< Matrix< DDRMat > >& aSensitivities,
            bool                        aWithSensitivities )
    {
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        for ( uint iCircle = 0; iCircle < 3; iCircle++ )
        {
            real tX = aCircleADVs( 3 * iCircle );
            real tY = aCircleADVs( 3 * iCircle + 1 );
            real tR = aCircleADVs( 3 * iCircle + 2 );
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            tFieldParameterLists.set( "center_x", tX, tX, tX );
            tFieldParameterLists.set( "center_y", tY, tY, tY );
            tFieldParameterLists.set( "radius", tR, tR, tR );
            tFieldParameterLists.set( "use_analytic_intersection", aAnalytic( iCircle ) >= 1 );
            if ( aAnalytic( iCircle ) == 2 )
            {
                // one root iteration: the edge falls back to the linear crossing
                tFieldParameterLists.set( "analytic_intersection_max_iterations", 1u );
            }
        }

        Geometry_Engine_Parameters tGeometryEngineParameters;
        ADV_Manager                tADVManager;
        Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
        tGeometryEngineParameters.mADVManager = tADVManager;
        tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
        Geometry_Engine_Test tGeometryEngine( aMesh, tGeometryEngineParameters );
        REQUIRE( tGeometryEngine.get_advs().size() == 9 );

        Cell_Data     tCell        = get_cell_data( aMesh, 0 );
        Node_Manager& tNodeManager = tGeometryEngine.get_node_manager();

        // corner k of the docstring = the mesh corner with that parametric position
        Vector< Matrix< DDRMat > > tCornerPositions = { { { -1.0, -1.0 } }, { { 1.0, -1.0 } }, { { 1.0, 1.0 } }, { { -1.0, 1.0 } } };
        Vector< uint >             tCornerNodes( 4 );
        for ( uint iCorner = 0; iCorner < 4; iCorner++ )
        {
            for ( uint iNode = 0; iNode < 4; iNode++ )
            {
                if ( tCell.mCornerXi( iNode )( 0 ) == tCornerPositions( iCorner )( 0 ) and tCell.mCornerXi( iNode )( 1 ) == tCornerPositions( iCorner )( 1 ) )
                {
                    tCornerNodes( iCorner ) = tCell.mNodeIndices( iNode );
                }
            }
        }

        // { first, second, circle }; a node index of 4 or more is a derived node (edge i creates node 4 + i)
        Vector< Vector< uint > > tEdges = { { 0, 1, 0 }, { 0, 3, 0 }, { 4, 5, 1 }, { 6, 2, 2 } };
        uint                     tActiveCircle = 0;
        for ( uint iEdge = 0; iEdge < tEdges.size(); iEdge++ )
        {
            while ( tActiveCircle < tEdges( iEdge )( 2 ) )
            {
                tGeometryEngine.advance_geometry_index();
                tActiveCircle++;
            }
            uint tFirst  = tEdges( iEdge )( 0 ) < 4 ? tCornerNodes( tEdges( iEdge )( 0 ) ) : tEdges( iEdge )( 0 );
            uint tSecond = tEdges( iEdge )( 1 ) < 4 ? tCornerNodes( tEdges( iEdge )( 1 ) ) : tEdges( iEdge )( 1 );

            // parametric and global coordinates coincide on the single [-1, 1]^2 element
            const Matrix< DDRMat >& tFirstX  = tNodeManager.get_node( tFirst ).get_global_coordinates();
            const Matrix< DDRMat >& tSecondX = tNodeManager.get_node( tSecond ).get_global_coordinates();
            Matrix< DDRMat >        tFirstXi = { { tFirstX( 0 ), tFirstX( 1 ) } };
            Matrix< DDRMat >        tSecondXi = { { tSecondX( 0 ), tSecondX( 1 ) } };

            bool tQueued = tGeometryEngine.queue_intersection(
                    tFirst, tSecond, tFirstXi, tSecondXi, tCell.mNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
            REQUIRE( tQueued );
            REQUIRE( not tGeometryEngine.queued_intersection_first_parent_on_interface() );
            REQUIRE( not tGeometryEngine.queued_intersection_second_parent_on_interface() );
            REQUIRE( std::abs( tGeometryEngine.get_queued_intersection_local_coordinate() ) < 0.9 );
            tGeometryEngine.admit_queued_intersection();
            REQUIRE( tNodeManager.get_total_number_of_nodes() == 5 + iEdge );
        }

        aCoordinates.resize( 4 );
        aSensitivities.resize( 4 );
        for ( uint iNode = 0; iNode < 4; iNode++ )
        {
            aCoordinates( iNode ) = tNodeManager.get_node( 4 + iNode ).get_global_coordinates();
            if ( not aWithSensitivities )
            {
                continue;
            }
            Matrix< DDRMat > tI;
            eye( 2, 2, tI );
            Matrix< DDRMat > tSensitivities;
            tNodeManager.append_dcoordinate_dadv_from_derived_node( 4 + iNode, tSensitivities, tI );
            Vector< sint > tADVIDs = tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( 4 + iNode );
            REQUIRE( tSensitivities.n_cols() == tADVIDs.size() );
            aSensitivities( iNode ).set_size( 2, 9, 0.0 );
            for ( uint iColumn = 0; iColumn < tADVIDs.size(); iColumn++ )
            {
                REQUIRE( tADVIDs( iColumn ) >= 0 );
                REQUIRE( tADVIDs( iColumn ) < 9 );
                for ( uint iDim = 0; iDim < 2; iDim++ )
                {
                    aSensitivities( iNode )( iDim, tADVIDs( iColumn ) ) += tSensitivities( iDim, iColumn );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Analytic Intersection Nested Sensitivities", "[gen], [pdv], [intersection], [analytic_intersection]" )
    {
        if ( par_size() == 1 )
        {
            Derived_Node::clear_override_linear_interpolation();
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 1, 1 );

            Vector< real > tADVs = { -1.2, -1.1, 1.0, -1.5, 0.5, 1.35, 1.2, 1.3, 1.0 };

            // every circle analytic (1), and every circle analytic with every root falling back (2). A linear node on these parents (key off on circle 1, or on all) does not pass: the linear
            // path takes d phi / d x at a derived parent as the chord slope of its locators, not the exact field gradient
            // (MB-083, a separate defect of the linear path)
            for ( const Vector< uint >& tAnalytic : { Vector< uint >{ 1, 1, 1 }, Vector< uint >{ 2, 2, 2 } } )
            {
                Vector< Matrix< DDRMat > > tCoordinates;
                Vector< Matrix< DDRMat > > tAnalyticSensitivities;
                build_nested_circle_intersections( tMesh, tADVs, tAnalytic, tCoordinates, tAnalyticSensitivities, true );

                // analytic circles: the nodes lie on their circle
                Vector< uint > tNodeCircle = { 0, 0, 1, 2 };
                for ( uint iNode = 0; iNode < 4; iNode++ )
                {
                    uint tCircle = tNodeCircle( iNode );
                    if ( tAnalytic( tCircle ) == 1 )    // 2: linear fallback, not on the circle
                    {
                        real tDistance = std::abs( std::hypot( tCoordinates( iNode )( 0 ) - tADVs( 3 * tCircle ), tCoordinates( iNode )( 1 ) - tADVs( 3 * tCircle + 1 ) ) - tADVs( 3 * tCircle + 2 ) );
                        CHECK( tDistance < 1e-12 );
                    }
                }

                // central finite differences, every perturbed chain rebuilt from scratch
                real tStep     = 1e-6;
                real tMaxError = 0.0;
                for ( uint iADV = 0; iADV < 9; iADV++ )
                {
                    Vector< real > tPlus  = tADVs;
                    Vector< real > tMinus = tADVs;
                    tPlus( iADV ) += tStep;
                    tMinus( iADV ) -= tStep;
                    Vector< Matrix< DDRMat > > tPlusCoordinates;
                    Vector< Matrix< DDRMat > > tMinusCoordinates;
                    Vector< Matrix< DDRMat > > tUnused;
                    build_nested_circle_intersections( tMesh, tPlus, tAnalytic, tPlusCoordinates, tUnused, false );
                    build_nested_circle_intersections( tMesh, tMinus, tAnalytic, tMinusCoordinates, tUnused, false );

                    for ( uint iNode = 0; iNode < 4; iNode++ )
                    {
                        for ( uint iDim = 0; iDim < 2; iDim++ )
                        {
                            real tFD    = ( tPlusCoordinates( iNode )( iDim ) - tMinusCoordinates( iNode )( iDim ) ) / ( 2.0 * tStep );
                            real tValue = tAnalyticSensitivities( iNode )( iDim, iADV );
                            real tError = std::abs( tValue - tFD ) / std::max( 1.0, std::abs( tFD ) );
                            tMaxError   = std::max( tMaxError, tError );
                            INFO( "node " << 4 + iNode << ", dimension " << iDim << ", ADV " << iADV << ": analytic " << tValue << ", finite difference " << tFD );
                            CHECK( tError < 1e-6 );
                        }
                    }
                }

                // the chain is non-trivial: node 7 depends on all three circles
                for ( uint iCircle = 0; iCircle < 3; iCircle++ )
                {
                    real tSum = 0.0;
                    for ( uint iParameter = 0; iParameter < 3; iParameter++ )
                    {
                        tSum += std::abs( tAnalyticSensitivities( 3 )( 0, 3 * iCircle + iParameter ) ) + std::abs( tAnalyticSensitivities( 3 )( 1, 3 * iCircle + iParameter ) );
                    }
                    CHECK( tSum > 1e-3 );
                }

                std::printf( "analytic intersection, nested chain (analytic %d %d %d): max |analytic - FD| / max( 1, |FD| ) = %.2e\n",
                        int( tAnalytic( 0 ) ),
                        int( tAnalytic( 1 ) ),
                        int( tAnalytic( 2 ) ),
                        tMaxError );
            }

            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Line field that is NaN in a band around its zero: the analytic root cannot be bracketed there.
     */
    class NaN_Band_Line : public Line
    {
      public:
        NaN_Band_Line( const ADV& aXCenter, const ADV& aYCenter, const ADV& aXNormal, const ADV& aYNormal, real aHalfWidth )
                : Line( aXCenter, aYCenter, aXNormal, aYNormal )
                , mHalfWidth( aHalfWidth )
        {
        }

        real get_field_value( const Matrix< DDRMat >& aCoordinates ) override
        {
            real tValue = Line::get_field_value( aCoordinates );
            return std::abs( tValue ) < mHalfWidth ? std::numeric_limits< real >::quiet_NaN() : tValue;
        }

      private:
        real mHalfWidth;
    };

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Analytic Intersection Default Off and Fallback", "[gen], [intersection], [analytic_intersection]" )
    {
        if ( par_size() == 1 )
        {
            Derived_Node::clear_override_linear_interpolation();

            SECTION( "The key is off by default" )
            {
                CHECK( not prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ).get< bool >( "use_analytic_intersection" ) );
                CHECK( not Level_Set_Parameters().mUseAnalyticIntersection );
            }

            SECTION( "Off: bit-identical to the linear crossing of the end values" )
            {
                Cut_Statistics tDefault  = cut_round( 2, 8, -1 );
                Cut_Statistics tOff      = cut_round( 2, 8, 0 );
                Cut_Statistics tAnalytic = cut_round( 2, 8, 1 );

                REQUIRE( tDefault.mNumIntersections > 0 );
                REQUIRE( tDefault.mLocalCoordinates.size() == tOff.mLocalCoordinates.size() );
                REQUIRE( tDefault.mLocalCoordinates.size() == tAnalytic.mLocalCoordinates.size() );
                real tMaxChange = 0.0;
                for ( uint iCut = 0; iCut < tDefault.mLocalCoordinates.size(); iCut++ )
                {
                    CHECK( tDefault.mLocalCoordinates( iCut ) == tOff.mLocalCoordinates( iCut ) );
                    tMaxChange = std::max( tMaxChange, std::abs( tAnalytic.mLocalCoordinates( iCut ) - tOff.mLocalCoordinates( iCut ) ) );
                }
                CHECK( tDefault.mMaxDistance == tOff.mMaxDistance );

                // the key changes the cut
                CHECK( tMaxChange > 1e-3 );

                // and the off path is exactly xi = ( 2 c - phi_2 - phi_1 ) / ( phi_2 - phi_1 ) of the end values
                Circle           tCircle( tCenterX, tCenterY, tRadius );
                Matrix< DDRMat > tFirst  = { { -1.0, 0.0 } };
                Matrix< DDRMat > tSecond = { { 0.0, 0.0 } };
                real             tPhi1   = tCircle.get_field_value( tFirst );
                real             tPhi2   = tCircle.get_field_value( tSecond );
                real             tLinear = ( 2.0 * 0.0 - tPhi2 - tPhi1 ) / ( tPhi2 - tPhi1 );

                // the node stores the parent basis weight 0.5 ( 1 - xi ) and returns 1 - 2 * weight
                real tStoredLinear = 1.0 - 2.0 * ( 0.5 * ( 1.0 - tLinear ) );

                mtk::Interpolation_Mesh*   tMesh                = create_simple_mesh( 2, 2 );
                Submodule_Parameter_Lists  tFieldParameterLists = create_round_parameters( 2, -1, false );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                ADV_Manager                tADVManager;
                Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

                // the edge from (-1, 0) to (0, 0) of the 2 x 2 mesh
                bool tFound = false;
                for ( uint iCell = 0; iCell < tMesh->get_num_elems() and not tFound; iCell++ )
                {
                    Cell_Data tCell = get_cell_data( tMesh, iCell );
                    for ( uint iFirst = 0; iFirst < 4; iFirst++ )
                    {
                        for ( uint iSecond = 0; iSecond < 4; iSecond++ )
                        {
                            Matrix< DDRMat > tX1 = tMesh->get_node_coordinate( tCell.mNodeIndices( iFirst ) );
                            Matrix< DDRMat > tX2 = tMesh->get_node_coordinate( tCell.mNodeIndices( iSecond ) );
                            if ( tX1( 0 ) == -1.0 and tX1( 1 ) == 0.0 and tX2( 0 ) == 0.0 and tX2( 1 ) == 0.0 and not tFound )
                            {
                                REQUIRE( tGeometryEngine.queue_intersection( tCell.mNodeIndices( iFirst ), tCell.mNodeIndices( iSecond ),
                                        tCell.mCornerXi( iFirst ), tCell.mCornerXi( iSecond ), tCell.mNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR ) );
                                CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() == tStoredLinear );
                                tFound = true;
                            }
                        }
                    }
                }
                CHECK( tFound );
                delete tMesh;
            }

            SECTION( "Analytic needs exact derived-node values and excludes multilinear / cell-order intersections" )
            {
                Level_Set_Parameters tParameters;
                tParameters.mUseAnalyticIntersection     = true;
                tParameters.mUseMultilinearInterpolation = true;
                CHECK_THROWS( Level_Set_Geometry( std::make_shared< Circle >( tCenterX, tCenterY, tRadius ), tParameters ) );

                tParameters.mUseMultilinearInterpolation = false;
                CHECK_NOTHROW( Level_Set_Geometry( std::make_shared< Circle >( tCenterX, tCenterY, tRadius ), tParameters ) );
            }

            SECTION( "A field that is not finite inside the edge falls back to the linear crossing and is counted" )
            {
                // line x = 0.3 (normal (1, 0)), NaN for |x - 0.3| < 0.05: every interior iterate near the root is NaN
                mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 1, 1 );

                Level_Set_Parameters tParameters;
                tParameters.mUseAnalyticIntersection = true;
                std::shared_ptr< Level_Set_Geometry > tGeometry =
                        std::make_shared< Level_Set_Geometry >( std::make_shared< NaN_Band_Line >( 0.3, 0.0, 1.0, 0.0, 0.05 ), tParameters );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                tGeometryEngineParameters.mGeometries = { std::shared_ptr< Geometry >( tGeometry ) };
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

                Cell_Data tCell = get_cell_data( tMesh, 0 );
                uint      tNumQueued = 0;
                for ( uint iFirst = 0; iFirst < 4; iFirst++ )
                {
                    for ( uint iSecond = iFirst + 1; iSecond < 4; iSecond++ )
                    {
                        // horizontal edges (y equal) cross x = 0.3
                        if ( tCell.mCornerXi( iFirst )( 1 ) != tCell.mCornerXi( iSecond )( 1 ) or tCell.mCornerXi( iFirst )( 0 ) == tCell.mCornerXi( iSecond )( 0 ) )
                        {
                            continue;
                        }
                        REQUIRE( tGeometryEngine.queue_intersection( tCell.mNodeIndices( iFirst ), tCell.mNodeIndices( iSecond ),
                                tCell.mCornerXi( iFirst ), tCell.mCornerXi( iSecond ), tCell.mNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR ) );

                        // the linear crossing of a line is exact
                        Matrix< DDRMat > tX = tGeometryEngine.get_queued_intersection_global_coordinates();
                        CHECK( tX( 0 ) == Approx( 0.3 ).margin( 1e-14 ) );
                        tNumQueued++;
                    }
                }
                REQUIRE( tNumQueued == 2 );
                CHECK( tGeometry->get_number_of_analytic_intersections() == 2 );
                CHECK( tGeometry->get_number_of_analytic_intersection_fallbacks() == 2 );

                delete tMesh;
            }

            SECTION( "Edges without a strict sign change take the linear path; a converged root is not a fallback" )
            {
                mtk::Interpolation_Mesh*   tMesh                = create_simple_mesh( 4, 4 );
                Submodule_Parameter_Lists  tFieldParameterLists = create_round_parameters( 2, 1, false );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                ADV_Manager                tADVManager;
                Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

                Cut_Statistics tStatistics = cut_all_cells( tMesh, tGeometryEngine );
                auto           tGeometry   = std::dynamic_pointer_cast< Level_Set_Geometry >( tGeometryEngineParameters.mGeometries( 0 ) );
                REQUIRE( tGeometry != nullptr );
                CHECK( tGeometry->get_number_of_analytic_intersections() >= tStatistics.mNumIntersections );
                CHECK( tGeometry->get_number_of_analytic_intersection_fallbacks() == 0 );

                delete tMesh;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Line field with a flat zero plateau around its zero: phi = v - w sign( v ) for |v| > w, else 0 (v the line field).
     * A root inside the plateau has d phi / d xi = 0.
     */
    class Plateau_Line : public Line
    {
      public:
        Plateau_Line( real aXCenter, real aYCenter, real aXNormal, real aYNormal, real aHalfWidth )
                : Line( aXCenter, aYCenter, aXNormal, aYNormal )
                , mHalfWidth( aHalfWidth )
        {
        }

        real get_field_value( const Matrix< DDRMat >& aCoordinates ) override
        {
            real tValue = Line::get_field_value( aCoordinates );
            return std::abs( tValue ) <= mHalfWidth ? 0.0 : tValue - std::copysign( mHalfWidth, tValue );
        }

      private:
        real mHalfWidth;
    };

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Analytic Intersection Fallback Sensitivities", "[gen], [pdv], [intersection], [analytic_intersection]" )
    {
        if ( par_size() == 1 )
        {
            Derived_Node::clear_override_linear_interpolation();

            SECTION( "Iteration limit reached: linear crossing, and the sensitivities of the linear crossing (central FD)" )
            {
                for ( uint iNumDims = 2; iNumDims <= 3; iNumDims++ )
                {
                    // one iteration: no edge of the circle / sphere converges
                    Cut_Statistics tLinear   = cut_round( iNumDims, 4, 0 );
                    Cut_Statistics tFallback = cut_round( iNumDims, 4, 1, 1 );
                    REQUIRE( tFallback.mNumIntersections == tLinear.mNumIntersections );
                    REQUIRE( tFallback.mLocalCoordinates.size() == tLinear.mLocalCoordinates.size() );
                    for ( uint iCut = 0; iCut < tLinear.mLocalCoordinates.size(); iCut++ )
                    {
                        CHECK( tFallback.mLocalCoordinates( iCut ) == tLinear.mLocalCoordinates( iCut ) );
                    }

                    // ADVs, admit every intersection, FD of the re-queued (again falling back) coordinates
                    mtk::Interpolation_Mesh*   tMesh                = iNumDims == 2 ? create_simple_mesh( 4, 4 ) : create_cube_mesh( 3 );
                    Submodule_Parameter_Lists  tFieldParameterLists = create_round_parameters( iNumDims, 1, true, 1 );
                    Geometry_Engine_Parameters tGeometryEngineParameters;
                    ADV_Manager                tADVManager;
                    Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
                    tGeometryEngineParameters.mADVManager = tADVManager;
                    tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
                    Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );
                    auto                 tGeometry = std::dynamic_pointer_cast< Level_Set_Geometry >( tGeometryEngineParameters.mGeometries( 0 ) );
                    REQUIRE( tGeometry != nullptr );

                    real tMaxError   = 0.0;
                    uint tNumChecked = 0;
                    cut_all_cells( tMesh, tGeometryEngine, true, &tMaxError, &tNumChecked );

                    std::printf( "analytic intersection, %uD iteration-limit fallback: %u intersection nodes, %u of %u analytic roots fell back, "
                                 "max |sensitivity - FD| / max( 1, |FD| ) = %.2e\n",
                            iNumDims,
                            tNumChecked,
                            tGeometry->get_number_of_analytic_intersection_fallbacks(),
                            tGeometry->get_number_of_analytic_intersections(),
                            tMaxError );
                    CHECK( tNumChecked > 8 );
                    CHECK( tGeometry->get_number_of_analytic_intersections() > 0 );
                    CHECK( tGeometry->get_number_of_analytic_intersection_fallbacks() == tGeometry->get_number_of_analytic_intersections() );
                    CHECK( tMaxError < 1e-6 );

                    delete tMesh;
                }
            }

            SECTION( "Field not finite inside the edge: linear crossing, and the sensitivities of the linear crossing (central FD)" )
            {
                // line x = 0.3 with all four parameters as ADVs, NaN for |phi| < 0.05
                mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 1, 1 );
                ADV_Manager              tADVManager( 4 );
                Vector< ADV >            tLineADVs = {
                    tADVManager.create_adv( { 0.3, 0.3, 0.3 } ),
                    tADVManager.create_adv( { 0.1, 0.1, 0.1 } ),
                    tADVManager.create_adv( { 1.0, 1.0, 1.0 } ),
                    tADVManager.create_adv( { 0.2, 0.2, 0.2 } )
                };
                Level_Set_Parameters tParameters;
                tParameters.mUseAnalyticIntersection = true;
                std::shared_ptr< Level_Set_Geometry > tGeometry = std::make_shared< Level_Set_Geometry >(
                        std::make_shared< NaN_Band_Line >( tLineADVs( 0 ), tLineADVs( 1 ), tLineADVs( 2 ), tLineADVs( 3 ), 0.05 ), tParameters );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                tGeometryEngineParameters.mADVManager = tADVManager;
                tGeometryEngineParameters.mGeometries = { std::shared_ptr< Geometry >( tGeometry ) };
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );
                Node_Manager&        tNodeManager = tGeometryEngine.get_node_manager();
                REQUIRE( tGeometryEngine.get_advs().size() == 4 );

                Cell_Data tCell      = get_cell_data( tMesh, 0 );
                uint      tNumEdges  = 0;
                real      tMaxError  = 0.0;
                real      tEpsilon   = 1e-6;
                Matrix< DDRMat > tI;
                eye( 2, 2, tI );
                for ( uint iFirst = 0; iFirst < 4; iFirst++ )
                {
                    for ( uint iSecond = iFirst + 1; iSecond < 4; iSecond++ )
                    {
                        // horizontal edges cross the line
                        if ( tCell.mCornerXi( iFirst )( 1 ) != tCell.mCornerXi( iSecond )( 1 ) or tCell.mCornerXi( iFirst )( 0 ) == tCell.mCornerXi( iSecond )( 0 ) )
                        {
                            continue;
                        }
                        uint             tFirstNode  = tCell.mNodeIndices( iFirst );
                        uint             tSecondNode = tCell.mNodeIndices( iSecond );
                        Matrix< DDRMat > tFirstXi    = tCell.mCornerXi( iFirst );
                        Matrix< DDRMat > tSecondXi   = tCell.mCornerXi( iSecond );
                        uint             tFallbacks  = tGeometry->get_number_of_analytic_intersection_fallbacks();
                        REQUIRE( tGeometryEngine.queue_intersection( tFirstNode, tSecondNode, tFirstXi, tSecondXi, tCell.mNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR ) );
                        CHECK( tGeometry->get_number_of_analytic_intersection_fallbacks() == tFallbacks + 1 );

                        uint tIntersectionIndex = tNodeManager.get_total_number_of_nodes();
                        tGeometryEngine.admit_queued_intersection();
                        const Intersection_Node_Analytic* tNode = dynamic_cast< const Intersection_Node_Analytic* >( &static_cast< const Node_Manager& >( tNodeManager ).get_derived_node( tIntersectionIndex ) );
                        REQUIRE( tNode != nullptr );
                        CHECK( tNode->is_linear_fallback() );

                        Matrix< DDRMat > tAnalytic( 0, 0 );
                        tNodeManager.append_dcoordinate_dadv_from_derived_node( tIntersectionIndex, tAnalytic, tI );
                        Vector< sint > tADVIDs = tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( tIntersectionIndex );
                        REQUIRE( tAnalytic.n_cols() == tADVIDs.size() );
                        Matrix< DDRMat > tPerADV( 2, 4, 0.0 );
                        for ( uint iColumn = 0; iColumn < tADVIDs.size(); iColumn++ )
                        {
                            REQUIRE( tADVIDs( iColumn ) >= 0 );
                            REQUIRE( tADVIDs( iColumn ) < 4 );
                            tPerADV( 0, tADVIDs( iColumn ) ) += tAnalytic( 0, iColumn );
                            tPerADV( 1, tADVIDs( iColumn ) ) += tAnalytic( 1, iColumn );
                        }

                        Vector< real > tNominal = tGeometryEngine.get_advs();
                        real           tSum     = 0.0;
                        for ( uint iADV = 0; iADV < 4; iADV++ )
                        {
                            Vector< real > tPerturbed = tNominal;
                            tPerturbed( iADV ) += tEpsilon;
                            tGeometryEngine.set_advs( tPerturbed );
                            tGeometryEngine.queue_intersection( tFirstNode, tSecondNode, tFirstXi, tSecondXi, tCell.mNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
                            Matrix< DDRMat > tPositive = tGeometryEngine.get_queued_intersection_global_coordinates();
                            tPerturbed( iADV ) -= 2.0 * tEpsilon;
                            tGeometryEngine.set_advs( tPerturbed );
                            tGeometryEngine.queue_intersection( tFirstNode, tSecondNode, tFirstXi, tSecondXi, tCell.mNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
                            Matrix< DDRMat > tNegative = tGeometryEngine.get_queued_intersection_global_coordinates();
                            tGeometryEngine.set_advs( tNominal );
                            for ( uint iDim = 0; iDim < 2; iDim++ )
                            {
                                real tFD    = ( tPositive( iDim ) - tNegative( iDim ) ) / ( 2.0 * tEpsilon );
                                real tError = std::abs( tPerADV( iDim, iADV ) - tFD ) / std::max( 1.0, std::abs( tFD ) );
                                tMaxError   = std::max( tMaxError, tError );
                                tSum += std::abs( tFD );
                                INFO( "edge " << tNumEdges << ", dimension " << iDim << ", ADV " << iADV << ": sensitivity " << tPerADV( iDim, iADV ) << ", finite difference " << tFD );
                                CHECK( tError < 1e-6 );
                            }
                        }
                        CHECK( tSum > 1e-3 );
                        tNumEdges++;
                    }
                }
                std::printf( "analytic intersection, NaN-band fallback: %u edges, max |sensitivity - FD| / max( 1, |FD| ) = %.2e\n", tNumEdges, tMaxError );
                CHECK( tNumEdges == 2 );

                delete tMesh;
            }

            SECTION( "d phi / d xi = 0 at the root (flat plateau): the coordinate exists, the sensitivity request errors" )
            {
                // line x = 0.3, flat (phi = 0) for |x - 0.3| <= 0.1
                mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 1, 1 );
                Level_Set_Parameters     tParameters;
                tParameters.mUseAnalyticIntersection = true;
                std::shared_ptr< Level_Set_Geometry > tGeometry =
                        std::make_shared< Level_Set_Geometry >( std::make_shared< Plateau_Line >( 0.3, 0.0, 1.0, 0.0, 0.1 ), tParameters );
                Geometry_Engine_Parameters tGeometryEngineParameters;
                tGeometryEngineParameters.mGeometries = { std::shared_ptr< Geometry >( tGeometry ) };
                Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );
                Node_Manager&        tNodeManager = tGeometryEngine.get_node_manager();

                Cell_Data tCell     = get_cell_data( tMesh, 0 );
                uint      tNumEdges = 0;
                Matrix< DDRMat > tI;
                eye( 2, 2, tI );
                for ( uint iFirst = 0; iFirst < 4; iFirst++ )
                {
                    for ( uint iSecond = iFirst + 1; iSecond < 4; iSecond++ )
                    {
                        if ( tCell.mCornerXi( iFirst )( 1 ) != tCell.mCornerXi( iSecond )( 1 ) or tCell.mCornerXi( iFirst )( 0 ) == tCell.mCornerXi( iSecond )( 0 ) )
                        {
                            continue;
                        }
                        REQUIRE( tGeometryEngine.queue_intersection( tCell.mNodeIndices( iFirst ), tCell.mNodeIndices( iSecond ),
                                tCell.mCornerXi( iFirst ), tCell.mCornerXi( iSecond ), tCell.mNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR ) );
                        Matrix< DDRMat > tX = tGeometryEngine.get_queued_intersection_global_coordinates();
                        CHECK( std::abs( tX( 0 ) - 0.3 ) <= 0.1 );

                        uint tIntersectionIndex = tNodeManager.get_total_number_of_nodes();
                        tGeometryEngine.admit_queued_intersection();
                        const Intersection_Node_Analytic* tNode = dynamic_cast< const Intersection_Node_Analytic* >( &static_cast< const Node_Manager& >( tNodeManager ).get_derived_node( tIntersectionIndex ) );
                        REQUIRE( tNode != nullptr );
                        CHECK( not tNode->is_linear_fallback() );
                        Matrix< DDRMat > tSensitivities( 0, 0 );
                        CHECK_THROWS( tNode->append_dcoordinate_dadv( tSensitivities, tI ) );
                        tNumEdges++;
                    }
                }
                CHECK( tNumEdges == 2 );
                CHECK( tGeometry->get_number_of_analytic_intersection_fallbacks() == 0 );

                delete tMesh;
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
