/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Intersection_Cell_Order.cpp
 *
 * Cell-order intersections (level-set key use_cell_order_intersection): XTK edge intersections and derived node
 * values on the background cell's own-order Lagrange interpolant instead of its corner multilinear one.
 * See doc/internal/XTK_HIGH_FIDELITY_CUT.md.
 *
 */

#include "catch.hpp"
#include <cmath>
#include <algorithm>
#include <cstdio>
#include "fn_eye.hpp"

#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "cl_GEN_Derived_Node.hpp"
#include "cl_GEN_Node_Manager.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    // corners of a QUAD cell in parametric coordinates (QUAD4/QUAD9 node order)
    static Vector< Matrix< DDRMat > > tCellCornerXi = {
        { { -1.0, -1.0 } },
        { { 1.0, -1.0 } },
        { { 1.0, 1.0 } },
        { { -1.0, 1.0 } }
    };

    // circle used by the precision tests
    static const real tCenterX = 0.1;
    static const real tCenterY = -0.05;
    static const real tRadius  = 0.63;

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Statistics of the edge intersections XTK's 2D regular subdivision would query on every cell of a mesh: the 4
     * outer edges and the 4 spokes from the corners to the cell center (Regular_Subdivision_4_TRIS).
     */
    struct Cut_Statistics
    {
        real           mMaxDistance      = 0.0;    // max distance of an intersection from the analytic circle
        uint           mNumIntersections = 0;      // number of intersected edges
        Vector< real > mCoordinates;               // all intersection coordinates (x0, y0, x1, y1, ...) in query order
    };

    //--------------------------------------------------------------------------------------------------------------

    static Submodule_Parameter_Lists
    create_circle_parameters( bool aMultilinear, bool aCellOrder, sint aDiscretizationIndex )
    {
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
        tFieldParameterLists.set( "center_x", tCenterX );
        tFieldParameterLists.set( "center_y", tCenterY );
        tFieldParameterLists.set( "radius", tRadius );
        tFieldParameterLists.set( "discretization_mesh_index", aDiscretizationIndex );
        tFieldParameterLists.set( "use_multilinear_interpolation", aMultilinear );
        tFieldParameterLists.set( "use_cell_order_intersection", aCellOrder );
        return tFieldParameterLists;
    }

    //--------------------------------------------------------------------------------------------------------------

    static mtk::Interpolation_Order
    get_cell_order( uint aNumCellNodes )
    {
        return aNumCellNodes == 4 ? mtk::Interpolation_Order::LINEAR
             : aNumCellNodes == 9 ? mtk::Interpolation_Order::QUADRATIC
                                  : mtk::Interpolation_Order::CUBIC;
    }

    //--------------------------------------------------------------------------------------------------------------

    static Cut_Statistics
    cut_all_cells( mtk::Interpolation_Mesh* aMesh, Geometry_Engine_Test& aGeometryEngine )
    {
        Cut_Statistics tStatistics;
        Node_Manager&  tNodeManager = aGeometryEngine.get_node_manager();

        for ( uint iCell = 0; iCell < aMesh->get_num_elems(); iCell++ )
        {
            // all nodes of the background cell
            Matrix< IndexMat > tSignedNodeIndices = aMesh->get_nodes_connected_to_element_loc_inds( iCell );
            uint               tNumCellNodes      = tSignedNodeIndices.numel();
            Matrix< DDUMat >   tNodeIndices( tNumCellNodes, 1 );
            for ( uint iNode = 0; iNode < tNumCellNodes; iNode++ )
            {
                tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
            }
            mtk::Interpolation_Order tOrder = get_cell_order( tNumCellNodes );

            // center node of the 4-TRIS template, a derived node of this cell
            uint tCenterNodeIndex = tNodeManager.get_total_number_of_nodes();
            aGeometryEngine.create_new_derived_nodes( { tSignedNodeIndices }, { { { 0.0, 0.0 } } }, mtk::Geometry_Type::QUAD, tOrder );

            // outer edges and spokes
            for ( uint iEdge = 0; iEdge < 8; iEdge++ )
            {
                uint             tFirstNode  = iEdge < 4 ? tNodeIndices( iEdge ) : tNodeIndices( iEdge - 4 );
                uint             tSecondNode = iEdge < 4 ? tNodeIndices( ( iEdge + 1 ) % 4 ) : tCenterNodeIndex;
                Matrix< DDRMat > tFirstXi    = iEdge < 4 ? tCellCornerXi( iEdge ) : tCellCornerXi( iEdge - 4 );
                Matrix< DDRMat > tSecondXi   = iEdge < 4 ? tCellCornerXi( ( iEdge + 1 ) % 4 ) : Matrix< DDRMat >( { { 0.0, 0.0 } } );

                bool tIntersected = aGeometryEngine.queue_intersection(
                        tFirstNode, tSecondNode, tFirstXi, tSecondXi, tNodeIndices, mtk::Geometry_Type::QUAD, tOrder );

                if ( tIntersected
                        and not aGeometryEngine.queued_intersection_first_parent_on_interface()
                        and not aGeometryEngine.queued_intersection_second_parent_on_interface() )
                {
                    Matrix< DDRMat > tCoordinates = aGeometryEngine.get_queued_intersection_global_coordinates();
                    real             tDistance    = std::abs( std::hypot( tCoordinates( 0 ) - tCenterX, tCoordinates( 1 ) - tCenterY ) - tRadius );
                    tStatistics.mMaxDistance      = std::max( tStatistics.mMaxDistance, tDistance );
                    tStatistics.mNumIntersections++;
                    tStatistics.mCoordinates.push_back( tCoordinates( 0 ) );
                    tStatistics.mCoordinates.push_back( tCoordinates( 1 ) );
                }
            }
        }
        return tStatistics;
    }

    //--------------------------------------------------------------------------------------------------------------

    static Cut_Statistics
    cut_circle( uint aNumElements, uint aOrder, bool aMultilinear, bool aCellOrder )
    {
        mtk::Interpolation_Mesh* tMesh = create_simple_mesh( aNumElements, aNumElements, aOrder, aOrder );

        Submodule_Parameter_Lists  tFieldParameterLists = create_circle_parameters( aMultilinear, aCellOrder, -2 );
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

    TEST_CASE( "Cell-Order Intersection Precision", "[gen], [intersection], [cell_order_intersection]" )
    {
        if ( par_size() == 1 )
        {
            // derived nodes are forced to linear interpolation by any earlier multilinear geometry in this process
            Derived_Node::clear_override_linear_interpolation();

            SECTION( "Order 2: intersections converge to the circle at O(h^3) instead of O(h^2)" )
            {
                // cell-order runs first: a multilinear geometry sets the global linear override
                Cut_Statistics tCellOrder4 = cut_circle( 4, 2, false, true );
                Cut_Statistics tCellOrder8 = cut_circle( 8, 2, false, true );
                Cut_Statistics tLinear4    = cut_circle( 4, 2, false, false );
                Cut_Statistics tLinear8    = cut_circle( 8, 2, false, false );
                Cut_Statistics tBilinear4  = cut_circle( 4, 2, true, false );
                Cut_Statistics tBilinear8  = cut_circle( 8, 2, true, false );
                Derived_Node::clear_override_linear_interpolation();

                std::printf(
                        "cell-order intersection: max distance to circle (h = 0.5 / 0.25): multilinear %.3e / %.3e, "
                        "linear %.3e / %.3e, cell order %.3e / %.3e\n",
                        tBilinear4.mMaxDistance,
                        tBilinear8.mMaxDistance,
                        tLinear4.mMaxDistance,
                        tLinear8.mMaxDistance,
                        tCellOrder4.mMaxDistance,
                        tCellOrder8.mMaxDistance );

                REQUIRE( tCellOrder4.mNumIntersections > 0 );
                REQUIRE( tCellOrder8.mNumIntersections > 0 );

                // 5x (h = 0.5, radius ~ h) and 10x (h = 0.25) closer to the circle than the corner-multilinear cut
                CHECK( tCellOrder4.mMaxDistance < 0.2 * tBilinear4.mMaxDistance );
                CHECK( tCellOrder8.mMaxDistance < 0.1 * tBilinear8.mMaxDistance );
                CHECK( tCellOrder8.mMaxDistance < tLinear8.mMaxDistance );

                // observed orders: ~2 for the corner-multilinear cut, ~3 for the QUAD9 interpolant
                real tBilinearRate  = std::log2( tBilinear4.mMaxDistance / tBilinear8.mMaxDistance );
                real tCellOrderRate = std::log2( tCellOrder4.mMaxDistance / tCellOrder8.mMaxDistance );
                std::printf( "cell-order intersection: observed rate multilinear %.2f, cell order %.2f\n", tBilinearRate, tCellOrderRate );
                CHECK( tBilinearRate < 2.5 );
                CHECK( tCellOrderRate > 2.5 );
            }

            SECTION( "Order 1: cell order is identical to multilinear" )
            {
                Cut_Statistics tCellOrder = cut_circle( 6, 1, false, true );
                Cut_Statistics tBilinear  = cut_circle( 6, 1, true, false );
                Derived_Node::clear_override_linear_interpolation();

                REQUIRE( tCellOrder.mNumIntersections == tBilinear.mNumIntersections );
                REQUIRE( tCellOrder.mCoordinates.size() == tBilinear.mCoordinates.size() );
                for ( uint iCoordinate = 0; iCoordinate < tCellOrder.mCoordinates.size(); iCoordinate++ )
                {
                    CHECK( tCellOrder.mCoordinates( iCoordinate ) == Approx( tBilinear.mCoordinates( iCoordinate ) ).margin( 1e-10 ) );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Cell-Order Intersection Sensitivities", "[gen], [pdv], [intersection], [cell_order_intersection]" )
    {
        if ( par_size() == 1 )
        {
            Derived_Node::clear_override_linear_interpolation();

            // QUAD9 cells, quadratic B-spline level set: the ADVs are the B-spline coefficients
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 3, 3, 2, 2 );

            Submodule_Parameter_Lists  tFieldParameterLists = create_circle_parameters( false, true, 0 );
            Geometry_Engine_Parameters tGeometryEngineParameters;
            ADV_Manager                tADVManager;
            Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );
            Node_Manager&        tNodeManager = tGeometryEngine.get_node_manager();

            real             tEpsilon = 1e-6;
            Matrix< DDRMat > tI;
            eye( 2, 2, tI );
            uint tNumChecked     = 0;
            real tMaxRelativeErr = 0.0;

            for ( uint iCell = 0; iCell < tMesh->get_num_elems(); iCell++ )
            {
                Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( iCell );
                Matrix< DDUMat >   tNodeIndices( 9, 1 );
                for ( uint iNode = 0; iNode < 9; iNode++ )
                {
                    tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
                }

                // center node of the 4-TRIS template
                uint tCenterNodeIndex = tNodeManager.get_total_number_of_nodes();
                tGeometryEngine.create_new_derived_nodes( { tSignedNodeIndices }, { { { 0.0, 0.0 } } }, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::QUADRATIC );

                for ( uint iEdge = 0; iEdge < 8; iEdge++ )
                {
                    uint             tFirstNode  = iEdge < 4 ? tNodeIndices( iEdge ) : tNodeIndices( iEdge - 4 );
                    uint             tSecondNode = iEdge < 4 ? tNodeIndices( ( iEdge + 1 ) % 4 ) : tCenterNodeIndex;
                    Matrix< DDRMat > tFirstXi    = iEdge < 4 ? tCellCornerXi( iEdge ) : tCellCornerXi( iEdge - 4 );
                    Matrix< DDRMat > tSecondXi   = iEdge < 4 ? tCellCornerXi( ( iEdge + 1 ) % 4 ) : Matrix< DDRMat >( { { 0.0, 0.0 } } );

                    bool tIntersected = tGeometryEngine.queue_intersection(
                            tFirstNode, tSecondNode, tFirstXi, tSecondXi, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::QUADRATIC );
                    if ( not tIntersected
                            or tGeometryEngine.queued_intersection_first_parent_on_interface()
                            or tGeometryEngine.queued_intersection_second_parent_on_interface() )
                    {
                        continue;
                    }

                    // admit and get the analytic d(coordinate)/d(ADV)
                    uint tIntersectionIndex = tNodeManager.get_total_number_of_nodes();
                    tGeometryEngine.admit_queued_intersection();
                    Matrix< DDRMat > tAnalytic( 0, 0 );
                    tNodeManager.append_dcoordinate_dadv_from_derived_node( tIntersectionIndex, tAnalytic, tI );
                    Vector< sint > tADVIDs = tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( tIntersectionIndex );
                    REQUIRE( tAnalytic.n_cols() == tADVIDs.size() );

                    // the analytic sensitivities hold one column per (cell node, ADV) contribution: sum them per ADV
                    Vector< sint >   tUniqueADVIDs;
                    Matrix< DDRMat > tAnalyticPerADV( 2, tADVIDs.size(), 0.0 );
                    for ( uint iColumn = 0; iColumn < tADVIDs.size(); iColumn++ )
                    {
                        uint tUniqueIndex = 0;
                        while ( tUniqueIndex < tUniqueADVIDs.size() and tUniqueADVIDs( tUniqueIndex ) != tADVIDs( iColumn ) )
                        {
                            tUniqueIndex++;
                        }
                        if ( tUniqueIndex == tUniqueADVIDs.size() )
                        {
                            tUniqueADVIDs.push_back( tADVIDs( iColumn ) );
                        }
                        tAnalyticPerADV( 0, tUniqueIndex ) += tAnalytic( 0, iColumn );
                        tAnalyticPerADV( 1, tUniqueIndex ) += tAnalytic( 1, iColumn );
                    }

                    // central finite differences: perturb each ADV and re-queue the same edge
                    Vector< real >   tADVs  = tGeometryEngine.get_advs();
                    real             tScale = 1e-3;
                    Matrix< DDRMat > tFiniteDifference( 2, tUniqueADVIDs.size(), 0.0 );
                    for ( uint iADV = 0; iADV < tUniqueADVIDs.size(); iADV++ )
                    {
                        sint tADVID = tUniqueADVIDs( iADV );

                        tADVs( tADVID - 1 ) += tEpsilon;
                        tGeometryEngine.set_advs( tADVs );
                        tGeometryEngine.queue_intersection(
                                tFirstNode, tSecondNode, tFirstXi, tSecondXi, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::QUADRATIC );
                        Matrix< DDRMat > tPositive = tGeometryEngine.get_queued_intersection_global_coordinates();

                        tADVs( tADVID - 1 ) -= 2.0 * tEpsilon;
                        tGeometryEngine.set_advs( tADVs );
                        tGeometryEngine.queue_intersection(
                                tFirstNode, tSecondNode, tFirstXi, tSecondXi, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::QUADRATIC );
                        Matrix< DDRMat > tNegative = tGeometryEngine.get_queued_intersection_global_coordinates();

                        tADVs( tADVID - 1 ) += tEpsilon;
                        tGeometryEngine.set_advs( tADVs );

                        tFiniteDifference( 0, iADV ) = ( tPositive( 0 ) - tNegative( 0 ) ) / ( 2.0 * tEpsilon );
                        tFiniteDifference( 1, iADV ) = ( tPositive( 1 ) - tNegative( 1 ) ) / ( 2.0 * tEpsilon );
                        tScale                       = std::max( { tScale, std::abs( tAnalyticPerADV( 0, iADV ) ), std::abs( tAnalyticPerADV( 1, iADV ) ) } );
                    }

                    // compare against the largest sensitivity of this node
                    for ( uint iADV = 0; iADV < tUniqueADVIDs.size(); iADV++ )
                    {
                        for ( uint iDim = 0; iDim < 2; iDim++ )
                        {
                            real tRelativeError = std::abs( tAnalyticPerADV( iDim, iADV ) - tFiniteDifference( iDim, iADV ) ) / tScale;
                            tMaxRelativeErr     = std::max( tMaxRelativeErr, tRelativeError );
                            CHECK( tRelativeError < 1e-6 );
                        }
                    }
                    tNumChecked++;
                }
            }

            std::printf( "cell-order intersection: %u intersection nodes, max |analytic - FD| / max|analytic| = %.2e\n", tNumChecked, tMaxRelativeErr );
            CHECK( tNumChecked > 8 );

            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
