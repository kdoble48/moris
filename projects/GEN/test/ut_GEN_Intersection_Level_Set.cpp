/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Intersection_Level_Set.cpp
 *
 */

#include "catch.hpp"
#include <algorithm>
#include <cmath>
#include "fn_eye.hpp"

#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Geometry_Engine_Test.hpp"
#include "cl_GEN_PDV_Host_Manager.hpp"
#include "cl_GEN_BSpline_Field.hpp"
#include "cl_GEN_Design_Factory.hpp"
#include "fn_check_equal.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "cl_GEN_Mesh_Field.hpp"
#include "cl_MTK_Mesh_Factory.hpp"

#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    static Vector< Matrix< DDRMat > > tQuadParametricCoordinates = {
        { { -1.0, -1.0 } },
        { { 1.0, -1.0 } },
        { { 1.0, 1.0 } },
        { { -1.0, 1.0 } }
    };

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    get_quad_local_coordinates( uint tNodeNumber )
    {
        real tRadians = M_PI * ( 5 + 2 * tNodeNumber ) / 4;
        return sqrt( 2 ) * Matrix< DDRMat >( { { cos( tRadians ), sin( tRadians ) } } );
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Linear Intersections", "[gen], [pdv], [intersection], [linear intersection]" )
    {
        if ( par_size() == 1 )
        {
            // Create mesh
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 2, 2 );

            // Circle
            Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            tFieldParameterLists.set( "center_x", -0.25 );
            tFieldParameterLists.set( "center_y", 0.0 );
            tFieldParameterLists.set( "radius", 0.7499999999 );
            tFieldParameterLists.set( "discretization_mesh_index", 0 );

            // Plane 1
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tFieldParameterLists.set( "center_x", 0.25, 0.25, 0.25 );
            tFieldParameterLists.set( "center_y", 0.0, 0.0, 0.0 );
            tFieldParameterLists.set( "normal_x", 1.0, 1.0, 1.0 );
            tFieldParameterLists.set( "normal_y", 0.0, 0.0, 0.0 );

            // Plane 2
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tFieldParameterLists.set( "center_x", 1.0 );
            tFieldParameterLists.set( "center_y", 0.0 );
            tFieldParameterLists.set( "normal_x", 1.0 );
            tFieldParameterLists.set( "normal_y", 0.0 );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            ADV_Manager tADVManager;
            Design_Factory tDesignFactory( tFieldParameterLists, tADVManager );
            tGeometryEngineParameters.mADVManager = tADVManager;
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

            // TODO ensure this writes the mesh/fields correctly instead of just relying on no errors being thrown
            tGeometryEngine.output_fields_on_mesh( tMesh, "intersection_test.exo" );

            // Solution for is_intersected() per geometry and per element
            Vector< Vector< bool > > tIsElementIntersected = {
                { true, true, true, true },      // Geometry 0
                { false, true, false, true },    // Geometry 1
                { false, true, false, true }
            };    // Geometry 2

            // Per geometry, per element, per edge
            Vector< Vector< Vector< bool > > > tIsEdgeIntersected = {
                // Geometry 0
                { { false, true, true, false },            // Element 0
                        { false, false, true, true },      // Element 1
                        { true, true, false, false },      // Element 2
                        { true, false, false, true } },    // Element 3
                // Geometry 1
                { { false, false, false, false },
                        { true, false, true, false },
                        { false, false, false, false },
                        { true, false, true, false } },
                // Geometry 2
                { { false, false, false, false },
                        { true, true, true, false },
                        { false, false, false, false },
                        { true, true, true, false } }
            };

            // Intersection coordinates
            real tFrac = 2.0 / ( 3.0 + sqrt( 17.0 ) );

            Matrix< DDRMat > tIntersectionLocalCoordinates = { { -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5, 0.0000000002, -0.0000000002, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0 } };

            Vector< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.5, 0.0 } },
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } },
                { { 0.5, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } },
                { { 0.25, -1.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 1.0 } },
                { { 0.25, -0.25 - tFrac / 4.0 } },
                { { 0.25, 0.25 + tFrac / 4.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 1.0 } }
            };

            // Check element intersections
            uint tIntersectionCount = 0;
            for ( uint tGeometryIndex = 0; tGeometryIndex < 3; tGeometryIndex++ )
            {
                for ( uint tElementIndex = 0; tElementIndex < 4; tElementIndex++ )
                {
                    // Node indices per element
                    Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );
                    Matrix< DDUMat >   tNodeIndices( 4, 1 );
                    for ( uint iNode = 0; iNode < 4; iNode++ )
                    {
                        tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
                    }

                    // Check edges for properly queued intersections
                    Matrix< DDRMat > tNodeCoordinates( 4, 2 );
                    for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                    {
                        // Node coordinates
                        Matrix< DDRMat > tFirstNodeCoordinates  = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
                        Matrix< DDRMat > tSecondNodeCoordinates = tMesh->get_node_coordinate( tNodeIndices( ( tNodeNumber + 1 ) % 4 ) );

                        // Queue intersection
                        bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                                tNodeIndices( tNodeNumber ),
                                tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                tQuadParametricCoordinates( tNodeNumber ),
                                tQuadParametricCoordinates( ( tNodeNumber + 1 ) % 4 ),
                                tNodeIndices,
                                mtk::Geometry_Type::QUAD,
                                mtk::Interpolation_Order::LINEAR );
                        REQUIRE( tIntersectionQueued == tIsEdgeIntersected( tGeometryIndex )( tElementIndex )( tNodeNumber ) );

                        // Check queued intersection
                        if ( tIntersectionQueued )
                        {
                            // Check parents
                            bool tFirstParentOnInterface  = false;
                            bool tSecondParentOnInterface = false;

                            // Parent nodes on interface
                            if ( tGeometryIndex == 2 and tFirstNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tFirstParentOnInterface = true;
                            }

                            if ( tGeometryIndex == 2 and tSecondNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tSecondParentOnInterface = true;
                            }
                            CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
                            CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

                            // See if local coordinate is a number
                            real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();
                            if ( not std::isnan( tLocalCoordinate ) )
                            {
                                // Check local coordinate; note reference values are for radius of exactly 0.75; thus higher margin needed
                                CHECK( tLocalCoordinate == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

                                // Check global coordinates
                                CHECK_EQUAL(
                                        tGeometryEngine.get_queued_intersection_global_coordinates(),
                                        tIntersectionGlobalCoordinates( tIntersectionCount ), );
                            }

                            // Admit intersection
                            tGeometryEngine.admit_queued_intersection();
                            tIntersectionCount++;
                        }

                        // Set node coordinates for element checking
                        tNodeCoordinates.set_row( tNodeNumber, tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) ) );
                    }

                    // Check with solution
                    bool tIsIntersected = tGeometryEngine.is_intersected( tGeometryIndex, tSignedNodeIndices, tNodeCoordinates );
                    CHECK( tIsIntersected == tIsElementIntersected( tGeometryIndex )( tElementIndex ) );
                }

                // Intersection on intersection
                if ( tGeometryIndex == 1 )
                {
                    // Queue intersection on intersection 1
                    bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                            9,
                            11,
                            { { -1.0, -tFrac } },
                            { { 0.0, 1.0 } },
                            { { 1, 4, 5, 2 } },
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Order::LINEAR );

                    // Check intersection on intersection 1
                    REQUIRE( tIntersectionQueued == true );
                    CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == false );
                    CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == false );
                    CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ) );
                    CHECK_EQUAL(
                            tGeometryEngine.get_queued_intersection_global_coordinates(),
                            tIntersectionGlobalCoordinates( tIntersectionCount ), );

                    // Admit intersection on intersection 1
                    tGeometryEngine.admit_queued_intersection();
                    tIntersectionCount++;

                    // Queue intersection on intersection 2
                    tIntersectionQueued = tGeometryEngine.queue_intersection(
                            11,
                            14,
                            { { 0.0, -1.0 } },
                            { { -1.0, tFrac } },
                            { { 2, 5, 8, 6 } },
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Order::LINEAR );

                    // Check intersection on intersection 1
                    REQUIRE( tIntersectionQueued == true );
                    CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == false );
                    CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == false );
                    CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ) );
                    CHECK_EQUAL(
                            tGeometryEngine.get_queued_intersection_global_coordinates(),
                            tIntersectionGlobalCoordinates( tIntersectionCount ), );

                    // Admit intersection on intersection 1
                    tGeometryEngine.admit_queued_intersection();
                    tIntersectionCount++;
                }

                // Advance geometry index
                if ( tGeometryIndex < 2 )
                {
                    tGeometryEngine.advance_geometry_index();
                }
            }

            // Check total number of intersections
            CHECK( tIntersectionCount == 20 );

            // Test the new child nodes on the level set field (geometry 0)
            CHECK( tGeometryEngine.get_geometric_region( 0, 9, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 10, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 11, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 12, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 13, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 14, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 15, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 16, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 17, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 18, { {} } ) == Geometric_Region::NEGATIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 19, { {} } ) == Geometric_Region::NEGATIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 20, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 21, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 22, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 23, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 25, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 26, { {} } ) == Geometric_Region::POSITIVE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 28, { {} } ) == Geometric_Region::POSITIVE );

            // Get the PDV host manager
            auto tPDVHostManager = dynamic_cast< PDV_Host_Manager* >( tGeometryEngine.get_design_variable_interface() );

            // Test that the new intersections have been added to the PDV host manager, but ONLY for the circle
            Vector< Matrix< DDRMat > > tPDVValues( 0 );
            Vector< Vector< bool > >     tIsActive( 0 );
            tPDVHostManager->get_ig_pdv_value(
                    { { 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28 } },
                    { PDV_Type::X_COORDINATE, PDV_Type::Y_COORDINATE },
                    tPDVValues,
                    tIsActive );

            // Background nodes
            for ( uint tNodeIndex = 0; tNodeIndex < 9; tNodeIndex++ )
            {
                // Check if not active
                CHECK( tIsActive( 0 )( tNodeIndex ) == false );
                CHECK( tIsActive( 1 )( tNodeIndex ) == false );
            }

            // Nodes on the circle and first plane interfaces (depend on ADVs, active)
            for ( uint tNodeIndex = 9; tNodeIndex < 23; tNodeIndex++ )
            {
                // Check if active
                CHECK( tIsActive( 0 )( tNodeIndex ) == true );
                CHECK( tIsActive( 1 )( tNodeIndex ) == true );

                // Check for node coordinates as PDV values
                CHECK( tPDVValues( 0 )( tNodeIndex ) == Approx( tIntersectionGlobalCoordinates( tNodeIndex - 9 )( 0 ) ) );
                CHECK( tPDVValues( 1 )( tNodeIndex ) == Approx( tIntersectionGlobalCoordinates( tNodeIndex - 9 )( 1 ) ) );
            }

            // Nodes on the second plane interface (inactive)
            for ( uint tNodeIndex = 23; tNodeIndex < 27; tNodeIndex++ )
            {
                // Check if not active
                CHECK( tIsActive( 0 )( tNodeIndex ) == false );
                CHECK( tIsActive( 1 )( tNodeIndex ) == false );
            }

            // Check sensitivities
            Vector< Matrix< DDRMat > > tIntersectionSensitivities = {
                { { 0.0, 0.0 }, { ( 9 + sqrt( 17 ) ) / 16, ( 3 * sqrt( 17 ) - 5 ) / 16 } },
                { { 0.0, 2.0 }, { 0.0, 0.0 } },
                { { -0.5, -0.5 }, { 0.0, 0.0 } },
                { { 0.0, 0.0 }, { ( 3 * sqrt( 17 ) - 5 ) / 16, ( 9 + sqrt( 17 ) ) / 16 } },
                { { 2.0, 0.0 }, { 0.0, 0.0 } },
                { { 0.0, 0.0 }, { -( 3 * sqrt( 17 ) - 5 ) / 16, -( 9 + sqrt( 17 ) ) / 16 } },
                { { -0.5, -0.5 }, { 0.0, 0.0 } },
                { { 0.0, 0.0 }, { -( 9 + sqrt( 17 ) ) / 16, -( 3 * sqrt( 17 ) - 5 ) / 16 } },
                { { 0.75, 0.0, 0.1875, 0.75, 0.25, 0.0, -0.1875, 0.25 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.25, 0.0, -0.1875, 0.0, 0.75, 0.0, 0.1875, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.75, 0.0, 0.1875, 0.0, 0.25, 0.0, -0.1875, 0.0 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.25, 0.0, -0.1875, -0.25, 0.75, 0.0, 0.1875, -0.75 }, { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 } },
                { { 0.5, 0.0, 0.125, ( 1 + sqrt( 17 ) ) / 16, 0.5, 0.0, -0.125, 0.0, 0.0, 0.0, 0.0, 0.0 }, { ( 1 + sqrt( 17 ) ) / 8, 0.0, ( 1 + sqrt( 17 ) ) / 32, ( 9 + sqrt( 17 ) ) / 32, ( 1 + sqrt( 17 ) ) / 8, 0.0, -( 1 + sqrt( 17 ) ) / 32, 0.0, ( 9 + sqrt( 17 ) ) / 32, ( 3 * sqrt( 17 ) - 5 ) / 32, ( 1 + sqrt( 17 ) ) / 16, ( 1 + sqrt( 17 ) ) / 16 } },
                { { 0.5, 0.0, -0.125, 0.0, 0.5, 0.0, 0.125, -( 1 + sqrt( 17 ) ) / 16, 0.0, 0.0, 0.0, 0.0 }, { -( 1 + sqrt( 17 ) ) / 8, 0.0, ( 1 + sqrt( 17 ) ) / 32, 0.0, -( 1 + sqrt( 17 ) ) / 8, 0.0, -( 1 + sqrt( 17 ) ) / 32, ( 9 + sqrt( 17 ) ) / 32, -( 1 + sqrt( 17 ) ) / 16, -( 1 + sqrt( 17 ) ) / 16, -( 3 * sqrt( 17 ) - 5 ) / 32, -( 9 + sqrt( 17 ) ) / 32 } }
            };

            Vector< Vector< sint > > tIntersectionIDs = {
                { 6, 9 },
                { 9, 8 },
                { 10, 9 },
                { 9, 6 },
                { 8, 9 },
                { 9, 12 },
                { 9, 10 },
                { 12, 9 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3 },
                { 0, 1, 2, 3, 0, 1, 2, 3, 6, 9, 10, 9 },
                { 0, 1, 2, 3, 0, 1, 2, 3, 10, 9, 9, 12 }
            };

            Matrix< DDRMat > tHostADVSensitivities;
            Matrix< DDRMat > tI;
            Node_Manager&    tNodeManager = tGeometryEngine.get_node_manager();
            for ( uint iNodeIndex = 9; iNodeIndex < 23; iNodeIndex++ )
            {
                tHostADVSensitivities.set_size( 0.0, 0.0 );
                eye( 2, 2, tI );
                tNodeManager.append_dcoordinate_dadv_from_derived_node( iNodeIndex, tHostADVSensitivities, tI );
                CHECK_EQUAL(
                        tHostADVSensitivities,
                        tIntersectionSensitivities( iNodeIndex - 9 ),
                        1E8, );
                CHECK_EQUAL(
                        tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( iNodeIndex ),
                        tIntersectionIDs( iNodeIndex - 9 ), );

                // MB-022: the sensitivities are linear in the incoming chain factor, D(A) = A * D(I), at every
                // nesting depth. A non-commuting A also checks the order of the recursive product (A * J, not J * A).
                // Nodes 21 and 22 have intersection-node parents; before the fix they ignored A in the parent terms.
                Matrix< DDRMat > tChainFactor = { { 1.0, 2.0 }, { -3.0, 0.5 } };
                Matrix< DDRMat > tFactoredSensitivities;
                tNodeManager.append_dcoordinate_dadv_from_derived_node( iNodeIndex, tFactoredSensitivities, tChainFactor );
                Matrix< DDRMat > tExpectedFactoredSensitivities = tChainFactor * tHostADVSensitivities;
                CHECK_EQUAL(
                        tFactoredSensitivities,
                        tExpectedFactoredSensitivities,
                        1E8, );
            }

            //------------------------------------------------------------------------------------------------------
            // Start second check
            //------------------------------------------------------------------------------------------------------

            // Set new ADVs, level set field now has no intersections
            mtk::Mesh_Pair tMeshPair( tMesh, create_integration_mesh_from_interpolation_mesh( mtk::MeshType::HMR, tMesh ) );
            Vector< real > tADVs = { 0.25, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 };
            tGeometryEngine.set_advs( tADVs );
            tGeometryEngine.reset_mesh_information( tMeshPair.get_interpolation_mesh() );

            // Solution for is_intersected() per geometry and per element
            tIsElementIntersected = {
                { false, false, false, false },    // Geometry 0
                { false, true, false, true },      // Geometry 1
                { false, true, false, true }
            };    // Geometry 2

            // Per geometry, per element, per edge
            tIsEdgeIntersected = {
                // Geometry 0
                { { false, false, false, false },            // Element 0
                        { false, false, false, false },      // Element 1
                        { false, false, false, false },      // Element 2
                        { false, false, false, false } },    // Element 3
                // Geometry 1
                { { false, false, false, false },
                        { true, false, true, false },
                        { false, false, false, false },
                        { true, false, true, false } },
                // Geometry 2
                { { false, false, false, false },
                        { true, true, true, false },
                        { false, false, false, false },
                        { true, true, true, false } }
            };

            // Intersection coordinates
            tIntersectionLocalCoordinates  = { { -0.5, 0.5, -0.5, 0.5, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0 } };
            tIntersectionGlobalCoordinates = {
                { { 0.25, -1.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 0.0 } },
                { { 0.25, 1.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, -1.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 0.0 } },
                { { 1.0, 1.0 } }
            };

            // Check element intersections
            tIntersectionCount = 0;
            for ( uint tGeometryIndex = 0; tGeometryIndex < 3; tGeometryIndex++ )
            {
                for ( uint tElementIndex = 0; tElementIndex < 4; tElementIndex++ )
                {
                    // Node indices per element
                    Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );
                    Matrix< DDUMat >   tNodeIndices( 4, 1 );
                    for ( uint iNode = 0; iNode < 4; iNode++ )
                    {
                        tNodeIndices( iNode ) = tSignedNodeIndices( iNode );
                    }

                    // Check edges for properly queued intersections
                    Matrix< DDRMat > tNodeCoordinates( 4, 2 );
                    for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                    {
                        // Node coordinates
                        Matrix< DDRMat > tFirstNodeCoordinates  = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
                        Matrix< DDRMat > tSecondNodeCoordinates = tMesh->get_node_coordinate( tNodeIndices( ( tNodeNumber + 1 ) % 4 ) );

                        // Queue intersection
                        bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                                tNodeIndices( tNodeNumber ),
                                tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                tQuadParametricCoordinates( tNodeNumber ),
                                tQuadParametricCoordinates( ( tNodeNumber + 1 ) % 4 ),
                                tNodeIndices,
                                mtk::Geometry_Type::QUAD,
                                mtk::Interpolation_Order::LINEAR );
                        REQUIRE( tIntersectionQueued == tIsEdgeIntersected( tGeometryIndex )( tElementIndex )( tNodeNumber ) );

                        // Check queued intersection
                        if ( tIntersectionQueued )
                        {
                            // Check parents
                            bool tFirstParentOnInterface  = false;
                            bool tSecondParentOnInterface = false;

                            // Parent nodes on interface
                            if ( tGeometryIndex == 2 and tFirstNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tFirstParentOnInterface = true;
                            }
                            if ( tGeometryIndex == 2 and tSecondNodeCoordinates( 0 ) == Approx( 1.0 ) )
                            {
                                tSecondParentOnInterface = true;
                            }
                            CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() == tFirstParentOnInterface );
                            CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() == tSecondParentOnInterface );

                            // See if local coordinate is a number
                            real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();
                            if ( not std::isnan( tLocalCoordinate ) )
                            {
                                // Check local coordinate
                                CHECK( tLocalCoordinate == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ) );

                                // Check global coordinates
                                CHECK_EQUAL(
                                        tGeometryEngine.get_queued_intersection_global_coordinates(),
                                        tIntersectionGlobalCoordinates( tIntersectionCount ), );
                            }

                            // Admit intersection
                            tGeometryEngine.admit_queued_intersection();
                            tIntersectionCount++;
                        }

                        // Set node coordinates for element checking
                        tNodeCoordinates.set_row( tNodeNumber, tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) ) );
                    }

                    // Check with solution
                    CHECK( tGeometryEngine.is_intersected( tGeometryIndex, tSignedNodeIndices, tNodeCoordinates ) == tIsElementIntersected( tGeometryIndex )( tElementIndex ) );
                }

                // Advance geometry index
                if ( tGeometryIndex < 2 )
                {
                    tGeometryEngine.advance_geometry_index();
                }
            }

            // Check total number of intersections
            CHECK( tIntersectionCount == 10 );

            // Test the new child nodes on the level set field (geometry 0)
            for ( uint iNodeIndex = 9; iNodeIndex < 19; iNodeIndex++ )
            {
                CHECK( tGeometryEngine.get_geometric_region( 0, iNodeIndex, { {} } ) == Geometric_Region::POSITIVE );
            }

            // Clean up
            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Bilinear Intersections", "[gen], [pdv], [intersection], [bilinear intersection]" )
    {
        if ( par_size() == 1 )
        {
            // Create mesh
            mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 2, 2 );

            // Set up circle
            Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            tFieldParameterLists.set( "center_x", -0.25 );
            tFieldParameterLists.set( "center_y", 0.0 );
            tFieldParameterLists.set( "radius", 0.7499999999 );
            tFieldParameterLists.set( "discretization_mesh_index", 0 );
            tFieldParameterLists.set( "use_multilinear_interpolation", true );

            // Create geometry engine
            Geometry_Engine_Parameters tGeometryEngineParameters;
            ADV_Manager tADVManager;
            Design_Factory tDesignFactory( tFieldParameterLists, tADVManager );
            tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
            Geometry_Engine_Test tGeometryEngine( tMesh, tGeometryEngineParameters );

            // Solution for is_intersected() per geometry and per element
            Vector< bool > tIsElementIntersected = { true, true, true, true };

            // Per element, per edge
            Vector< Vector< bool > > tIsEdgeIntersected = {
                { false, true, true, false },    // Element 0
                { false, false, true, true },    // Element 1
                { true, true, false, false },    // Element 2
                { true, false, false, true }
            };    // Element 3

            // Intersection coordinates
            real tFrac = 2.0 / ( 3.0 + sqrt( 17.0 ) );

            Matrix< DDRMat > tIntersectionLocalCoordinates = {
                { -tFrac, 1.0, 0.0, tFrac, -1.0, tFrac, 0.0, -tFrac, -0.5, 0.5, -0.5, 0.5 }
            };

            Vector< Matrix< DDRMat > > tIntersectionGlobalCoordinates = {
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.5, 0.0 } },
                { { 0.0, -0.5 - ( tFrac / 2.0 ) } },
                { { -1.0, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } },
                { { 0.5, 0.0 } },
                { { 0.0, 0.5 + ( tFrac / 2.0 ) } }
            };

            // Initialize sensitivity variables
            real             tEpsilon = 1E-12;
            Matrix< DDRMat > tHostADVSensitivities;
            Matrix< DDRMat > tI;
            eye( 2, 2, tI );

            // Check element intersections
            uint tIntersectionCount = 0;
            for ( uint tElementIndex = 0; tElementIndex < 4; tElementIndex++ )
            {
                // Get element info
                Matrix< IndexMat > tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( tElementIndex );
                Matrix< DDUMat >   tNodeIndices( 4, 1 );
                for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                {
                    tNodeIndices( tNodeNumber ) = tSignedNodeIndices( tNodeNumber );
                }

                // Check edges for properly queued intersections
                for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
                {
                    // Queue intersection
                    bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                            tNodeIndices( tNodeNumber ),
                            tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                            get_quad_local_coordinates( tNodeNumber ),
                            get_quad_local_coordinates( ( tNodeNumber + 1 ) % 4 ),
                            tNodeIndices,
                            mtk::Geometry_Type::QUAD,
                            mtk::Interpolation_Order::LINEAR );
                    REQUIRE( tIntersectionQueued == tIsEdgeIntersected( tElementIndex )( tNodeNumber ) );

                    // Check queued intersection
                    if ( tIntersectionQueued )
                    {
                        // Check parents
                        CHECK( not tGeometryEngine.queued_intersection_first_parent_on_interface() );
                        CHECK( not tGeometryEngine.queued_intersection_second_parent_on_interface() );

                        // Check local coordinates
                        CHECK( tGeometryEngine.get_queued_intersection_local_coordinate() == Approx( tIntersectionLocalCoordinates( tIntersectionCount ) ).margin( 1e-9 ) );

                        // Check global coordinates
                        CHECK_EQUAL( tGeometryEngine.get_queued_intersection_global_coordinates(), tIntersectionGlobalCoordinates( tIntersectionCount ), );

                        // Admit intersection
                        tGeometryEngine.admit_queued_intersection();

                        // Get node manager
                        Node_Manager& tNodeManager = tGeometryEngine.get_node_manager();

                        // Check sensitivities
                        tHostADVSensitivities.set_size( 0.0, 0.0 );
                        tNodeManager.append_dcoordinate_dadv_from_derived_node( 9 + tIntersectionCount, tHostADVSensitivities, tI );
                        Vector< sint > tADVIDs = tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( 9 + tIntersectionCount );

                        // Finite difference sensitivities by queueing dummy nodes
                        Matrix< DDRMat > tFDSensitivities( tHostADVSensitivities.n_rows(), tHostADVSensitivities.n_cols(), 0.0 );
                        Vector< real > tADVs = tGeometryEngine.get_advs();
                        for ( uint iADVIndex = 0; iADVIndex < tADVIDs.size(); iADVIndex++ )
                        {
                            // Get ADV ID
                            sint tADVID = tADVIDs( iADVIndex );

                            // Positive perturbation
                            tADVs( tADVID - 1 ) += tEpsilon;
                            tGeometryEngine.set_advs( tADVs );
                            tGeometryEngine.queue_intersection(
                                    tNodeIndices( tNodeNumber ),
                                    tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                    get_quad_local_coordinates( tNodeNumber ),
                                    get_quad_local_coordinates( ( tNodeNumber + 1 ) % 4 ),
                                    tNodeIndices,
                                    mtk::Geometry_Type::QUAD,
                                    mtk::Interpolation_Order::LINEAR );
                            Matrix< DDRMat > tPositiveGlobalCoordinates = tGeometryEngine.get_queued_intersection_global_coordinates();

                            // Negative perturbation
                            tADVs( tADVID - 1 ) -= 2.0 * tEpsilon;
                            tGeometryEngine.set_advs( tADVs );
                            tGeometryEngine.queue_intersection(
                                    tNodeIndices( tNodeNumber ),
                                    tNodeIndices( ( tNodeNumber + 1 ) % 4 ),
                                    get_quad_local_coordinates( tNodeNumber ),
                                    get_quad_local_coordinates( ( tNodeNumber + 1 ) % 4 ),
                                    tNodeIndices,
                                    mtk::Geometry_Type::QUAD,
                                    mtk::Interpolation_Order::LINEAR );
                            Matrix< DDRMat > tNegativeGlobalCoordinates = tGeometryEngine.get_queued_intersection_global_coordinates();

                            // Reset
                            tADVs( tADVID - 1 ) += tEpsilon;
                            tGeometryEngine.set_advs( tADVs );
                            tFDSensitivities( 0, iADVIndex ) = ( tPositiveGlobalCoordinates( 0 ) - tNegativeGlobalCoordinates( 0 ) ) / ( 2.0 * tEpsilon );
                            tFDSensitivities( 1, iADVIndex ) = ( tPositiveGlobalCoordinates( 1 ) - tNegativeGlobalCoordinates( 1 ) ) / ( 2.0 * tEpsilon );
                        }

                        // Check sensitivities, needs a generous error factor
                        CHECK_EQUAL( tHostADVSensitivities, tFDSensitivities, 1E12, );

                        // Increment intersection count
                        tIntersectionCount++;
                    }
                }
            }

            // Check total number of intersections
            CHECK( tIntersectionCount == 8 );

            // Test the new child nodes on the level set field
            CHECK( tGeometryEngine.get_geometric_region( 0, 9, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 10, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 11, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 12, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 13, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 14, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 15, { {} } ) == Geometric_Region::INTERFACE );
            CHECK( tGeometryEngine.get_geometric_region( 0, 16, { {} } ) == Geometric_Region::INTERFACE );

            // Get full element info for element 0
            Matrix< IndexMat >       tSignedNodeIndices = tMesh->get_nodes_connected_to_element_loc_inds( 0 );
            Matrix< DDUMat >         tNodeIndices( 4, 1 );
            Vector< Matrix< DDRMat > > tNodeCoordinates( 4 );
            for ( uint tNodeNumber = 0; tNodeNumber < 4; tNodeNumber++ )
            {
                tNodeIndices( tNodeNumber )     = tSignedNodeIndices( tNodeNumber );
                tNodeCoordinates( tNodeNumber ) = tMesh->get_node_coordinate( tNodeIndices( tNodeNumber ) );
            }

            // Queue custom intersection 1 and check for bilinear intersection
            bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                    0, 2, { { -1.0, -1.0 } }, { { 1.0, 1.0 } }, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
            REQUIRE( tIntersectionQueued );

            // Queue custom intersection 2 and check for no bilinear intersection
            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    1, 3, { { 1.0, -1.0 } }, { { -1.0, 1.0 } }, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
            REQUIRE( not tIntersectionQueued );

            // Queue custom intersection 3 and check for bilinear intersection
            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    9, 10, { { 1.0, tFrac } }, { { -1.0, 1.0 } }, tNodeIndices, mtk::Geometry_Type::QUAD, mtk::Interpolation_Order::LINEAR );
            REQUIRE( tIntersectionQueued );

            // Clean up
            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Bilinear Intersections with Nonzero Threshold", "[gen], [pdv], [intersection], [bilinear intersection], [nonzero_bilinear_threshold]" )
    {
        if ( par_size() == 1 )
        {
            // This tests a pathological case where the curved bi-linear interpolation results in an edge with the same signed parent nodes to have an intersection point
            // My mesh is one element where the geometry is a user defined nodal field (best way to control everything in this problem).
            // The isocontour threshold is away from 0. If this test fails, it is probably because someone modified the logic in cl_GEN_Intersection_Node.cpp which determins
            // if an entity is intersected.

            // Generate mesh
            uint             aNumDim = 2;
            Matrix< DDRMat > aCoords( 6, 2 );
            aCoords( 0, 0 ) = 0.0, aCoords( 0, 1 ) = 0.0;
            aCoords( 1, 0 ) = 0.25, aCoords( 1, 1 ) = 0.0;
            aCoords( 2, 0 ) = 0.25, aCoords( 2, 1 ) = 0.25;
            aCoords( 3, 0 ) = 0.0, aCoords( 3, 1 ) = 1.0;
            aCoords( 4, 0 ) = 0.0, aCoords( 4, 1 ) = 0.0;
            Matrix< IdMat > aElemConn( 1, 4 );

            // 0D to 3D connectivity (node to element)
            aElemConn( 0, 0 ) = 1;
            aElemConn( 0, 1 ) = 2;
            aElemConn( 0, 2 ) = 3;
            aElemConn( 0, 3 ) = 4;

            Matrix< IdMat > aElemLocaltoGlobal = { { 1 } };

            // Create MORIS mesh using MTK database
            moris::mtk::MtkMeshData aMeshData;
            aMeshData.CreateAllEdgesAndFaces    = true;
            aMeshData.SpatialDim                = &aNumDim;
            aMeshData.ElemConn( 0 )             = &aElemConn;
            aMeshData.NodeCoords                = &aCoords;
            aMeshData.LocaltoGlobalElemMap( 0 ) = &aElemLocaltoGlobal;
            moris::mtk::Scalar_Field_Info< DDRMat > tLSF;
            std::string                             tLSFName = "lsf1";
            tLSF.set_field_name( tLSFName );
            tLSF.set_field_entity_rank( mtk::EntityRank::NODE );

            // Add to mesh input field container
            moris::mtk::MtkFieldsInfo tFieldsInfo;
            add_field_for_mesh_input( &tLSF, tFieldsInfo );
            aMeshData.FieldsInfo = &tFieldsInfo;

            moris::mtk::Interpolation_Mesh* tMeshData = moris::mtk::create_interpolation_mesh( mtk::MeshType::STK, aMeshData );

            moris::uint                    tNumNodes = tMeshData->get_num_entities( mtk::EntityRank::NODE );
            moris::Matrix< moris::DDRMat > tLevelsetVal( tNumNodes, 1, -1.3 );

            moris_id tIndexOfNodeId1 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 1, mtk::EntityRank::NODE );
            moris_id tIndexOfNodeId2 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 2, mtk::EntityRank::NODE );
            moris_id tIndexOfNodeId3 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 3, mtk::EntityRank::NODE );
            moris_id tIndexOfNodeId4 = tMeshData->get_loc_entity_ind_from_entity_glb_id( 4, mtk::EntityRank::NODE );

            tLevelsetVal( tIndexOfNodeId1 ) = 5.71398452074828e-01;
            tLevelsetVal( tIndexOfNodeId2 ) = 5.06776630434012e-01;
            tLevelsetVal( tIndexOfNodeId3 ) = 4.25145951405591e-01;
            tLevelsetVal( tIndexOfNodeId4 ) = 5.00394812706282e-01;

            tMeshData->add_mesh_field_real_scalar_data_loc_inds( tLSFName, mtk::EntityRank::NODE, tLevelsetVal );

            Vector< std::shared_ptr< gen::Geometry > > tGeometry( 1 );
            Level_Set_Parameters                     tLevelSetParameters;
            tLevelSetParameters.mUseMultilinearInterpolation = true;
            tLevelSetParameters.mIsocontourThreshold         = 0.5;
            tLevelSetParameters.mIsocontourTolerance         = 1E-13;
            auto tField                                      = std::make_shared< Mesh_Field >( tMeshData, tLSFName );
            tGeometry( 0 )                                   = std::make_shared< Level_Set_Geometry >( tField, tLevelSetParameters );

            Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;

            Geometry_Engine tGeometryEngine( tMeshData, tGeometryEngineParameters );

            // get the cell
            moris::mtk::Cell& tCell = tMeshData->get_mtk_cell( 0 );

            // verify that it says that the cell is intersected
            bool tIsIntersected = tGeometryEngine.is_intersected( 0, tCell.get_vertex_inds(), tCell.get_vertex_coords() );

            CHECK( tIsIntersected );

            Vector< Matrix< IndexMat > > tVertexIndices = { tCell.get_vertex_inds() };

            Vector< Matrix< DDRMat > > tLocalCoords( 1 );

            tLocalCoords( 0 ) = { { +0.000000000000000e+00, +0.000000000000000e+00 } };

            tGeometryEngine.create_new_derived_nodes(
                    tVertexIndices,
                    tLocalCoords,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            Matrix< DDRMat >   tVertexCoords  = tCell.get_vertex_coords();
            Matrix< IndexMat > tVertexInds    = tCell.get_vertex_inds();
            Matrix< DDUMat >   tVertexIndsDDU = { { (uint)tVertexInds( 0 ), (uint)tVertexInds( 1 ), (uint)tVertexInds( 2 ), (uint)tVertexInds( 3 ) } };

            Vector< Matrix< DDRMat > > tBGCellCoords( 4 );
            tBGCellCoords( 0 ) = tVertexCoords.get_row( 0 );
            tBGCellCoords( 1 ) = tVertexCoords.get_row( 1 );
            tBGCellCoords( 2 ) = tVertexCoords.get_row( 2 );
            tBGCellCoords( 3 ) = tVertexCoords.get_row( 3 );

            Matrix< DDRMat > tLocalCoordsMat = {
                { -1.000000000000000e+00, -1.000000000000000e+00 },
                { +1.000000000000000e+00, -1.000000000000000e+00 },
                { +1.000000000000000e+00, +1.000000000000000e+00 },
                { -1.000000000000000e+00, +1.000000000000000e+00 },
                { +0.000000000000000e+00, +0.000000000000000e+00 }
            };

            // check that the cell is intersected
            moris_index tNodeIndex1 = 0;
            moris_index tNodeIndex2 = 1;

            bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 1;
            tNodeIndex2 = 2;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( tIntersectionQueued );

            tNodeIndex1 = 2;
            tNodeIndex2 = 3;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );
            CHECK( tIntersectionQueued );

            tNodeIndex1 = 0;
            tNodeIndex2 = 3;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );
            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 0;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 1;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( !tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 2;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );

            CHECK( tIntersectionQueued );

            // check that the cell is intersected
            tNodeIndex1 = 3;
            tNodeIndex2 = 4;

            tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tNodeIndex1,
                    tNodeIndex2,
                    tLocalCoordsMat.get_row( tNodeIndex1 ),
                    tLocalCoordsMat.get_row( tNodeIndex2 ),
                    tVertexIndsDDU,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );
            CHECK( !tIntersectionQueued );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Bilinear Intersections with Zero Corner Values", "[gen], [pdv], [intersection], [bilinear intersection], [MB-020]" )
    {
        if ( par_size() == 1 )
        {
            // Regression test for MB-020: the multilinear intersection scales the corner level-set values by the
            // inverse of their norm. If every corner value of the background cell is exactly zero (the multilinear
            // field vanishes on the whole cell, e.g. an order-2 field that is nonzero only at mid-side and interior
            // nodes), the scaling was 1/0 = inf, every scaled value became NaN, and the Newton search aborted with
            // "Newton did not converge: Reference residual nan". The local coordinate must stay finite instead.

            // One QUAD4 element on the unit square
            uint             tNumDim = 2;
            Matrix< DDRMat > tCoords = {
                { 0.0, 0.0 },
                { 1.0, 0.0 },
                { 1.0, 1.0 },
                { 0.0, 1.0 }
            };
            Matrix< IdMat > tElemConn          = { { 1, 2, 3, 4 } };
            Matrix< IdMat > tElemLocaltoGlobal = { { 1 } };

            // Nodal level-set field
            mtk::Scalar_Field_Info< DDRMat > tLSF;
            std::string                      tLSFName = "lsf_zero";
            tLSF.set_field_name( tLSFName );
            tLSF.set_field_entity_rank( mtk::EntityRank::NODE );
            mtk::MtkFieldsInfo tFieldsInfo;
            add_field_for_mesh_input( &tLSF, tFieldsInfo );

            // Create mesh
            mtk::MtkMeshData tMeshData;
            tMeshData.CreateAllEdgesAndFaces    = true;
            tMeshData.SpatialDim                = &tNumDim;
            tMeshData.ElemConn( 0 )             = &tElemConn;
            tMeshData.NodeCoords                = &tCoords;
            tMeshData.LocaltoGlobalElemMap( 0 ) = &tElemLocaltoGlobal;
            tMeshData.FieldsInfo                = &tFieldsInfo;
            mtk::Interpolation_Mesh* tMesh      = mtk::create_interpolation_mesh( mtk::MeshType::STK, tMeshData );

            // Every corner value exactly zero
            uint             tNumNodes = tMesh->get_num_entities( mtk::EntityRank::NODE );
            Matrix< DDRMat > tLevelSetValues( tNumNodes, 1, 0.0 );
            tMesh->add_mesh_field_real_scalar_data_loc_inds( tLSFName, mtk::EntityRank::NODE, tLevelSetValues );

            // Multilinear level-set geometry on the zero isocontour
            Level_Set_Parameters tLevelSetParameters;
            tLevelSetParameters.mUseMultilinearInterpolation  = true;
            tLevelSetParameters.mIsocontourThreshold          = 0.0;
            std::shared_ptr< Mesh_Field >         tField      = std::make_shared< Mesh_Field >( tMesh, tLSFName );
            Vector< std::shared_ptr< Geometry > > tGeometries = { std::make_shared< Level_Set_Geometry >( tField, tLevelSetParameters ) };

            Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometries;
            Geometry_Engine tGeometryEngine( tMesh, tGeometryEngineParameters );

            // Background cell
            mtk::Cell&         tCell          = tMesh->get_mtk_cell( 0 );
            Matrix< IndexMat > tVertexInds    = tCell.get_vertex_inds();
            Matrix< DDUMat >   tVertexIndsDDU = { { (uint)tVertexInds( 0 ), (uint)tVertexInds( 1 ), (uint)tVertexInds( 2 ), (uint)tVertexInds( 3 ) } };

            // Queue an intersection on every edge and both diagonals: each one calls
            // Level_Set_Geometry::compute_intersection_local_coordinate() on the all-zero cell
            Matrix< DDUMat > tParentPairs = {
                { 0, 1 },
                { 1, 2 },
                { 2, 3 },
                { 3, 0 },
                { 0, 2 },
                { 1, 3 }
            };
            for ( uint iPair = 0; iPair < tParentPairs.n_rows(); iPair++ )
            {
                uint tFirstNode  = tParentPairs( iPair, 0 );
                uint tSecondNode = tParentPairs( iPair, 1 );

                REQUIRE_NOTHROW( tGeometryEngine.queue_intersection(
                        tVertexIndsDDU( tFirstNode ),
                        tVertexIndsDDU( tSecondNode ),
                        tQuadParametricCoordinates( tFirstNode ),
                        tQuadParametricCoordinates( tSecondNode ),
                        tVertexIndsDDU,
                        mtk::Geometry_Type::QUAD,
                        mtk::Interpolation_Order::LINEAR ) );

                // Local coordinate along the parent edge is finite and in [-1, 1]
                real tLocalCoordinate = tGeometryEngine.get_queued_intersection_local_coordinate();
                CHECK( std::isfinite( tLocalCoordinate ) );
                CHECK( tLocalCoordinate >= -1.0 );
                CHECK( tLocalCoordinate <= 1.0 );

                // Both parents lie on the zero isocontour
                CHECK( tGeometryEngine.queued_intersection_first_parent_on_interface() );
                CHECK( tGeometryEngine.queued_intersection_second_parent_on_interface() );
            }

            // Clean up
            delete tMesh;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Builds a geometry engine with ADV-driven lines on a single [-1, 1]^2 quad and queues a chain of (nested)
     * intersections, then returns the coordinates, local coordinates and dense coordinate sensitivities of every
     * derived node. Each line is phi = n_x ( x - c_x ) + n_y ( y - c_y ) with all four parameters as ADVs, so the ADV
     * vector is { c_x, c_y, n_x, n_y } for line 0, then line 1, and so on.
     *
     * @param aMesh Single-element [-1, 1]^2 interpolation mesh
     * @param aLineParameters Line parameters, 4 per line; these are the ADVs
     * @param aEdges Queue sequence { first parent, second parent, cutting line }, in increasing line order. A parent
     * below 4 is a corner of the background element (0: (-1, -1), 1: (1, -1), 2: (1, 1), 3: (-1, 1)); a parent of 4 or
     * more is the derived node with that index (edge i creates node 4 + i).
     * @param aCoordinates Global coordinates of each derived node
     * @param aLocalCoordinates Local coordinate of each derived node on its parent edge
     * @param aSensitivities Dense dX/dADV ( 2 x number of ADVs ) of each derived node
     */
    static void
    build_nested_line_intersections(
            mtk::Interpolation_Mesh*        aMesh,
            const Vector< real >&           aLineParameters,
            const Vector< Vector< uint > >& aEdges,
            Vector< Matrix< DDRMat > >&     aCoordinates,
            Vector< real >&                 aLocalCoordinates,
            Vector< Matrix< DDRMat > >&     aSensitivities )
    {
        // Lines, every parameter an ADV (lower bound = value = upper bound)
        uint                      tNumLines = aLineParameters.size() / 4;
        Submodule_Parameter_Lists tFieldParameterLists( "FIELDS" );
        for ( uint iLine = 0; iLine < tNumLines; iLine++ )
        {
            real tCenterX = aLineParameters( 4 * iLine );
            real tCenterY = aLineParameters( 4 * iLine + 1 );
            real tNormalX = aLineParameters( 4 * iLine + 2 );
            real tNormalY = aLineParameters( 4 * iLine + 3 );
            tFieldParameterLists.add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tFieldParameterLists.set( "center_x", tCenterX, tCenterX, tCenterX );
            tFieldParameterLists.set( "center_y", tCenterY, tCenterY, tCenterY );
            tFieldParameterLists.set( "normal_x", tNormalX, tNormalX, tNormalX );
            tFieldParameterLists.set( "normal_y", tNormalY, tNormalY, tNormalY );
        }

        // Create geometry engine
        Geometry_Engine_Parameters tGeometryEngineParameters;
        ADV_Manager                tADVManager;
        Design_Factory             tDesignFactory( tFieldParameterLists, tADVManager );
        tGeometryEngineParameters.mADVManager = tADVManager;
        tGeometryEngineParameters.mGeometries = tDesignFactory.get_geometries();
        Geometry_Engine_Test tGeometryEngine( aMesh, tGeometryEngineParameters );

        // The ADV vector is ordered line by line as { c_x, c_y, n_x, n_y }
        Vector< real >& tADVs = tGeometryEngine.get_advs();
        REQUIRE( tADVs.size() == aLineParameters.size() );
        for ( uint iADV = 0; iADV < tADVs.size(); iADV++ )
        {
            REQUIRE( tADVs( iADV ) == aLineParameters( iADV ) );
        }

        // Background element; parametric and global coordinates coincide on the single [-1, 1]^2 element
        Matrix< IndexMat > tSignedNodeIndices = aMesh->get_nodes_connected_to_element_loc_inds( 0 );
        Matrix< DDUMat >   tElementNodes( 4, 1 );
        for ( uint iNode = 0; iNode < 4; iNode++ )
        {
            tElementNodes( iNode )            = tSignedNodeIndices( iNode );
            Matrix< DDRMat > tNodeCoordinates = aMesh->get_node_coordinate( tElementNodes( iNode ) );
            REQUIRE( tNodeCoordinates( 0 ) == Approx( tQuadParametricCoordinates( iNode )( 0 ) ) );
            REQUIRE( tNodeCoordinates( 1 ) == Approx( tQuadParametricCoordinates( iNode )( 1 ) ) );
        }
        Node_Manager& tNodeManager = tGeometryEngine.get_node_manager();
        REQUIRE( tNodeManager.get_total_number_of_nodes() == 4 );

        uint tActiveLine = 0;
        aLocalCoordinates.resize( aEdges.size() );
        for ( uint iEdge = 0; iEdge < aEdges.size(); iEdge++ )
        {
            // Advance to the line that cuts this edge
            REQUIRE( aEdges( iEdge )( 2 ) >= tActiveLine );
            REQUIRE( aEdges( iEdge )( 2 ) < tNumLines );
            while ( tActiveLine < aEdges( iEdge )( 2 ) )
            {
                tGeometryEngine.advance_geometry_index();
                tActiveLine++;
            }

            // Parent node indices: element corners or derived nodes
            uint tFirstNodeIndex  = aEdges( iEdge )( 0 ) < 4 ? tElementNodes( aEdges( iEdge )( 0 ) ) : aEdges( iEdge )( 0 );
            uint tSecondNodeIndex = aEdges( iEdge )( 1 ) < 4 ? tElementNodes( aEdges( iEdge )( 1 ) ) : aEdges( iEdge )( 1 );

            // Parametric coordinates of the edge nodes (equal to the global coordinates on this element)
            const Matrix< DDRMat >& tFirstCoordinates  = tNodeManager.get_node( tFirstNodeIndex ).get_global_coordinates();
            const Matrix< DDRMat >& tSecondCoordinates = tNodeManager.get_node( tSecondNodeIndex ).get_global_coordinates();
            Matrix< DDRMat >        tFirstParametric   = { { tFirstCoordinates( 0 ), tFirstCoordinates( 1 ) } };
            Matrix< DDRMat >        tSecondParametric  = { { tSecondCoordinates( 0 ), tSecondCoordinates( 1 ) } };

            // Fixed branch: every edge is cut in its interior, never at a parent
            bool tIntersectionQueued = tGeometryEngine.queue_intersection(
                    tFirstNodeIndex,
                    tSecondNodeIndex,
                    tFirstParametric,
                    tSecondParametric,
                    tElementNodes,
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Order::LINEAR );
            REQUIRE( tIntersectionQueued );
            REQUIRE( not tGeometryEngine.queued_intersection_first_parent_on_interface() );
            REQUIRE( not tGeometryEngine.queued_intersection_second_parent_on_interface() );
            aLocalCoordinates( iEdge ) = tGeometryEngine.get_queued_intersection_local_coordinate();
            REQUIRE( std::abs( aLocalCoordinates( iEdge ) ) < 0.9 );

            // Admit; the new node takes the next index
            tGeometryEngine.admit_queued_intersection();
            REQUIRE( tNodeManager.get_total_number_of_nodes() == 5 + iEdge );
        }

        // Coordinates and dense dX/dADV of every derived node
        aCoordinates.resize( aEdges.size() );
        aSensitivities.resize( aEdges.size() );
        for ( uint iEdge = 0; iEdge < aEdges.size(); iEdge++ )
        {
            uint tNodeIndex       = 4 + iEdge;
            aCoordinates( iEdge ) = tNodeManager.get_node( tNodeIndex ).get_global_coordinates();

            Matrix< DDRMat > tI;
            eye( 2, 2, tI );
            Matrix< DDRMat > tSensitivities;
            tNodeManager.append_dcoordinate_dadv_from_derived_node( tNodeIndex, tSensitivities, tI );
            Vector< sint > tADVIDs = tNodeManager.get_coordinate_determining_adv_ids_from_derived_node( tNodeIndex );
            REQUIRE( tSensitivities.n_rows() == 2 );
            REQUIRE( tSensitivities.n_cols() == tADVIDs.size() );

            // Accumulate repeated ADV IDs into a dense matrix
            aSensitivities( iEdge ).set_size( 2, aLineParameters.size(), 0.0 );
            for ( uint iColumn = 0; iColumn < tADVIDs.size(); iColumn++ )
            {
                REQUIRE( tADVIDs( iColumn ) >= 0 );
                REQUIRE( tADVIDs( iColumn ) < static_cast< sint >( aLineParameters.size() ) );
                for ( uint iDimension = 0; iDimension < 2; iDimension++ )
                {
                    aSensitivities( iEdge )( iDimension, tADVIDs( iColumn ) ) += tSensitivities( iDimension, iColumn );
                }
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    /**
     * Checks the analytic coordinate sensitivities of every derived node of a nested line-intersection chain against
     * central finite differences of the intersection coordinates. Every perturbed state rebuilds the whole chain in a
     * fresh geometry engine and must stay on the nominal branch (same parent edges, nearby local coordinates).
     *
     * @param aLineParameters Line parameters, 4 per line; these are the ADVs
     * @param aEdges Queue sequence, see build_nested_line_intersections()
     * @param aDependsOnLine Per derived node and line, 1 if the node coordinates must depend on that line's ADVs
     * @return Maximum relative error per derived node
     */
    static Vector< real >
    check_nested_line_intersections_against_finite_differences(
            const Vector< real >&           aLineParameters,
            const Vector< Vector< uint > >& aEdges,
            const Vector< Vector< uint > >& aDependsOnLine )
    {
        // Single [-1, 1]^2 element
        mtk::Interpolation_Mesh* tMesh = create_simple_mesh( 1, 1 );

        uint tNumADVs  = aLineParameters.size();
        uint tNumLines = tNumADVs / 4;
        uint tNumNodes = aEdges.size();

        // Analytic sensitivities at the nominal design
        Vector< Matrix< DDRMat > > tCoordinates;
        Vector< real >             tLocalCoordinates;
        Vector< Matrix< DDRMat > > tAnalytic;
        build_nested_line_intersections( tMesh, aLineParameters, aEdges, tCoordinates, tLocalCoordinates, tAnalytic );

        // Central finite differences of the intersection coordinates
        real                       tStep = 1.0e-6;
        Vector< Matrix< DDRMat > > tFiniteDifference( tNumNodes );
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            tFiniteDifference( iNode ).set_size( 2, tNumADVs, 0.0 );
        }
        for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
        {
            Vector< real > tPlusADVs  = aLineParameters;
            Vector< real > tMinusADVs = aLineParameters;
            tPlusADVs( iADV ) += tStep;
            tMinusADVs( iADV ) -= tStep;

            Vector< Matrix< DDRMat > > tPlusCoordinates;
            Vector< Matrix< DDRMat > > tMinusCoordinates;
            Vector< real >             tPlusLocalCoordinates;
            Vector< real >             tMinusLocalCoordinates;
            Vector< Matrix< DDRMat > > tUnusedSensitivities;
            build_nested_line_intersections( tMesh, tPlusADVs, aEdges, tPlusCoordinates, tPlusLocalCoordinates, tUnusedSensitivities );
            build_nested_line_intersections( tMesh, tMinusADVs, aEdges, tMinusCoordinates, tMinusLocalCoordinates, tUnusedSensitivities );

            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                // Same branch: the perturbed intersections stay next to the nominal ones on the same edges
                REQUIRE( std::abs( tPlusLocalCoordinates( iNode ) - tLocalCoordinates( iNode ) ) < 1.0e-4 );
                REQUIRE( std::abs( tMinusLocalCoordinates( iNode ) - tLocalCoordinates( iNode ) ) < 1.0e-4 );

                for ( uint iDimension = 0; iDimension < 2; iDimension++ )
                {
                    tFiniteDifference( iNode )( iDimension, iADV ) =
                            ( tPlusCoordinates( iNode )( iDimension ) - tMinusCoordinates( iNode )( iDimension ) ) / ( 2.0 * tStep );
                }
            }
        }

        // Compare analytic to finite difference, entry by entry
        real           tTolerance = 1.0e-6;
        Vector< real > tMaxErrors( tNumNodes, 0.0 );
        for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
        {
            for ( uint iADV = 0; iADV < tNumADVs; iADV++ )
            {
                for ( uint iDimension = 0; iDimension < 2; iDimension++ )
                {
                    real tAnalyticValue = tAnalytic( iNode )( iDimension, iADV );
                    real tFDValue       = tFiniteDifference( iNode )( iDimension, iADV );
                    real tError         = std::abs( tAnalyticValue - tFDValue ) / std::max( 1.0, std::abs( tFDValue ) );
                    tMaxErrors( iNode ) = std::max( tMaxErrors( iNode ), tError );

                    INFO( "node " << 4 + iNode << ", dimension " << iDimension << ", ADV " << iADV
                                  << ": analytic " << tAnalyticValue << ", finite difference " << tFDValue );
                    CHECK( tError < tTolerance );
                }
            }

            // The chain is non-trivial: each node feels the ADVs of exactly the lines of itself and its ancestors
            for ( uint iLine = 0; iLine < tNumLines; iLine++ )
            {
                real tLineSum = 0.0;
                for ( uint iParameter = 0; iParameter < 4; iParameter++ )
                {
                    for ( uint iDimension = 0; iDimension < 2; iDimension++ )
                    {
                        tLineSum += std::abs( tFiniteDifference( iNode )( iDimension, 4 * iLine + iParameter ) );
                    }
                }
                INFO( "node " << 4 + iNode << ", line " << iLine << ": sum |dX/dADV| = " << tLineSum );
                CHECK( ( tLineSum > 1.0e-3 ) == ( aDependsOnLine( iNode )( iLine ) == 1 ) );
            }

            INFO( "node " << 4 + iNode << ": max relative error " << tMaxErrors( iNode ) );
            CHECK( tMaxErrors( iNode ) < tTolerance );
        }

        // Clean up
        delete tMesh;

        return tMaxErrors;
    }

    //--------------------------------------------------------------------------------------------------------------

    TEST_CASE( "Nested Linear Intersections Finite Difference", "[gen], [pdv], [intersection], [linear intersection], [MB-022]" )
    {
        if ( par_size() == 1 )
        {
            // Three levels of moving intersections, all on lines whose four parameters are ADVs:
            //   node 4 (P1): line 0 on corner 0 -> corner 1      level 1
            //   node 5 (P2): line 0 on corner 2 -> corner 3      level 1
            //   node 6 (Q):  line 1 on P1 -> P2                  level 2, both parents moving
            //   node 7 (R):  line 2 on Q -> corner 2             level 3, moving first parent
            // R reaches the ADVs of line 0 only through Q's parents, so those columns are a product of three chain
            // factors; dropping the incoming factor in the parent recursion (MB-022) breaks them.
            // Q lies on line 0 and moves sideways with it. Line 2's normal is parallel to Q's parent edge P1 -> P2
            // (both are (-0.2, 1)), so the edge-projected field gradient used for d(xi)/d(parent coordinates) is the
            // full gradient and MB-030 does not enter this check (see the [MB-030] case below for that defect).
            // Line parameters { c_x, c_y, n_x, n_y } per line (the ADVs):
            //   line 0: nearly vertical, cuts the bottom and top edges
            //   line 1: nearly horizontal, cuts P1 -> P2
            //   line 2: normal along P1 -> P2, cuts Q -> (1, 1)
            Vector< real >           tLineParameters = { 0.1, 0.0, 1.0, 0.2, 0.0, 0.05, 0.15, 1.0, 0.5, 0.5, -0.2, 1.0 };
            Vector< Vector< uint > > tEdges          = {
                { 0, 1, 0 },
                { 2, 3, 0 },
                { 4, 5, 1 },
                { 6, 2, 2 }
            };
            Vector< Vector< uint > > tDependsOnLine = {
                { 1, 0, 0 },
                { 1, 0, 0 },
                { 1, 1, 0 },
                { 1, 1, 1 }
            };
            check_nested_line_intersections_against_finite_differences( tLineParameters, tEdges, tDependsOnLine );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    // Hidden (known failure, MB-030): the field gradient used for d(xi)/d(parent coordinates) of an intersection-node
    // parent is the projection of grad(phi) onto that parent's own edge (Level_Set_Geometry::get_dfield_dcoordinates),
    // which is exact only while the parent moves along its edge. Here Q's parent edge P1 -> P2 lies on line 0 and moves
    // sideways with line 0, so R (level 3) is wrong in the line-0 columns and T (level 4) in the line-0 and line-1
    // columns. Run with GEN-test.exe "[MB-030]"; remove the [.] tag once MB-030 is fixed.
    TEST_CASE( "Nested Linear Intersections Finite Difference Off-Edge Parent", "[.], [gen], [pdv], [intersection], [MB-030]" )
    {
        if ( par_size() == 1 )
        {
            //   node 4 (P1): line 0 on corner 0 -> corner 1     level 1
            //   node 5 (P2): line 0 on corner 2 -> corner 3     level 1
            //   node 6 (Q):  line 1 on P1 -> P2                 level 2 (Q lies on line 0)
            //   node 7 (R):  line 2 on Q -> corner 2            level 3
            //   node 8 (T):  line 3 on P1 -> R                  level 4
            // Line parameters { c_x, c_y, n_x, n_y } per line (the ADVs):
            //   line 0: nearly vertical
            //   line 1: nearly horizontal
            //   line 2: cuts Q -> (1, 1)
            //   line 3: cuts P1 -> R
            Vector< real >           tLineParameters = { 0.1, 0.0, 1.0, 0.2, 0.0, 0.05, 0.15, 1.0, 0.6, 0.5, 1.0, 0.3, 0.5, -0.2, 1.0, -0.4 };
            Vector< Vector< uint > > tEdges          = {
                { 0, 1, 0 },
                { 2, 3, 0 },
                { 4, 5, 1 },
                { 6, 2, 2 },
                { 4, 7, 3 }
            };
            Vector< Vector< uint > > tDependsOnLine = {
                { 1, 0, 0, 0 },
                { 1, 0, 0, 0 },
                { 1, 1, 0, 0 },
                { 1, 1, 1, 0 },
                { 1, 1, 1, 1 }
            };
            check_nested_line_intersections_against_finite_differences( tLineParameters, tEdges, tDependsOnLine );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
