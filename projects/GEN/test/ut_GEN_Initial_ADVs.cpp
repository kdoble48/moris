/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * ut_GEN_Initial_ADVs.cpp
 *
 * MB-086: the GEN-driven refinement of the initial mesh runs before the ADVs are distributed. With an OPT restart file
 * the workflow hands the restart ADVs to GEN first (Geometry_Engine::set_initial_primitive_advs()), so the refinement and
 * the first ADV distribution both see the restart design.
 *
 */

#include "catch.hpp"

#include "cl_GEN_Geometry_Engine.hpp"
#include "fn_GEN_create_simple_mesh.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"

namespace moris::gen
{
    TEST_CASE( "MB-086 initial primitive ADVs", "[gen], [MB-086]" )
    {
        if ( par_size() == 1 )
        {
            // circle phi = |x - c| - r with the centre x-coordinate a design variable (initial value 0)
            Module_Parameter_Lists tParameterLists( Module_Type::GEN );
            tParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            tParameterLists.set( "center_x", Design_Variable( -1.0, 0.0, 1.0 ) );
            tParameterLists.set( "center_y", 0.0 );
            tParameterLists.set( "radius", 0.5 );
            tParameterLists.set( "number_of_refinements", 1 );
            tParameterLists.set( "refinement_mesh_index", 0 );

            Geometry_Engine tGeometryEngine( tParameterLists );

            REQUIRE( not tGeometryEngine.advs_distributed() );
            REQUIRE( tGeometryEngine.get_num_primitive_advs() == 1 );

            // the point (0.7, 0) lies outside the deck's circle (c = 0) ...
            Matrix< DDRMat > tPoint = { { 0.7, 0.0 } };
            CHECK( tGeometryEngine.get_geometric_region( 0, 0, tPoint ) == Geometric_Region::POSITIVE );

            // ... and inside the restart design's circle (c = 0.4): what the refinement of the initial mesh queries
            Vector< real > tRestartADVs = { 0.4 };
            CHECK( tGeometryEngine.set_initial_primitive_advs( tRestartADVs ) == Approx( 0.4 ) );
            CHECK( tGeometryEngine.get_geometric_region( 0, 0, tPoint ) == Geometric_Region::NEGATIVE );

            // the first distribution of the ADVs starts from the restart values, not from the deck's
            mtk::Interpolation_Mesh* tMesh            = create_simple_mesh( 4, 4, 1, 1 );
            mtk::Integration_Mesh*   tIntegrationMesh = mtk::create_integration_mesh_from_interpolation_mesh( mtk::MeshType::HMR, tMesh );
            mtk::Mesh_Pair           tMeshPair( tMesh, tIntegrationMesh, true );
            tGeometryEngine.distribute_advs( tMeshPair );

            REQUIRE( tGeometryEngine.advs_distributed() );
            Vector< real >& tADVs = tGeometryEngine.get_advs();
            REQUIRE( tADVs.size() == 1 );
            CHECK( tADVs( 0 ) == Approx( 0.4 ) );
            CHECK( tGeometryEngine.get_geometric_region( 0, 0, tPoint ) == Geometric_Region::NEGATIVE );
        }
    }

    TEST_CASE( "MB-089 feature re-detection parameter", "[gen], [MB-089]" )
    {
        if ( par_size() == 1 )
        {
            // three circles: the first two with sub-cell feature refinement, the last two with the re-detection key
            Module_Parameter_Lists tParameterLists( Module_Type::GEN );
            for ( uint iGeom = 0; iGeom < 3; iGeom++ )
            {
                tParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
                tParameterLists.set( "center_x", 0.0 );
                tParameterLists.set( "center_y", 0.0 );
                tParameterLists.set( "radius", 0.25 + 0.1 * iGeom );
                tParameterLists.set( "number_of_refinements", 1 );
                tParameterLists.set( "refinement_mesh_index", 0 );
                if ( iGeom < 2 )
                {
                    tParameterLists.set( "feature_refinement_samples", 4u );
                }
                if ( iGeom > 0 )
                {
                    tParameterLists.set( "feature_refinement_redetect_edges", 30 );
                }
            }

            Geometry_Engine tGeometryEngine( tParameterLists );

            // default off; on with feature refinement; ignored without it (nothing to redo)
            CHECK( tGeometryEngine.get_feature_refinement_redetect_edges( 0 ) == -1 );
            CHECK( tGeometryEngine.get_feature_refinement_redetect_edges( 1 ) == 30 );
            CHECK( tGeometryEngine.get_feature_refinement_redetect_edges( 2 ) == -1 );
        }
    }
}    // namespace moris::gen
