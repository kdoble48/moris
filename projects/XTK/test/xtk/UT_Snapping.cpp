/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_Snapping.cpp
 *
 */

#include "catch.hpp"
#include <string>
#include "cl_Vector.hpp"

// XTKL: Geometry  Include
#include "cl_Logger.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_HMR.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "fn_PRM_FEM_Parameters.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_GEN_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_Module_Parameter_Lists.hpp"

namespace moris::xtk
{

    TEST_CASE( "Plane Aligned with Background 2d", "[Snapping_2d]" )
    {
        if ( par_size() == 1 )
        {
            // Intersecting plane test parameters
            real tPlaneAngle                  = 90;      // CCW, degrees
            real tIntersectionLocalCoordinate = -0.5;    // (-1, 1), where the top middle vertical edge is intersected
            // bool tSnapFromRight = true;
            // bool tSnapBoth = false;
            bool tBilinear = false;

            // XTK parameter list
            Module_Parameter_Lists tXTKParams( Module_Type::XTK );
            tXTKParams( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
            tXTKParams( 0 )( 0 ).set( "decompose", true );
            tXTKParams( 0 )( 0 ).set( "decomposition_type", "conformal" );
            tXTKParams( 0 )( 0 ).set( "enrich", true );
            tXTKParams( 0 )( 0 ).set( "basis_rank", "bspline" );
            tXTKParams( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
            tXTKParams( 0 )( 0 ).set( "ghost_stab", true );
            tXTKParams( 0 )( 0 ).set( "multigrid", false );
            tXTKParams( 0 )( 0 ).set( "verbose", true );
            tXTKParams( 0 )( 0 ).set( "print_enriched_ig_mesh", true );
            tXTKParams( 0 )( 0 ).set( "exodus_output_XTK_ig_mesh", true );
            tXTKParams( 0 )( 0 ).set( "high_to_low_dbl_side_sets", true );
            tXTKParams( 0 )( 0 ).set( "probe_bg_cells", "1,2" );

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // GEN Parameters
            Module_Parameter_Lists tGENParams( Module_Type::GEN );

            // Calculations
            real tOffset  = 0.5 * ( 1 + tIntersectionLocalCoordinate );
            real tXNormal = sin( -tPlaneAngle * M_PI / 180 );
            real tYNormal = cos( -tPlaneAngle * M_PI / 180 );

            tGENParams( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
            tGENParams( 0 )( 0 ).set( "output_mesh_file", "GEN_Snapping.exo" );

            // Geometry parameter lists
            moris::uint tGeoCounter        = 0;
            tGENParams( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_threshold", 1e-16 );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "intersection_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "center_x", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "center_y", tOffset );
            tGENParams( 1 )( tGeoCounter ).set( "normal_x", tXNormal );
            tGENParams( 1 )( tGeoCounter ).set( "normal_y", tYNormal );
            tGENParams( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tBilinear );
            tGeoCounter++;

            // Geometry parameter lists
            tGENParams( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_threshold", 1e-16 );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "intersection_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "center_x", -0.5 );
            tGENParams( 1 )( tGeoCounter ).set( "center_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tBilinear );
            tGeoCounter++;

            // Geometry parameter lists
            tGENParams( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_threshold", 1e-16 );
            tGENParams( 1 )( tGeoCounter ).set( "isocontour_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "intersection_tolerance", 1e-12 );
            tGENParams( 1 )( tGeoCounter ).set( "center_x", 1.0 );
            tGENParams( 1 )( tGeoCounter ).set( "center_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_x", 1.0 );
            tGENParams( 1 )( tGeoCounter ).set( "normal_y", 0.0 );
            tGENParams( 1 )( tGeoCounter ).set( "use_multilinear_interpolation", tBilinear );
            tGeoCounter++;

            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // HMR parameters
            std::string tInterpolationOrder = "1";

            int                                      tRefineBuffer = 1;
            Module_Parameter_Lists tHMRParams( Module_Type::HMR );
            tHMRParams( 0 )( 0 ).set( "number_of_elements_per_dimension", 2, 1 );
            tHMRParams( 0 )( 0 ).set( "domain_dimensions", 2.0, 2.0 );
            tHMRParams( 0 )( 0 ).set( "domain_offset", -1.0, -1.0 );
            tHMRParams( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
            tHMRParams( 0 )( 0 ).set( "lagrange_orders", tInterpolationOrder );
            tHMRParams( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
            tHMRParams( 0 )( 0 ).set( "bspline_orders", tInterpolationOrder );
            tHMRParams( 0 )( 0 ).set( "bspline_pattern", std::string( "0" ) );
            tHMRParams( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
            tHMRParams( 0 )( 0 ).set( "refinement_buffer", tRefineBuffer );
            tHMRParams( 0 )( 0 ).set( "staircase_buffer", tRefineBuffer );
            // ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
            // PSEUDO Workflow
            // HMR initialize
            std::shared_ptr< hmr::HMR > pHMR = std::make_shared< hmr::HMR >( tHMRParams );

            // Geometry engine initialize
            std::shared_ptr< gen::Geometry_Engine > pGEN = std::make_shared< gen::Geometry_Engine >( tGENParams, nullptr );

            // Initialize  Background Mesh Mesh
            std::shared_ptr< mtk::Mesh_Manager > pBGMTK = std::make_shared< mtk::Mesh_Manager >();

            // Initialize output mesh from XTK
            std::shared_ptr< mtk::Mesh_Manager > pOutputMTK = std::make_shared< mtk::Mesh_Manager >();

            // XTK initialize
            std::shared_ptr< xtk::Model > pXTK = std::make_shared< xtk::Model >( tXTKParams( 0 )( 0 ) );

            // Set performer to HMR
            pHMR->set_performer( pBGMTK );

            // Set XTK Cooperations
            pXTK->set_geometry_engine( pGEN.get() );
            pXTK->set_input_performer( pBGMTK );
            pXTK->set_output_performer( pOutputMTK );

            // uniform initial refinement
            pHMR->perform_initial_refinement();

            // HMR finalize
            pHMR->perform();
            //
            pGEN->distribute_advs( pBGMTK->get_mesh_pair( 0 ) );

            // Output GEN fields, if requested
            pGEN->output_fields( pBGMTK->get_interpolation_mesh( 0 ) );

            // XTK perform - decompose - enrich - ghost - multigrid
            pXTK->perform_decomposition();
            pXTK->perform_enrichment();
        }
    }

    // Regression MB-001 (moris/doc/internal/BUG_LEDGER.md): a fixed band exactly one background cell thick, with both
    // faces on grid lines, puts every vertex of its cells on the interface. GEN returns MORIS_INDEX_MAX for those cells;
    // XTK used to substitute a MORIS_INDEX_MAX - 1 sentinel that Enriched_Integration_Mesh::setup_blockset_with_cell_clusters
    // then used as a bulk-phase index (segfault). Resolving the phase is not enough either: the band's boundary runs along
    // cell edges, where no interface facets are created, so it would be silently decoupled (MB-010). XTK must stop with
    // a clear error.
    TEST_CASE( "Band One Cell Thick on Grid Lines 2d", "[Snapping_2d]" )
    {
        if ( par_size() == 1 )
        {
            Module_Parameter_Lists tXTKParams( Module_Type::XTK );
            tXTKParams( 0 ).add_parameter_list( prm::create_xtk_parameter_list() );
            tXTKParams( 0 )( 0 ).set( "decompose", true );
            tXTKParams( 0 )( 0 ).set( "decomposition_type", "conformal" );
            tXTKParams( 0 )( 0 ).set( "enrich", true );
            tXTKParams( 0 )( 0 ).set( "basis_rank", "bspline" );
            tXTKParams( 0 )( 0 ).set( "enrich_mesh_indices", "0" );
            tXTKParams( 0 )( 0 ).set( "ghost_stab", true );
            tXTKParams( 0 )( 0 ).set( "multigrid", false );

            // band 0 <= y <= 0.5 = min( y, 0.5 - y ): exactly one cell of the 4 x 4 mesh on [-1,1]^2 (h = 0.5)
            Module_Parameter_Lists tGENParams( Module_Type::GEN );
            tGENParams( 0 ).add_parameter_list( prm::create_gen_parameter_list() );
            tGENParams( 2 ).add_parameter_list( prm::create_field_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 2 )( 0 ).set( "name", "Band bottom" );
            tGENParams( 2 )( 0 ).set( "center_x", 0.0 );
            tGENParams( 2 )( 0 ).set( "center_y", 0.0 );
            tGENParams( 2 )( 0 ).set( "normal_x", 0.0 );
            tGENParams( 2 )( 0 ).set( "normal_y", 1.0 );
            tGENParams( 2 ).add_parameter_list( prm::create_field_parameter_list( gen::Field_Type::LINE ) );
            tGENParams( 2 )( 1 ).set( "name", "Band top" );
            tGENParams( 2 )( 1 ).set( "center_x", 0.0 );
            tGENParams( 2 )( 1 ).set( "center_y", 0.5 );
            tGENParams( 2 )( 1 ).set( "normal_x", 0.0 );
            tGENParams( 2 )( 1 ).set( "normal_y", -1.0 );
            tGENParams( 1 ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::COMBINED_FIELDS ) );
            tGENParams( 1 )( 0 ).set( "dependencies", "Band bottom", "Band top" );
            tGENParams( 1 )( 0 ).set( "use_minimum", true );
            tGENParams( 1 )( 0 ).set( "isocontour_threshold", 0.0 );
            tGENParams( 1 )( 0 ).set( "isocontour_tolerance", 1e-12 );
            tGENParams( 1 )( 0 ).set( "intersection_tolerance", 1e-12 );

            Module_Parameter_Lists tHMRParams( Module_Type::HMR );
            tHMRParams( 0 )( 0 ).set( "number_of_elements_per_dimension", 4, 4 );
            tHMRParams( 0 )( 0 ).set( "domain_dimensions", 2.0, 2.0 );
            tHMRParams( 0 )( 0 ).set( "domain_offset", -1.0, -1.0 );
            tHMRParams( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
            tHMRParams( 0 )( 0 ).set( "lagrange_orders", "1" );
            tHMRParams( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
            tHMRParams( 0 )( 0 ).set( "bspline_orders", "1" );
            tHMRParams( 0 )( 0 ).set( "bspline_pattern", std::string( "0" ) );
            tHMRParams( 0 )( 0 ).set( "lagrange_to_bspline", "0" );
            tHMRParams( 0 )( 0 ).set( "refinement_buffer", 1 );
            tHMRParams( 0 )( 0 ).set( "staircase_buffer", 1 );

            std::shared_ptr< hmr::HMR >             pHMR       = std::make_shared< hmr::HMR >( tHMRParams );
            std::shared_ptr< gen::Geometry_Engine > pGEN       = std::make_shared< gen::Geometry_Engine >( tGENParams, nullptr );
            std::shared_ptr< mtk::Mesh_Manager >    pBGMTK     = std::make_shared< mtk::Mesh_Manager >();
            std::shared_ptr< mtk::Mesh_Manager >    pOutputMTK = std::make_shared< mtk::Mesh_Manager >();
            std::shared_ptr< xtk::Model >           pXTK       = std::make_shared< xtk::Model >( tXTKParams( 0 )( 0 ) );

            pHMR->set_performer( pBGMTK );
            pXTK->set_geometry_engine( pGEN.get() );
            pXTK->set_input_performer( pBGMTK );
            pXTK->set_output_performer( pOutputMTK );
            pHMR->perform_initial_refinement();
            pHMR->perform();
            pGEN->distribute_advs( pBGMTK->get_mesh_pair( 0 ) );

            // the band cells have no vertex off the interface: GEN cannot determine their phase
            mtk::Interpolation_Mesh* tIPMesh  = pBGMTK->get_mesh_pair( 0 ).get_interpolation_mesh();
            uint                     tNumBand = 0;
            for ( uint iCell = 0; iCell < tIPMesh->get_num_elems(); iCell++ )
            {
                const mtk::Cell& tCell     = tIPMesh->get_mtk_cell( iCell );
                Matrix< DDRMat > tCentroid = tCell.compute_cell_centroid();
                if ( tCentroid( 1 ) > 0.0 and tCentroid( 1 ) < 0.5 )
                {
                    tNumBand++;
                    CHECK( pGEN->get_element_phase_index( tCell ) == MORIS_INDEX_MAX );
                }
            }
            REQUIRE( tNumBand == 4 );

            // used to segfault later in Enriched_Integration_Mesh::setup_blockset_with_cell_clusters; must now stop cleanly
            CHECK_THROWS( pXTK->perform_decomposition() );
        }
    }

}    // namespace moris::xtk
