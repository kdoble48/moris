/*
 * Copyright (c) 2022 University of Colorado
 *Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Unenrichment.cpp
 *
 */

#include "catch.hpp"
#include "paths.hpp"
#include "cl_Matrix.hpp"

// implementations to test
#include "cl_MTK_Mesh_Factory.hpp"

#define protected public
#define private public
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Enriched_Interpolation_Mesh.hpp"
#undef protected
#undef private

#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Line.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Plane.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_MTK_Intersection_Detect.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "cl_XTK_Basis_Processor.hpp"
#include "fn_equal_to.hpp"

namespace moris::xtk
{
    TEST_CASE( "XTK Basis Extension", "[XTK],[XTK_Basis_Extension]" )
    {
        // This test is designed to test the XTK basis extention for a simple problem, the mesh is ouput to an exodus file
        // This is done by comparing the t-matrices at a given lagrange node before and after the extention
        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;

            // HMR parameter list
            moris::Module_Parameter_Lists tParameters( Module_Type::HMR );

            tParameters.set( "number_of_elements_per_dimension", 3, 3 );
            tParameters.set( "domain_dimensions", 1.0, 1.0 );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", "2" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "2" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "refinement_buffer", 1 );
            tParameters.set( "staircase_buffer", 1 );

            // create the hmr mesh
            hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            tHMR.finalize();

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            // define the sphere such that it is non interacting
            auto                                              tCircle   = std::make_shared< moris::gen::Circle >( 3.0, 0.5, 2.5 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tCircle ) };

            // define ge engine
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            // xtk parameter list
            moris::Parameter_List tXTKParameters = moris::prm::create_xtk_parameter_list();
            tXTKParameters.insert( "has_parameter_list", true );
            tXTKParameters.set( "decompose", true );
            tXTKParameters.set( "decomposition_type", "conformal" );
            tXTKParameters.set( "enrich", true );
            tXTKParameters.set( "basis_rank", "bspline" );
            tXTKParameters.set( "enrich_mesh_indices", "0" );
            tXTKParameters.set( "ghost_stab", true );
            tXTKParameters.set( "multigrid", false );
            tXTKParameters.set( "verbose", true );
            tXTKParameters.set( "print_enriched_ig_mesh", false );
            tXTKParameters.set( "exodus_output_XTK_ig_mesh", false );
            tXTKParameters.set( "high_to_low_dbl_side_sets", true );
            tXTKParameters.set( "use_SPG_based_enrichment", true );

            // dimension of the
            size_t tModelDimension = 2;

            // These are necessary set function since there is a disconnect between the functions in the xtk model
            Model tXTKModel( tModelDimension, tInterpolationMesh, &tGeometryEngine );
            tXTKModel.mVerbose            = false;
            tXTKModel.mParameterList      = tXTKParameters;
            tXTKModel.mBsplineMeshIndices = { { 0 } };

            // Specify decomposition Method and Cut Mesh --- ------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            // perform basis enrichment
            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0, false, true );

            Enriched_Interpolation_Mesh& tEnrichedIPMesh = tXTKModel.get_enriched_interp_mesh( 0 );

            Matrix< IndexMat >      tBasisIndicesBeforeExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_indices();
            const Matrix< DDRMat >* tBasisWeightsBeforeExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_weights();

            Matrix< DDRMat >   tExpectedWeights = { { 0.25, 0.25, 0.25 } };
            Matrix< IndexMat > tExpectedIndices = { { 5, 3, 13, 11 } };

            // define a generic lambda function to check if two values are equal within a threshold ( works for moris_index and real types)
            auto isEqualLambda = []( const auto& a, const auto& b ) {
                moris::real const & error_factor = 1.0;
                return moris::equal_to( a, b, error_factor );
            };

            CHECK( std::equal( tExpectedWeights.begin(), tExpectedWeights.end(), tBasisWeightsBeforeExtension->begin(), isEqualLambda ) );
            CHECK( std::equal( tExpectedIndices.begin(), tExpectedIndices.end(), tBasisIndicesBeforeExtension.begin(), isEqualLambda ) );

            // create basis processor object
            Basis_Processor tBasisProcessor( &tXTKModel );
            tBasisProcessor.perform_basis_extention();

            // get the t-matrix info after extention
            Matrix< IndexMat >      tBasisIndicesAfterExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_indices();
            const Matrix< DDRMat >* tBasisWeightsAfterExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_weights();

            // update the expected values after extention ( they are transposed )
            tExpectedWeights = { { 0.25, 1, -0.75, 1, 0.25, -0.75 } };
            tExpectedIndices = { { 0, 3, 1, 11, 8, 9 } };

            CHECK( std::equal( tExpectedWeights.begin(), tExpectedWeights.end(), tBasisWeightsAfterExtension->begin(), isEqualLambda ) );
            CHECK( std::equal( tExpectedIndices.begin(), tExpectedIndices.end(), tBasisIndicesAfterExtension.begin(), isEqualLambda ) );

        }    // end if: serial test run

    }        // end TEST_CASE: "XTK Basis Extension"

    // ----------------------------------------------------------------------------------

    TEST_CASE( "XTK Cell Agglomeration", "[XTK],[XTK_Cell_Agglomeration]" )
    {
        // This test is designed to test the XTK cell agglomeration for a simple problem, the mesh is ouput to an exodus file
        // This is done by comparing the t-matrices at a given lagrange node before and after the agglomeration
        if ( par_size() == 1 )
        {
            uint tLagrangeMeshIndex = 0;

            // HMR parameter list
            moris::Module_Parameter_Lists tParameters( Module_Type::HMR );

            tParameters.set( "number_of_elements_per_dimension", 3, 3 );
            tParameters.set( "domain_dimensions", 1.0, 1.0 );
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", "1" );
            tParameters.set( "lagrange_pattern", "0" );
            tParameters.set( "bspline_orders", "1" );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "refinement_buffer", 1 );
            tParameters.set( "staircase_buffer", 1 );

            // create the hmr mesh
            hmr::HMR tHMR( tParameters );

            tHMR.perform_initial_refinement();

            tHMR.finalize();

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( tLagrangeMeshIndex );

            // define the sphere such that it is non interacting
            auto                                              tCircle   = std::make_shared< moris::gen::Circle >( 3.0, 0.5, 2.5 );
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tCircle ) };

            // define ge engine
            moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = tGeometry;
            moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh, tGeometryEngineParameters );

            // xtk parameter list
            moris::Parameter_List tXTKParameters = moris::prm::create_xtk_parameter_list();
            tXTKParameters.insert( "has_parameter_list", true );
            tXTKParameters.set( "decompose", true );
            tXTKParameters.set( "decomposition_type", "conformal" );
            tXTKParameters.set( "enrich", true );
            tXTKParameters.set( "basis_rank", "bspline" );
            tXTKParameters.set( "enrich_mesh_indices", "0" );
            tXTKParameters.set( "ghost_stab", true );
            tXTKParameters.set( "multigrid", false );
            tXTKParameters.set( "verbose", true );
            tXTKParameters.set( "print_enriched_ig_mesh", false );
            tXTKParameters.set( "exodus_output_XTK_ig_mesh", false );
            tXTKParameters.set( "high_to_low_dbl_side_sets", true );
            tXTKParameters.set( "use_SPG_based_enrichment", true );

            // dimension of the
            size_t tModelDimension = 2;

            // These are necessary set function since there is a disconnect between the functions in the xtk model
            Model tXTKModel( tModelDimension, tInterpolationMesh, &tGeometryEngine );
            tXTKModel.mVerbose            = false;
            tXTKModel.mParameterList      = tXTKParameters;
            tXTKModel.mBsplineMeshIndices = { { 0 } };

            // Specify decomposition Method and Cut Mesh --- ------------------------------------
            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            tXTKModel.decompose( tDecompositionMethods );

            // perform basis enrichment based on the SPGs
            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0, false, true );

            // get the enriched interpolation mesh
            Enriched_Interpolation_Mesh& tEnrichedIPMesh = tXTKModel.get_enriched_interp_mesh( 0 );

            // get t-matrices before extention
            Matrix< IndexMat >      tBasisIndicesBeforeExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_indices();
            const Matrix< DDRMat >* tBasisWeightsBeforeExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_weights();

            // initialize the expected values fore the t-matrices
            Matrix< DDRMat >   tExpectedWeights = { { 1.0 } };
            Matrix< IndexMat > tExpectedIndices = { { 8 } };

            // define a generic lambda function to check if two values are equal within a threshold ( works for moris_index and real types)
            auto isEqualLambda = []( const auto& aA, const auto& aB ) {
                moris::real const & tErrorFactor = 1.0;
                return moris::equal_to( aA, aB, tErrorFactor );
            };

            // check if the t-matrices are correct before extention
            CHECK( std::equal( tExpectedWeights.begin(), tExpectedWeights.end(), tBasisWeightsBeforeExtension->begin(), isEqualLambda ) );
            CHECK( std::equal( tExpectedIndices.begin(), tExpectedIndices.end(), tBasisIndicesBeforeExtension.begin(), isEqualLambda ) );

            // create basis processor object and perform cell agglomeration
            Basis_Processor tBasisProcessor( &tXTKModel );
            tBasisProcessor.perform_cell_agglomeration();

            // get the t-matrix info after extention
            Matrix< IndexMat >      tBasisIndicesAfterExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_indices();
            const Matrix< DDRMat >* tBasisWeightsAfterExtension = tEnrichedIPMesh.get_mtk_vertex( 9 ).get_interpolation( 0 )->get_weights();

            // update the expected values after extention ( they are transposed )
            tExpectedWeights = { { 2.0, -1.0 } };
            tExpectedIndices = { { 10, 11 } };

            // check if the t-matrices are correct after extention
            CHECK( std::equal( tExpectedWeights.begin(), tExpectedWeights.end(), tBasisWeightsAfterExtension->begin(), isEqualLambda ) );
            CHECK( std::equal( tExpectedIndices.begin(), tExpectedIndices.end(), tBasisIndicesAfterExtension.begin(), isEqualLambda ) );

        }    // end if: serial test run

    }        // end TEST_CASE: "XTK Cell Agglomeration"

    // ----------------------------------------------------------------------------------

    /**
     * Builds a 2D mesh of aNumElemsPerDim x aNumElemsPerDim linear elements on the unit square, cuts it with vertical
     * lines at the given x positions, enriches it based on subphase groups and runs the cell agglomeration.
     */
    static void
    cut_with_vertical_lines_and_agglomerate(
            int                    aNumElemsPerDim,
            Vector< real > const & aLineXPositions )
    {
        // HMR parameter list
        moris::Module_Parameter_Lists tParameters( Module_Type::HMR );

        tParameters.set( "number_of_elements_per_dimension", aNumElemsPerDim, aNumElemsPerDim );
        tParameters.set( "domain_dimensions", 1.0, 1.0 );
        tParameters.set( "lagrange_output_meshes", "0" );

        tParameters.set( "lagrange_orders", "1" );
        tParameters.set( "lagrange_pattern", "0" );
        tParameters.set( "bspline_orders", "1" );
        tParameters.set( "bspline_pattern", "0" );

        tParameters.set( "lagrange_to_bspline", "0" );

        tParameters.set( "refinement_buffer", 1 );
        tParameters.set( "staircase_buffer", 1 );

        // create the hmr mesh
        hmr::HMR tHMR( tParameters );
        tHMR.perform_initial_refinement();
        tHMR.finalize();

        std::shared_ptr< moris::hmr::Interpolation_Mesh_HMR > tInterpolationMesh( tHMR.create_interpolation_mesh( 0 ) );

        // vertical lines (none of them coincides with an element edge)
        Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry;
        for ( real tLineXPosition : aLineXPositions )
        {
            tGeometry.push_back( std::make_shared< gen::Level_Set_Geometry >( std::make_shared< moris::gen::Line >( tLineXPosition, 0.5, 1.0, 0.0 ) ) );
        }

        // define the geometry engine
        moris::gen::Geometry_Engine_Parameters tGeometryEngineParameters;
        tGeometryEngineParameters.mGeometries = tGeometry;
        moris::gen::Geometry_Engine tGeometryEngine( tInterpolationMesh.get(), tGeometryEngineParameters );

        // xtk parameter list
        moris::Parameter_List tXTKParameters = moris::prm::create_xtk_parameter_list();
        tXTKParameters.insert( "has_parameter_list", true );
        tXTKParameters.set( "decompose", true );
        tXTKParameters.set( "decomposition_type", "conformal" );
        tXTKParameters.set( "enrich", true );
        tXTKParameters.set( "basis_rank", "bspline" );
        tXTKParameters.set( "enrich_mesh_indices", "0" );
        tXTKParameters.set( "ghost_stab", false );
        tXTKParameters.set( "multigrid", false );
        tXTKParameters.set( "verbose", false );
        tXTKParameters.set( "print_enriched_ig_mesh", false );
        tXTKParameters.set( "exodus_output_XTK_ig_mesh", false );
        tXTKParameters.set( "high_to_low_dbl_side_sets", true );
        tXTKParameters.set( "use_SPG_based_enrichment", true );

        // decompose and enrich
        Model tXTKModel( 2, tInterpolationMesh.get(), &tGeometryEngine );
        tXTKModel.mVerbose            = false;
        tXTKModel.mParameterList      = tXTKParameters;
        tXTKModel.mBsplineMeshIndices = { { 0 } };

        Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
        tXTKModel.decompose( tDecompositionMethods );
        tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0, false, true );

        // agglomerate the cut cells
        Basis_Processor tBasisProcessor( &tXTKModel );
        tBasisProcessor.perform_cell_agglomeration();
    }

    // ----------------------------------------------------------------------------------

    TEST_CASE( "XTK Cell Agglomeration Without Root Cells", "[XTK],[XTK_Cell_Agglomeration]" )
    {
        // Cell agglomeration takes the basis of every cut subphase group (SPG) from an uncut SPG of the same phase that is
        // reachable through the SPG neighbourhood (its root). An SPG without a root must stop the run with an error instead
        // of being left with root -1 (MB-110: the search gave up silently and override_t_matrices then read SPG index -1).
        if ( par_size() == 1 )
        {
            SECTION( "every B-spline cell is cut: no SPG has a root" )
            {
                // 2 x 2 elements: the line at x = 0.45 cuts the left column, the line at x = 0.55 the right one
                REQUIRE_THROWS_WITH( cut_with_vertical_lines_and_agglomerate( 2, { 0.45, 0.55 } ),
                        Catch::Matchers::Contains( "8 subphase group(s) found no root B-spline cell" ) );
            }

            SECTION( "a strip thinner than a cell: its SPGs have no root" )
            {
                // 3 x 3 elements: the outer phases reach uncut cells, the strip 0.45 < x < 0.55 lies inside the cut middle column
                REQUIRE_THROWS_WITH( cut_with_vertical_lines_and_agglomerate( 3, { 0.45, 0.55 } ),
                        Catch::Matchers::Contains( "3 subphase group(s) found no root B-spline cell" ) );
            }

            SECTION( "every SPG reaches an uncut cell" )
            {
                // 3 x 3 elements cut by a single line through the middle column
                REQUIRE_NOTHROW( cut_with_vertical_lines_and_agglomerate( 3, { 0.45 } ) );
            }
        }
    }    // end TEST_CASE: "XTK Cell Agglomeration Without Root Cells"

    // ----------------------------------------------------------------------------------

}    // namespace moris::xtk
