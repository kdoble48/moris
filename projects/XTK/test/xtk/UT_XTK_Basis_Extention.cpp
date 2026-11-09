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
#include "cl_XTK_Basis_Processor.hpp"
#undef protected
#undef private

#include "cl_XTK_Enriched_Integration_Mesh.hpp"

#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR.hpp"
#include "cl_GEN_Circle.hpp"
#include "cl_GEN_Sphere.hpp"
#include "cl_GEN_Plane.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "fn_PRM_XTK_Parameters.hpp"
#include "cl_MTK_Intersection_Detect.hpp"
#include "cl_MTK_Writer_Exodus.hpp"
#include "fn_equal_to.hpp"
#include "cl_HMR_Lagrange_Mesh_Base.hpp"
#include "cl_HMR_BSpline_Mesh_Base.hpp"
#include "cl_HMR_Basis.hpp"
#include "cl_HMR_Background_Mesh_Base.hpp"

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

    namespace
    {
        /**
         * @brief linear test field u = 0.3 + 1.7 x - 0.9 y (+ 0.6 z in 3D) used to check that the T-matrices reproduce linear functions
         */
        real
        linear_test_field( Matrix< DDRMat > const & aCoords )
        {
            real tValue = 0.3 + 1.7 * aCoords( 0 ) - 0.9 * aCoords( 1 );
            if ( aCoords.numel() > 2 )
            {
                tValue += 0.6 * aCoords( 2 );
            }
            return tValue;
        }

        // ----------------------------------------------------------------------------------

        /**
         * @brief Greville point of a B-spline basis on a uniform level-0 B-spline mesh in a serial run (the formula of
         * hmr::BSpline_Mesh::calculate_basis_coordinates(), whose result get_coords() only returns in debug builds)
         */
        Matrix< DDRMat >
        greville_point_of_basis(
                hmr::BSpline_Mesh_Base* aBSplineMesh,
                moris_index             aBasisIndex,
                real                    aElementLength )
        {
            hmr::Basis const * tBasis       = aBSplineMesh->get_basis_by_index( aBasisIndex );
            luint const *      tIJK         = tBasis->get_ijk();
            real               tOrder       = (real)aBSplineMesh->Mesh_Base::get_order();
            Matrix< DDRMat >   tDomainShift = aBSplineMesh->get_background_mesh()->get_domain_offset();

            uint tSpatialDim = tDomainShift.numel();

            REQUIRE( tBasis->get_level() == 0 );

            Matrix< DDRMat > tCoords( 1, tSpatialDim );
            for ( uint iDim = 0; iDim < tSpatialDim; iDim++ )
            {
                tCoords( iDim ) = ( 0.5 * ( tOrder + 1.0 ) - tOrder + (real)tIJK[ iDim ] ) * aElementLength + tDomainShift( iDim );
            }
            return tCoords;
        }

        // ----------------------------------------------------------------------------------

        /**
         * @brief evaluates the linear test field through the T-matrix of an enriched interpolation vertex, using the value of
         * the field at the Greville point of each B-spline basis as its coefficient (exact for uniform B-splines)
         */
        real
        interpolate_linear_test_field(
                Matrix< IndexMat > const &    aBasisIndices,
                Matrix< DDRMat > const &      aBasisWeights,
                Vector< moris_index > const & aNonEnrichedBasisIndices,
                hmr::BSpline_Mesh_Base*       aBSplineMesh,
                real                          aElementLength )
        {
            real tValue = 0.0;
            for ( uint iBasis = 0; iBasis < aBasisIndices.numel(); iBasis++ )
            {
                moris_index      tNonEnrichedIndex = aNonEnrichedBasisIndices( aBasisIndices( iBasis ) );
                Matrix< DDRMat > tBasisCoords      = greville_point_of_basis( aBSplineMesh, tNonEnrichedIndex, aElementLength );
                tValue += aBasisWeights( iBasis ) * linear_test_field( tBasisCoords );
            }
            return tValue;
        }

        // ----------------------------------------------------------------------------------

        /**
         * @brief MB-108: runs cell agglomeration on a 6 x 6 (x 6) mesh cut by a circle (sphere) and checks every vertex of
         * the enriched interpolation mesh afterwards: partition of unity, reproduction of a linear field, no repeated basis in
         * a row and no follower basis left in any row.
         *
         * @param aSpatialDim spatial dimension, 2 or 3
         * @param aLagrangeOrder order of the Lagrange (interpolation) mesh
         * @param aBsplineOrder order of the B-spline mesh
         * @param aRefineLagrange put the Lagrange mesh on a pattern refined once beyond the B-spline pattern
         */
        void
        check_cell_agglomeration_t_matrices(
                uint aSpatialDim,
                uint aLagrangeOrder,
                uint aBsplineOrder,
                bool aRefineLagrange )
        {
            // HMR parameter list
            moris::Module_Parameter_Lists tParameters( Module_Type::HMR );

            if ( aSpatialDim == 3 )
            {
                tParameters.set( "number_of_elements_per_dimension", 6, 6, 6 );
                tParameters.set( "domain_dimensions", 1.0, 1.0, 1.0 );
            }
            else
            {
                tParameters.set( "number_of_elements_per_dimension", 6, 6 );
                tParameters.set( "domain_dimensions", 1.0, 1.0 );
            }
            tParameters.set( "lagrange_output_meshes", "0" );

            tParameters.set( "lagrange_orders", std::to_string( aLagrangeOrder ) );
            tParameters.set( "lagrange_pattern", aRefineLagrange ? "1" : "0" );
            tParameters.set( "bspline_orders", std::to_string( aBsplineOrder ) );
            tParameters.set( "bspline_pattern", "0" );

            tParameters.set( "lagrange_to_bspline", "0" );

            tParameters.set( "refinement_buffer", 1 );
            tParameters.set( "staircase_buffer", 1 );

            if ( aRefineLagrange )
            {
                tParameters.set( "pattern_initial_refinement", 0, 1 );
            }

            // create the hmr mesh
            hmr::HMR tHMR( tParameters );
            tHMR.perform_initial_refinement();
            tHMR.finalize();

            moris::hmr::Interpolation_Mesh_HMR* tInterpolationMesh = tHMR.create_interpolation_mesh( 0 );

            // circle (sphere) that cuts several cells, the inside keeps uncut cells that can serve as root cells
            std::shared_ptr< moris::gen::Field > tLevelSet;
            if ( aSpatialDim == 3 )
            {
                tLevelSet = std::make_shared< moris::gen::Sphere >( 0.51, 0.48, 0.47, 0.34 );
            }
            else
            {
                tLevelSet = std::make_shared< moris::gen::Circle >( 0.51, 0.48, 0.34 );
            }
            Vector< std::shared_ptr< moris::gen::Geometry > > tGeometry = { std::make_shared< gen::Level_Set_Geometry >( tLevelSet ) };

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
            tXTKParameters.set( "ghost_stab", false );
            tXTKParameters.set( "multigrid", false );
            tXTKParameters.set( "verbose", false );
            tXTKParameters.set( "print_enriched_ig_mesh", false );
            tXTKParameters.set( "exodus_output_XTK_ig_mesh", false );
            tXTKParameters.set( "high_to_low_dbl_side_sets", true );
            tXTKParameters.set( "use_SPG_based_enrichment", true );

            Model tXTKModel( aSpatialDim, tInterpolationMesh, &tGeometryEngine );
            tXTKModel.mVerbose            = false;
            tXTKModel.mParameterList      = tXTKParameters;
            tXTKModel.mBsplineMeshIndices = { { 0 } };

            Vector< enum Subdivision_Method > tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 };
            if ( aSpatialDim == 3 )
            {
                tDecompositionMethods = { Subdivision_Method::NC_REGULAR_SUBDIVISION_HEX8, Subdivision_Method::C_HIERARCHY_TET4 };
            }
            tXTKModel.decompose( tDecompositionMethods );

            tXTKModel.perform_basis_enrichment( mtk::EntityRank::BSPLINE, 0, false, true );

            Enriched_Interpolation_Mesh&  tEnrichedIPMesh          = tXTKModel.get_enriched_interp_mesh( 0 );
            Vector< moris_index > const & tNonEnrichedBasisIndices = tXTKModel.mEnrichment->mEnrichmentData( 0 ).mNonEnrBfIndForEnrBfInd;
            hmr::BSpline_Mesh_Base*       tBSplineMesh             = tInterpolationMesh->get_lagrange_mesh()->get_bspline_mesh( 0 );
            uint                          tNumVertices             = tEnrichedIPMesh.get_num_entities( mtk::EntityRank::NODE );
            real                          tElementLength           = 1.0 / 6.0;

            // save the T-matrix rows before the agglomeration and check that the reference field is reproduced by them
            Vector< Matrix< IndexMat > > tIndicesBefore( tNumVertices );
            Vector< Matrix< DDRMat > >   tWeightsBefore( tNumVertices );
            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                mtk::Vertex const &               tVertex        = tEnrichedIPMesh.get_mtk_vertex( iVertex );
                mtk::Vertex_Interpolation const * tInterpolation = tVertex.get_interpolation( 0 );
                tIndicesBefore( iVertex )                        = tInterpolation->get_indices();
                tWeightsBefore( iVertex )                        = *tInterpolation->get_weights();

                real tExpected = linear_test_field( tVertex.get_coords() );
                real tValue    = interpolate_linear_test_field( tIndicesBefore( iVertex ), tWeightsBefore( iVertex ), tNonEnrichedBasisIndices, tBSplineMesh, tElementLength );
                REQUIRE( std::abs( tValue - tExpected ) < 1.0e-12 );
            }

            // perform cell agglomeration (MB-108: used to abort with "Mat::operator(): index out of bounds")
            Basis_Processor tBasisProcessor( &tXTKModel );
            REQUIRE_NOTHROW( tBasisProcessor.perform_cell_agglomeration() );

            Vector< moris::moris_index > const & tFollowerBasis       = tBasisProcessor.mBasisData( 0 ).mFollowerBasis;
            Vector< moris::moris_index > const & tFollowerOwningCells = tBasisProcessor.mBasisData( 0 ).mFollowerBasisOwningCell;

            uint tMaxFollowersPerRow = 0;
            uint tNumChangedRows     = 0;
            real tMaxUnityError      = 0.0;
            real tMaxLinearError     = 0.0;
            uint tNumRepeatedBases   = 0;
            uint tNumFollowersLeft   = 0;

            for ( uint iVertex = 0; iVertex < tNumVertices; iVertex++ )
            {
                // count the follower bases with a root cell in the row before the agglomeration
                uint tNumFollowers = 0;
                for ( uint iBasis = 0; iBasis < tIndicesBefore( iVertex ).numel(); iBasis++ )
                {
                    moris_index tBasisIndex = tIndicesBefore( iVertex )( iBasis );
                    if ( tFollowerBasis( tBasisIndex ) == 1 and tFollowerOwningCells( tBasisIndex ) != gNoIndex )
                    {
                        tNumFollowers++;
                    }
                }
                tMaxFollowersPerRow = std::max( tMaxFollowersPerRow, tNumFollowers );

                mtk::Vertex const &               tVertex        = tEnrichedIPMesh.get_mtk_vertex( iVertex );
                mtk::Vertex_Interpolation const * tInterpolation = tVertex.get_interpolation( 0 );
                Matrix< IndexMat >                tIndices       = tInterpolation->get_indices();
                Matrix< DDRMat >                  tWeights       = *tInterpolation->get_weights();

                if ( tNumFollowers > 0 )
                {
                    tNumChangedRows++;
                }

                // partition of unity
                real tSum = 0.0;
                for ( uint iBasis = 0; iBasis < tWeights.numel(); iBasis++ )
                {
                    tSum += tWeights( iBasis );
                }
                tMaxUnityError = std::max( tMaxUnityError, std::abs( tSum - 1.0 ) );

                // reproduction of a linear field
                real tExpected  = linear_test_field( tVertex.get_coords() );
                real tValue     = interpolate_linear_test_field( tIndices, tWeights, tNonEnrichedBasisIndices, tBSplineMesh, tElementLength );
                tMaxLinearError = std::max( tMaxLinearError, std::abs( tValue - tExpected ) );

                // each basis appears once and no follower basis with a root cell is left
                for ( uint iBasis = 0; iBasis < tIndices.numel(); iBasis++ )
                {
                    for ( uint jBasis = iBasis + 1; jBasis < tIndices.numel(); jBasis++ )
                    {
                        if ( tIndices( iBasis ) == tIndices( jBasis ) )
                        {
                            tNumRepeatedBases++;
                        }
                    }

                    if ( tFollowerBasis( tIndices( iBasis ) ) == 1 and tFollowerOwningCells( tIndices( iBasis ) ) != gNoIndex )
                    {
                        tNumFollowersLeft++;
                    }
                }
            }

            // higher-order and refined Lagrange meshes have nodes with two or more follower bases of a root cell, which used to
            // overflow the row (MB-108)
            if ( aLagrangeOrder > 1 or aRefineLagrange )
            {
                CHECK( tMaxFollowersPerRow >= 2 );
            }
            CHECK( tNumChangedRows > 0 );
            CHECK( tMaxUnityError < 1.0e-12 );
            CHECK( tMaxLinearError < 1.0e-12 );
            CHECK( tNumRepeatedBases == 0 );
            CHECK( tNumFollowersLeft == 0 );
        }
    }    // namespace

    // ----------------------------------------------------------------------------------

    TEST_CASE( "XTK Cell Agglomeration Higher Order Lagrange", "[XTK],[XTK_Cell_Agglomeration]" )
    {
        // MB-108: cell agglomeration on Lagrange nodes that interpolate from two or more follower bases of a root cell
        if ( par_size() == 1 )
        {
            SECTION( "Lagrange order 2, linear B-splines" )
            {
                check_cell_agglomeration_t_matrices( 2, 2, 1, false );
            }

            SECTION( "Lagrange order 1 refined once beyond linear B-splines" )
            {
                check_cell_agglomeration_t_matrices( 2, 1, 1, true );
            }

            SECTION( "Lagrange order 2, quadratic B-splines" )
            {
                check_cell_agglomeration_t_matrices( 2, 2, 2, false );
            }

            SECTION( "3D, Lagrange order 1, linear B-splines" )
            {
                check_cell_agglomeration_t_matrices( 3, 1, 1, false );
            }

            SECTION( "3D, Lagrange order 2, linear B-splines" )
            {
                check_cell_agglomeration_t_matrices( 3, 2, 1, false );
            }
        }
    }

    // ----------------------------------------------------------------------------------

}    // namespace moris::xtk
