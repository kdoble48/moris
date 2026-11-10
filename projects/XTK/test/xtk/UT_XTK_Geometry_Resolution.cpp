/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Geometry_Resolution.cpp
 *
 * MB-084: the geometry resolution check of the cut (XTK geometry_resolution_samples) counts background edges that an
 * analytic interface crosses more often than the cut, against the edge-sign maps of the donut / U-channel benchmark
 * (runs/studies/moris_immersed_pathologies/donut_uchannel, Figures 1C and 2C: 6 edges crossed twice for the 0.3 h slot
 * in an element row, none for the slot on a node line, 2 for the ring at n = 4).
 *
 * MB-095: the check gives the same counts and measures on every number of ranks (each background cell counted by its
 * owner, each hidden edge once). The meshes put a processor boundary through the hidden edges (2 ranks split x,
 * 4 ranks split x and y).
 *
 */

#include "catch.hpp"

#include <cmath>
#include <iomanip>
#include <iostream>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_MTK_Mesh_Manager.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_User_Defined_Field.hpp"

namespace moris::xtk
{
    namespace
    {
        real gYs = 2.5123;

        real
        rect( real aX, real aY, real aX0, real aX1, real aY0, real aY1 )
        {
            return std::min( std::min( aX - aX0, aX1 - aX ), std::min( aY - aY0, aY1 - aY ) );
        }

        // U-channel, slot 0.3 wide centred at y = gYs (Func_Geom of du_deck.cpp, case 1, legs 1.2)
        real
        u_channel( const Matrix< DDRMat >& aPoint, const Vector< real >& )
        {
            real tOuter = rect( aPoint( 0 ), aPoint( 1 ), 0.9371, 8.0371, gYs - 1.35, gYs + 1.35 );
            real tSlot  = rect( aPoint( 0 ), aPoint( 1 ), 2.0371, 9.0371, gYs - 0.15, gYs + 0.15 );
            return std::max( std::min( tOuter, -tSlot ), 1.0371 - aPoint( 0 ) );
        }

        // donut (case 0): ring R = 1, t = 0.2 at (3.75, 1.25), bars 0.2 wide, pads x < 1.2913, x > 6.2087
        real
        donut( const Matrix< DDRMat >& aPoint, const Vector< real >& )
        {
            real tX    = aPoint( 0 );
            real tY    = aPoint( 1 );
            real tR    = std::hypot( tX - 3.75, tY - 1.25 );
            real tRing = std::min( 1.0 - tR, tR - 0.8 );
            real tBarL = rect( tX, tY, -1.0, 3.75 - 1.0 + 0.1, 1.15, 1.35 );
            real tBarR = rect( tX, tY, 3.75 + 1.0 - 0.1, 8.5, 1.15, 1.35 );
            return std::max( std::max( tRing, std::max( tBarL, tBarR ) ), std::max( 1.2913 - tX, tX - 6.2087 ) );
        }

        uint
        hidden_edges( gen::Field_Function aFunction, const Vector< uint >& aNumElements, const Vector< real >& aDimensions )
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
            hmr::HMR tHMR( tParameters );
            tHMR.finalize();
            hmr::Interpolation_Mesh_HMR* tMesh = tHMR.create_interpolation_mesh( 0 );

            gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = { std::make_shared< gen::Level_Set_Geometry >( std::make_shared< gen::User_Defined_Field >( aFunction ) ) };
            gen::Geometry_Engine tEngine( tMesh, tGeometryEngineParameters );

            Model tModel( 2, tMesh, &tEngine );
            tModel.mVerbose = false;
            tModel.get_parameter_list().insert( "geometry_resolution_samples", uint( 2 ) );    // the check is off by default
            tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 } );

            Cut_Integration_Mesh* tCutMesh = tModel.get_cut_integration_mesh();
            REQUIRE( tCutMesh->mResolutionHiddenEdges.size() == 1 );
            uint tCount = tCutMesh->mResolutionHiddenEdges( 0 );

            delete tMesh;
            return tCount;
        }

        // MB-095: the resolution check of one geometry on any number of ranks (HMR splits the mesh over the ranks)
        struct Resolution_Counts
        {
            uint mHiddenEdges     = 0;
            uint mUnseenCells     = 0;
            uint mTopologyCells   = 0;
            real mSampledMeasure  = 0.0;
            real mTemplateMeasure = 0.0;
        };

        Resolution_Counts
        resolution_counts( gen::Field_Function aFunction, uint aNumElementsX, uint aNumElementsY, real aLengthX, real aLengthY )
        {
            Module_Parameter_Lists tParameters( Module_Type::HMR );
            tParameters( 0 ).add_parameter_list( prm::create_hmr_parameter_list() );
            tParameters( 0 )( 0 ).set( "number_of_elements_per_dimension", aNumElementsX, aNumElementsY );
            tParameters( 0 )( 0 ).set( "domain_dimensions", aLengthX, aLengthY );
            tParameters( 0 )( 0 ).set( "domain_offset", 0.0, 0.0 );
            tParameters( 0 )( 0 ).set( "lagrange_output_meshes", "0" );
            tParameters( 0 )( 0 ).set( "lagrange_orders", "1" );
            tParameters( 0 )( 0 ).set( "lagrange_pattern", std::string( "0" ) );
            tParameters( 0 )( 0 ).set( "bspline_orders", "1" );
            tParameters( 0 )( 0 ).set( "bspline_pattern", "0" );
            tParameters( 0 )( 0 ).set( "lagrange_to_bspline", "0" );

            std::shared_ptr< hmr::HMR >          tHMR         = std::make_shared< hmr::HMR >( tParameters );
            std::shared_ptr< mtk::Mesh_Manager > tMeshManager = std::make_shared< mtk::Mesh_Manager >();
            tHMR->set_performer( tMeshManager );
            tHMR->perform_initial_refinement();
            tHMR->perform();
            mtk::Interpolation_Mesh* tMesh = tMeshManager->get_interpolation_mesh( 0 );

            gen::Geometry_Engine_Parameters tGeometryEngineParameters;
            tGeometryEngineParameters.mGeometries = { std::make_shared< gen::Level_Set_Geometry >( std::make_shared< gen::User_Defined_Field >( aFunction ) ) };
            gen::Geometry_Engine tEngine( tMesh, tGeometryEngineParameters );

            Model tModel( 2, tMesh, &tEngine );
            tModel.mVerbose = false;
            tModel.get_parameter_list().insert( "geometry_resolution_samples", uint( 2 ) );    // the check is off by default
            tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 } );

            Cut_Integration_Mesh* tCutMesh = tModel.get_cut_integration_mesh();
            REQUIRE( tCutMesh->mResolutionHiddenEdges.size() == 1 );
            Resolution_Counts tCounts;
            tCounts.mHiddenEdges     = tCutMesh->mResolutionHiddenEdges( 0 );
            tCounts.mUnseenCells     = tCutMesh->mResolutionUnseenCells( 0 );
            tCounts.mTopologyCells   = tCutMesh->mResolutionTopologyCells( 0 );
            tCounts.mSampledMeasure  = tCutMesh->mResolutionSampledMeasure( 0 );
            tCounts.mTemplateMeasure = tCutMesh->mResolutionTemplateMeasure( 0 );
            return tCounts;
        }
    }    // namespace

    TEST_CASE( "MB-084 resolution check counts edges crossed twice", "[XTK][MB-084][MB-084_CHECK]" )
    {
        if ( par_size() == 1 )
        {
            // U-channel, coarse h = 1: slot centred in an element row -> the 6 vertical edges x = 2 ... 7 through it
            gYs = 2.5123;
            CHECK( hidden_edges( &u_channel, { 9, 5 }, { 9.0, 5.0 } ) == 6 );

            // the same slot on a node line: every edge changes sign
            gYs = 2.0123;
            CHECK( hidden_edges( &u_channel, { 9, 5 }, { 9.0, 5.0 } ) == 0 );
            gYs = 2.5123;

            // donut at n = 4 (h = 0.625): the 2 vertical edges x = 3.75 through the ring top and bottom
            CHECK( hidden_edges( &donut, { 12, 4 }, { 7.5, 2.5 } ) == 2 );
        }
    }

    TEST_CASE( "MB-095 resolution check counts do not depend on the number of ranks", "[XTK][MB-095][MB-084_CHECK]" )
    {
        // the expected values are the serial ones (np 1)
        const uint tHiddenChannel   = 12;    // the vertical edges x = 2.5 ... 8 through the slot
        const uint tTopologyChannel = 12;
        const real tSampledChannel  = 23.375;
        const real tTemplateChannel = 24.125;
        const uint tHiddenDonut     = 2;
        const uint tTopologyDonut   = 4;
        const real tSampledDonut    = 9.375;
        const real tTemplateDonut   = 8.984375;

        // U-channel on h = 0.5, slot inside the element row y = 2.5 ... 3: 2 ranks split it at x = 4.5, 4 ranks also at y = 2.5
        gYs                        = 2.7623;
        Resolution_Counts tChannel = resolution_counts( &u_channel, 18, 10, 9.0, 5.0 );
        gYs                        = 2.5123;

        // donut at n = 4 (h = 0.625): the 2 hidden edges x = 3.75 lie on the processor boundary of 2 and 4 ranks
        Resolution_Counts tDonut = resolution_counts( &donut, 12, 4, 7.5, 2.5 );

        if ( par_rank() == 0 )
        {
            std::cout << std::setprecision( 17 ) << "MB-095 np " << par_size() << ": U-channel hidden " << tChannel.mHiddenEdges << " unseen "
                      << tChannel.mUnseenCells << " topology " << tChannel.mTopologyCells << " sampled " << tChannel.mSampledMeasure << " template "
                      << tChannel.mTemplateMeasure << " | donut hidden " << tDonut.mHiddenEdges << " unseen " << tDonut.mUnseenCells << " topology "
                      << tDonut.mTopologyCells << " sampled " << tDonut.mSampledMeasure << " template " << tDonut.mTemplateMeasure << "\n"
                      << std::setprecision( 6 ) << std::flush;
        }

        CHECK( tChannel.mHiddenEdges == tHiddenChannel );
        CHECK( tChannel.mUnseenCells == 0u );
        CHECK( tChannel.mTopologyCells == tTopologyChannel );
        CHECK( tChannel.mSampledMeasure == Approx( tSampledChannel ).epsilon( 1e-12 ) );
        CHECK( tChannel.mTemplateMeasure == Approx( tTemplateChannel ).epsilon( 1e-12 ) );
        CHECK( tDonut.mHiddenEdges == tHiddenDonut );
        CHECK( tDonut.mUnseenCells == 0u );
        CHECK( tDonut.mTopologyCells == tTopologyDonut );
        CHECK( tDonut.mSampledMeasure == Approx( tSampledDonut ).epsilon( 1e-12 ) );
        CHECK( tDonut.mTemplateMeasure == Approx( tTemplateDonut ).epsilon( 1e-12 ) );
    }
}    // namespace moris::xtk
