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
 */

#include "catch.hpp"

#include <cmath>

#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
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
            tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 } );

            Cut_Integration_Mesh* tCutMesh = tModel.get_cut_integration_mesh();
            REQUIRE( tCutMesh->mResolutionHiddenEdges.size() == 1 );
            uint tCount = tCutMesh->mResolutionHiddenEdges( 0 );

            delete tMesh;
            return tCount;
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
}    // namespace moris::xtk
