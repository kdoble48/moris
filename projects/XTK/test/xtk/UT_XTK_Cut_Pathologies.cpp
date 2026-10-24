/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_XTK_Cut_Pathologies.cpp
 *
 * Characterization sweep of the conformal decomposition (regular subdivision + node hierarchy) on a single HEX8
 * background cell [0,1]^3. The level set is the trilinear interpolant of 8 corner values; every combination of corner
 * values from a small value set is decomposed and checked against invariants that must hold for any input:
 *
 *   - the decomposition does not throw,
 *   - the integration cells tile the cell (total volume 1),
 *   - every integration cell has a bulk phase,
 *   - no integration cell has both a POSITIVE and a NEGATIVE vertex for any geometry,
 *   - every facet has one or two cells, and facets with one cell lie on the background cell boundary,
 *   - (interpolated and analytic fields) the negative volume equals the volume of the piecewise-linear interpolant on
 *     the 24-tet subdivision (face and cell centre values are corner averages), computed independently below.
 *
 * Three GEN field regimes are swept: STORED (discretization_mesh_index = -1: derived values interpolated from their
 * parents, the regime of discretized level sets), ANALYTIC (derived values evaluated exactly) and MULTILINEAR
 * (use_multilinear_interpolation). With two geometries the STORED regime must not depend on the geometry order: cutting
 * by one geometry only subdivides the piecewise-linear representation of the other.
 *
 * Three plane sweeps report how cuts accumulate (1 to 6 planes, both orders), two parallel planes approaching, and two
 * planes crossing at a shrinking angle. They print the quantities that matter for an immersed discretization, where
 * the integration cells are quadrature domains only: cell count, smallest cell and subphase volume, collapsed cells,
 * and the number and smallest area of interface facets.
 *
 * The multi-rank sweeps are in UT_XTK_Cut_Pathologies_MPI.cpp; the helpers are in cl_XTK_Cut_Pathology_*.hpp.
 *
 */

#include "catch.hpp"

#include "cl_XTK_Cut_Pathology_Harness.hpp"

namespace moris::xtk
{
    using namespace cut_pathology_test;

    TEST_CASE( "XTK pathology sweep: single hex, corner signs", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_SIGNS]" )
    {
        if ( par_size() == 1 )
        {
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
            {
                Sweep_Summary tSummary = single_geometry_sweep( "signs", { -1.0, 0.0, 1.0 }, iRegime );
                print_summary( "single hex, corner values {-1,0,1}, " + regime_name( iRegime ), tSummary );
                check_invariants( tSummary, iRegime != Regime::MULTILINEAR, false );
            }
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, near-zero corner values", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_NEAR_ZERO]" )
    {
        if ( par_size() == 1 )
        {
            // 1e-13 lies inside the isocontour tolerance (1e-12, values label as INTERFACE), 1e-11 just outside it
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
            {
                for ( real tDelta : { 1.0e-13, 1.0e-11 } )
                {
                    Sweep_Summary      tSummary = single_geometry_sweep( "near_zero", { -1.0, -tDelta, tDelta, 1.0 }, iRegime );
                    std::ostringstream tTitle;
                    tTitle << "single hex, corner values {-1,-d,+d,1}, d = " << tDelta << ", " << regime_name( iRegime );
                    print_summary( tTitle.str(), tSummary );
                    check_invariants( tSummary, iRegime != Regime::MULTILINEAR, false );
                }
            }
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, two geometries in both orders", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_TWO_GEOMETRIES]" )
    {
        if ( par_size() == 1 )
        {
            for ( Regime iRegime : { Regime::STORED, Regime::ANALYTIC, Regime::MULTILINEAR } )
            {
                Sweep_Summary tSummary = two_geometry_sweep( "two_geometries", { -1.0, 0.0, 1.0 }, iRegime );
                print_summary( "single hex, geometry A corner values {-1,0,1} x 8 second geometries, both orders, " + regime_name( iRegime ), tSummary );
                check_invariants( tSummary, iRegime == Regime::STORED, iRegime == Regime::STORED );
            }
        }
    }

    // one case of the two-geometry sweep in a fresh process, for triage:
    //   MORIS_XTK_PATHOLOGY_ONE="<stored|analytic|multilinear>|<second geometry name>|<AB|BA>|c0;c1;...;c7"
    TEST_CASE( "XTK pathology: one two-geometry case from the environment", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_ONE]" )
    {
        if ( par_size() == 1 )
        {
            const char* tSpec = std::getenv( "MORIS_XTK_PATHOLOGY_ONE" );
            REQUIRE( tSpec != nullptr );

            std::stringstream          tStream( tSpec );
            std::string                tRegimeName, tSecondName, tOrder, tCornerList;
            std::getline( tStream, tRegimeName, '|' );
            std::getline( tStream, tSecondName, '|' );
            std::getline( tStream, tOrder, '|' );
            std::getline( tStream, tCornerList, '|' );

            Regime tRegime = tRegimeName == "stored" ? Regime::STORED : tRegimeName == "analytic" ? Regime::ANALYTIC : Regime::MULTILINEAR;

            std::stringstream tCornerStream( tCornerList );
            std::string       tValue;
            for ( uint iCorner = 0; iCorner < 8; iCorner++ )
            {
                REQUIRE( std::getline( tCornerStream, tValue, ';' ) );
                gCornersA[ iCorner ] = std::stod( tValue );
            }

            bool tFound = false;
            for ( const std::pair< std::string, Corner_Values >& iSecond : second_geometries() )
            {
                if ( iSecond.first == tSecondName )
                {
                    gCornersB = iSecond.second;
                    tFound    = true;
                }
            }
            REQUIRE( tFound );

            Unit_Hex tHex;
            Vector< std::shared_ptr< gen::Geometry > > tGeometries =
                    tOrder == "AB" ? Vector< std::shared_ptr< gen::Geometry > >{ create_geometry( &field_a, tRegime ), create_geometry( &field_b, tRegime ) }
                                   : Vector< std::shared_ptr< gen::Geometry > >{ create_geometry( &field_b, tRegime ), create_geometry( &field_a, tRegime ) };
            Case_Result tResult = decompose_case( tHex.mMesh, tGeometries );

            std::cout << "case " << tSpec << ": threw " << tResult.mThrew << " (" << tResult.mMessage.substr( 0, 160 ) << "), cells "
                      << tResult.mNumCells << ", total volume " << tResult.mTotalVolume << ", volume<=0 " << tResult.mNumNonPositive
                      << ", mixed " << tResult.mNumMixedSign << ", bad facets " << tResult.mNumBadFacets << "\n";
            CHECK( not tResult.mThrew );
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, 1 to 6 planes accumulating", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_LADDER]" )
    {
        if ( par_size() == 1 )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            // six generic planes through points near the cell centre (off the grid and off the subdivision vertices)
            const Point tNormals[ 6 ] = { { 1.0, 0.31, 0.17 }, { 0.23, 1.0, -0.41 }, { -0.37, 0.19, 1.0 },
                { 1.0, -0.71, 0.53 }, { 0.61, 0.83, -0.29 }, { -0.47, 1.0, 0.67 } };
            const Point tPoints[ 6 ] = { { 0.47, 0.52, 0.49 }, { 0.53, 0.46, 0.51 }, { 0.50, 0.55, 0.45 },
                { 0.44, 0.50, 0.56 }, { 0.56, 0.48, 0.47 }, { 0.49, 0.44, 0.53 } };

            Sweep_Summary tSummary;
            print_metrics_header( "planes accumulating in one hex (first k of 6, then the same k in reverse order)" );
            for ( uint tK = 1; tK <= 6; tK++ )
            {
                Vector< Plane > tForward;
                for ( uint iPlane = 0; iPlane < tK; iPlane++ )
                {
                    tForward.push_back( plane_through( tNormals[ iPlane ], tPoints[ iPlane ], 0.0 ) );
                }
                Vector< Plane > tReverse;
                for ( uint iPlane = tK; iPlane > 0; iPlane-- )
                {
                    tReverse.push_back( tForward( iPlane - 1 ) );
                }

                Case_Result tForwardResult = decompose_planes( tMesh, tForward );
                Case_Result tReverseResult = decompose_planes( tMesh, tReverse );
                accumulate_invariants( tSummary, tForwardResult, {} );
                accumulate_invariants( tSummary, tReverseResult, {} );
                print_metrics_row( "k = " + std::to_string( tK ) + " forward", tForwardResult );
                print_metrics_row( "k = " + std::to_string( tK ) + " reverse", tReverseResult );
            }
            print_summary( "planes accumulating: invariants", tSummary );
            check_invariants( tSummary, false, false );
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, two parallel planes approaching", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_PARALLEL]" )
    {
        if ( par_size() == 1 )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            // A = n.(x - p), B = A - delta: the layer { A >= 0, B < 0 } has thickness delta (h = 1)
            const Point tNormal = { 1.0, 0.3, 0.2 };
            const Point tPoint  = { 0.45, 0.5, 0.5 };

            Sweep_Summary tSummary;
            print_metrics_header( "two parallel planes a distance delta apart", "     layer V  reference V" );
            for ( real tDelta : { 1.0e-1, 1.0e-3, 1.0e-6, 1.0e-9, 1.0e-11, 1.0e-12, 1.0e-13, 1.0e-14 } )
            {
                Plane tA = plane_through( tNormal, tPoint, 0.0 );
                Plane tB = plane_through( tNormal, tPoint, tDelta );

                Case_Result tResult = decompose_planes( tMesh, { tA, tB } );
                accumulate_invariants( tSummary, tResult, {} );

                // layer: A positive (bit 0 set), B negative (bit 1 clear)
                real tReference = reference_negative_volume( plane_corner_values( tB ) ) - reference_negative_volume( plane_corner_values( tA ) );
                std::ostringstream tExtra;
                tExtra << std::setw( 13 ) << std::setprecision( 3 ) << ( tResult.mThrew ? 0.0 : tResult.mSignVolume( 1 ) )
                       << std::setw( 13 ) << tReference;
                std::ostringstream tLabel;
                tLabel << "delta = " << tDelta;
                print_metrics_row( tLabel.str(), tResult, tExtra.str() );
            }
            print_summary( "parallel planes: invariants", tSummary );
            check_invariants( tSummary, false, false );
        }
    }

    TEST_CASE( "XTK pathology sweep: single hex, two planes crossing at a shrinking angle", "[.][XTK][XTK_PATHOLOGY][XTK_PATHOLOGY_HEX_ANGLE]" )
    {
        if ( par_size() == 1 )
        {
            Unit_Hex                     tHex;
            hmr::Interpolation_Mesh_HMR* tMesh = tHex.mMesh;

            // two vertical planes through a point off the subdivision vertices, crossing at angle theta
            const Point tPoint = { 0.513, 0.507, 0.5 };
            const real  tAlpha = 0.3;

            Sweep_Summary tSummary;
            print_metrics_header( "two planes crossing at angle theta" );
            for ( real tThetaDegrees : { 90.0, 30.0, 5.0, 1.0, 0.1, 0.01 } )
            {
                real  tTheta = tThetaDegrees * M_PI / 180.0;
                Plane tA     = plane_through( { std::cos( tAlpha ), std::sin( tAlpha ), 0.0 }, tPoint, 0.0 );
                Plane tB     = plane_through( { std::cos( tAlpha + tTheta ), std::sin( tAlpha + tTheta ), 0.0 }, tPoint, 0.0 );

                Case_Result tResult = decompose_planes( tMesh, { tA, tB } );
                accumulate_invariants( tSummary, tResult, {} );

                std::ostringstream tLabel;
                tLabel << "theta = " << tThetaDegrees << " deg";
                print_metrics_row( tLabel.str(), tResult );
            }
            print_summary( "crossing planes: invariants", tSummary );
            check_invariants( tSummary, false, false );
        }
    }

}    // namespace moris::xtk
