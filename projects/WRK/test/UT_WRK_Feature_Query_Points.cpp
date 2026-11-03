/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Feature_Query_Points.cpp
 *
 * Geometry-aware query points (GEN level-set keys feature_query_points, feature_query_check, feature_query_refine):
 * a plate x < 0.3122 joined by a disk whose outline meets the plate face at a thin void wedge (opening 18 degrees, as
 * the ring-plate junctions of the enrichment fixture). The plate face lies 0.0003 left of the node column x = 5/16, so
 * every Lagrange cell up to level 3 that holds the wedge tip welds it (both ends of its horizontal edges are solid); a
 * void point 0.002 into the wedge is represented as solid until level 4. The sub-cell feature detector does not flag
 * the wedge on the coarse cell (its lattice misses it, or rates it below the depth tolerance) and stops at once; the
 * query point keeps the refinement going until the cut represents it.
 *
 */

#include "catch.hpp"

#include <cmath>

#include "cl_Communication_Tools.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Parameters.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_MTK_Cell.hpp"
#include "cl_GEN_Geometry_Engine.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_User_Defined_Field.hpp"
#include "cl_GEN_Cell_Feature_Detector.hpp"
#include "cl_WRK_GEN_Performer.hpp"
#include "cl_WRK_perform_refinement.hpp"
#include "cl_XTK_Model.hpp"
#include "cl_XTK_Cut_Integration_Mesh.hpp"
#include "fn_PRM_GEN_Parameters.hpp"

namespace moris::wrk
{
    namespace
    {
        const real gPlate  = 0.3122;                     // plate face (0.0003 left of the node column 5/16)
        const real gRadius = 0.33;                       // disk radius
        const real gCentre = gPlate + 0.95 * gRadius;    // disk centre x (y = 0.5): the outline meets the face at 18 degrees

        real
        plate_and_disk( const Matrix< DDRMat >& aPoint, const Vector< real >& )
        {
            return std::max( gRadius - std::hypot( aPoint( 0 ) - gCentre, aPoint( 1 ) - 0.5 ), gPlate - aPoint( 0 ) );
        }

        // void point 0.002 into the lower wedge along its bisector, and a solid point in the plate next to the tip
        std::string
        query_points()
        {
            real tTipY = 0.5 - gRadius * std::sqrt( 1.0 - 0.95 * 0.95 );
            real tBx   = 0.31225 * 0.5 / std::hypot( 0.31225 * 0.5, 0.5 * ( 1.0 + 0.95 ) );    // bisector of (0, -1) and the outline tangent
            real tBy   = -std::sqrt( 1.0 - tBx * tBx );
            char tString[ 256 ];
            std::snprintf( tString, sizeof( tString ), "%.10f,%.10f:-1; %.10f,%.10f:1", gPlate + 0.002 * tBx, tTipY + 0.002 * tBy, gPlate - 0.002, tTipY - 0.002 );
            return tString;
        }

        struct Run_Result
        {
            uint mNumCells    = 0;      // active Lagrange cells after the refinement
            real mFingerprint = 0.0;    // sum over the cells of level x (centre x + 2 centre y): equal meshes give equal sums
            uint mMaxLevel    = 0;
            uint mUnresolved  = 0;    // query points the refinement left on the wrong side (WRK)
            uint mMismatches  = 0;    // query points the cut represents on the wrong side (XTK check)
            uint mLocated     = 0;
        };

        gen::Level_Set_Parameters
        level_set_parameters( const std::string& aPoints, bool aRefine, uint aCheck, uint aPasses )
        {
            Parameter_List tList = prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED );
            tList.set( "use_analytic_intersection", true );
            tList.set( "number_of_refinements", aPasses );
            tList.set( "refinement_mesh_index", 0u );
            tList.set( "feature_refinement_samples", 8u );
            tList.set( "feature_refinement_lipschitz", 1.0 );
            tList.set( "feature_refinement_depth_tolerance", 0.02 );
            tList.set( "feature_query_points", aPoints );
            tList.set( "feature_query_check", aCheck );
            tList.set( "feature_query_refine", aRefine );
            return gen::Level_Set_Parameters( tList );
        }

        hmr::Parameters
        hmr_parameters()
        {
            hmr::Parameters tParameters;
            tParameters.set_number_of_elements_per_dimension( { 1, 1 } );
            tParameters.set_domain_dimensions( { 1.0, 1.0 } );
            tParameters.set_domain_offset( { 0.0, 0.0 } );
            tParameters.set_bspline_truncation( true );
            tParameters.set_lagrange_orders( { 1 } );
            tParameters.set_lagrange_patterns( { 0 } );
            tParameters.set_bspline_orders( { 1 } );
            tParameters.set_bspline_patterns( { 1 } );
            tParameters.set_refinement_buffer( 0 );
            tParameters.set_staircase_buffer( 1 );
            return tParameters;
        }

        /**
         * GEN feature refinement of the Lagrange mesh (WRK), then, if requested, the XTK cut and its query-point check.
         */
        Run_Result
        refine_and_cut( const std::string& aPoints, bool aRefine, uint aCheck, uint aPasses, bool aCut )
        {
            Run_Result tResult;

            // GEN for the refinement. The workflow builds it without a mesh (every node then counts as a background node
            // and the analytic field is evaluated at its coordinates); the mesh constructor used here needs a mesh, so it
            // gets a uniform 32 x 32 copy, which has more nodes than any refined mesh of this test (same effect)
            hmr::Parameters tCoarseParameters = hmr_parameters();
            tCoarseParameters.set_number_of_elements_per_dimension( { 32, 32 } );
            hmr::HMR tCoarseHMR( tCoarseParameters );
            tCoarseHMR.finalize();
            hmr::Interpolation_Mesh_HMR*    tCoarseMesh = tCoarseHMR.create_interpolation_mesh( 0 );
            gen::Geometry_Engine_Parameters tGenParameters;
            tGenParameters.mGeometries                                = { std::make_shared< gen::Level_Set_Geometry >(
                    std::make_shared< gen::User_Defined_Field >( &plate_and_disk ), level_set_parameters( aPoints, aRefine, aCheck, aPasses ) ) };
            std::shared_ptr< gen::Geometry_Engine > tRefinementEngine = std::make_shared< gen::Geometry_Engine >( tCoarseMesh, tGenParameters );

            hmr::Parameters             tParameters = hmr_parameters();
            std::shared_ptr< hmr::HMR > tHMR        = std::make_shared< hmr::HMR >( tParameters );
            tHMR->perform_initial_refinement();
            Refinement_Mini_Performer tRefinement;
            tRefinement.perform_refinement_old( tHMR, { std::make_shared< Gen_Performer >( tRefinementEngine ) } );
            tResult.mUnresolved = tRefinement.get_num_unresolved_query_points();
            tHMR->finalize();

            hmr::Interpolation_Mesh_HMR* tMesh = tHMR->create_interpolation_mesh( 0 );
            tResult.mNumCells                  = tMesh->get_num_elems();
            for ( uint iCell = 0; iCell < tResult.mNumCells; iCell++ )
            {
                const mtk::Cell& tCell   = tMesh->get_mtk_cell( (moris_index)iCell );
                Matrix< DDRMat > tCoords = tCell.get_vertex_coords();
                real             tX      = 0.0;
                real             tY      = 0.0;
                for ( uint iVertex = 0; iVertex < tCoords.n_rows(); iVertex++ )
                {
                    tX += tCoords( iVertex, 0 ) / tCoords.n_rows();
                    tY += tCoords( iVertex, 1 ) / tCoords.n_rows();
                }
                tResult.mFingerprint += real( tCell.get_level() + 1 ) * ( tX + 2.0 * tY );
                tResult.mMaxLevel = std::max( tResult.mMaxLevel, uint( tCell.get_level() ) );
            }

            if ( aCut )
            {
                gen::Geometry_Engine_Parameters tCutParameters;
                tCutParameters.mGeometries = { std::make_shared< gen::Level_Set_Geometry >(
                        std::make_shared< gen::User_Defined_Field >( &plate_and_disk ), level_set_parameters( aPoints, aRefine, aCheck, aPasses ) ) };
                gen::Geometry_Engine tEngine( tMesh, tCutParameters );
                xtk::Model           tModel( 2, tMesh, &tEngine );
                tModel.mVerbose = false;
                tModel.decompose( { Subdivision_Method::NC_REGULAR_SUBDIVISION_QUAD4, Subdivision_Method::C_TRI3 } );

                xtk::Cut_Integration_Mesh* tCutMesh = tModel.get_cut_integration_mesh();
                if ( tCutMesh->mFeatureQueryMismatches.size() > 0 )
                {
                    tResult.mMismatches = tCutMesh->mFeatureQueryMismatches( 0 );
                    tResult.mLocated    = tCutMesh->mFeatureQueryLocated( 0 );
                }
            }

            delete tMesh;
            delete tCoarseMesh;
            return tResult;
        }
    }    // namespace

    TEST_CASE( "Feature query points: parser and represented side", "[WRK][FEATURE_QUERY]" )
    {
        // syntax: points separated by ';', "x,y[,z][:phase[:radius]]"
        Vector< gen::Feature_Query_Point > tPoints = gen::parse_feature_query_points( " 0.1, 0.2 ; 0.3,0.4,0.5:-1:0.02; 1,2:1 ;" );
        REQUIRE( tPoints.size() == 3 );
        CHECK( tPoints( 0 ).mCoordinates.numel() == 2 );
        CHECK( tPoints( 0 ).mPhase == 0 );
        CHECK( tPoints( 0 ).mRadius == 0.0 );
        CHECK( tPoints( 1 ).mCoordinates.numel() == 3 );
        CHECK( tPoints( 1 ).mCoordinates( 2 ) == Approx( 0.5 ) );
        CHECK( tPoints( 1 ).mPhase == -1 );
        CHECK( tPoints( 1 ).mRadius == Approx( 0.02 ) );
        CHECK( tPoints( 2 ).mPhase == 1 );
        CHECK( gen::parse_feature_query_points( "" ).size() == 0 );
        CHECK_THROWS( gen::parse_feature_query_points( "0.1,0.2:2" ) );
        CHECK_THROWS( gen::parse_feature_query_points( "0.1" ) );
        CHECK_THROWS( gen::parse_feature_query_points( "0.1,0.2:1:-1" ) );
        CHECK_THROWS( gen::parse_feature_query_points( "0.1,abc" ) );

        // rows of a feature_query_function: x, y, phase, radius
        Matrix< DDRMat > tRows = { { 0.1, 0.2, -1.0, 0.0 }, { 0.3, 0.4, 1.0, 0.5 } };
        CHECK( gen::feature_query_points_from_rows( tRows ).size() == 2 );
        CHECK( gen::feature_query_points_from_rows( tRows )( 1 ).mRadius == Approx( 0.5 ) );

        // represented side at the void point in the cell holding the tip: welded (solid) at level 3, void at level 4
        Vector< gen::Feature_Query_Point > tQuery = gen::parse_feature_query_points( query_points() );
        gen::Point_Field_Function          tPhi   = []( const Matrix< DDRMat >& aPoint ) -> real { return plate_and_disk( aPoint, {} ); };
        const Matrix< DDRMat >&            tVoid  = tQuery( 0 ).mCoordinates;
        CHECK( tPhi( tVoid ) < 0.0 );
        for ( uint iLevel = 0; iLevel <= 5; iLevel++ )
        {
            real             tH       = 1.0 / std::pow( 2.0, iLevel );
            real             tX       = std::floor( tVoid( 0 ) / tH ) * tH;
            real             tY       = std::floor( tVoid( 1 ) / tH ) * tH;
            Matrix< DDRMat > tCorners = { { tX, tY }, { tX + tH, tY }, { tX + tH, tY + tH }, { tX, tY + tH } };
            CHECK( gen::cell_contains_point( tCorners, tVoid ) );
            CHECK( gen::represented_region_at_point( tPhi, tCorners, 1.0e-12, tVoid, true ) == ( iLevel < 4 ? 1 : -1 ) );
        }
    }

    TEST_CASE( "Feature query points: check and refinement of a welded wedge tip", "[WRK][FEATURE_QUERY]" )
    {
        if ( par_size() == 1 )
        {
            const uint tPasses = 6;

            // detector only: it stops below level 4 (it does not flag the wedge), the cut welds the tip
            Run_Result tNoPoints  = refine_and_cut( "", false, 1, tPasses, false );
            Run_Result tCheckOnly = refine_and_cut( query_points(), false, 1, tPasses, true );
            CHECK( tCheckOnly.mMaxLevel < 4 );
            CHECK( tCheckOnly.mLocated == 2 );
            CHECK( tCheckOnly.mMismatches == 1 );
            CHECK( tCheckOnly.mUnresolved == 0 );    // refinement on the points is off: nothing tracked

            // default off: query points without feature_query_refine leave the refinement bit-identical
            CHECK( tCheckOnly.mNumCells == tNoPoints.mNumCells );
            CHECK( tCheckOnly.mFingerprint == tNoPoints.mFingerprint );

            // refinement on the points: the cut represents the void point; one level more than the detector alone
            Run_Result tRefined = refine_and_cut( query_points(), true, 1, tPasses, true );
            CHECK( tRefined.mMaxLevel == 4 );
            CHECK( tRefined.mMismatches == 0 );
            CHECK( tRefined.mUnresolved == 0 );
            CHECK( tRefined.mNumCells > tCheckOnly.mNumCells );

            // number_of_refinements too small for the point: logged UNRESOLVED, the cut still welds it
            Run_Result tCapped = refine_and_cut( query_points(), true, 1, 3, true );
            CHECK( tCapped.mUnresolved == 1 );
            CHECK( tCapped.mMismatches == 1 );

            // check off: no check results
            Run_Result tOff = refine_and_cut( query_points(), false, 0, tPasses, true );
            CHECK( tOff.mLocated == 0 );
            CHECK( tOff.mMismatches == 0 );

            // error mode: a mismatch aborts the cut
            CHECK_THROWS( refine_and_cut( query_points(), false, 2, tPasses, true ) );
            CHECK_NOTHROW( refine_and_cut( query_points(), true, 2, tPasses, true ) );
        }
    }
}    // namespace moris::wrk
