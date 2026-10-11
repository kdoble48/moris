/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Reinitialize_Preserve.cpp
 *
 * MB-027: the in-place SDF redistance is not idempotent; applied to a fixed design it moves the cut, so the criteria
 * drift with every application. sdf_preserve_interface (Reinitialize_Performer::freeze_interface_band) keeps every
 * design coefficient that can reach a possibly cut cell. This test applies a redistance-like operator (sign-preserving
 * nodal distance targets, L2 projection onto the B-spline basis, as in compute_sdf_from_interface) eight times to a
 * curved (elliptic, not a distance field) interface at B-spline orders 1 and 2 and checks that with the freeze band the
 * cut cells, every nodal value on them and every edge intersection point stay bit-identical, while the far field is
 * redistanced; and that without it the cut moves (negative control).
 */

#include <cmath>
#include <map>
#include <utility>

#include "catch.hpp"
#include "cl_Communication_Tools.hpp"
#include "cl_HMR.hpp"
#include "cl_HMR_Mesh_Interpolation.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Mesh_Factory.hpp"
#include "cl_MTK_Mesh_Pair.hpp"
#include "cl_MTK_Interpolation_Mesh.hpp"
#include "cl_MTK_Integration_Mesh.hpp"
#include "cl_MTK_Field_Discrete.hpp"
#include "cl_MTK_Mapper.hpp"
#include "cl_MTK_Vertex.hpp"
#include "fn_PRM_HMR_Parameters.hpp"
#include "cl_WRK_Reinitialize_Performer.hpp"

namespace moris::wrk
{
    namespace preserve_interface_test
    {
        // ellipse level set, deliberately not a distance field (|grad| varies from ~0.7 to ~10)
        real
        ellipse( real aX, real aY )
        {
            return 0.3 * ( std::pow( aX / 0.55, 2 ) + std::pow( aY / 0.4, 2 ) - 1.0 );
        }

        // first-order distance to the ellipse, |f| / |grad f|
        real
        distance( real aX, real aY )
        {
            real tGx = 0.3 * 2.0 * aX / ( 0.55 * 0.55 );
            real tGy = 0.3 * 2.0 * aY / ( 0.4 * 0.4 );
            return std::abs( ellipse( aX, aY ) ) / std::max( std::sqrt( tGx * tGx + tGy * tGy ), 1.0e-3 );
        }

        /** L2 projection of nodal values onto the B-spline coefficients (the step compute_sdf_from_interface uses) */
        Matrix< DDRMat >
        project( const mtk::Mesh_Pair& aMeshPair, const Matrix< DDRMat >& aNodalValues )
        {
            mtk::Field_Discrete tField( aMeshPair, 0 );
            tField.unlock_field();
            tField.set_values( aNodalValues );
            mtk::Mapper tMapper;
            tField.unlock_field();
            tMapper.map_input_field_to_output_field_2( &tField );
            return tField.get_coefficients();
        }

        /** cut geometry of a coefficient vector: nodal values, cut cells, corner-edge intersection points */
        struct Cut
        {
            Matrix< DDRMat >                                         mNodal;
            std::vector< bool >                                      mCutCell;
            std::map< std::pair< uint, uint >, std::pair< real, real > > mPoints;    // (cell, edge) -> point
        };

        Cut
        cut( mtk::Interpolation_Mesh* aMesh, const Matrix< DDRMat >& aCoeffs )
        {
            Cut tCut;
            tCut.mNodal = Reinitialize_Performer::nodal_values( aMesh, 0, aCoeffs );
            uint tNumCells = aMesh->get_num_elems();
            tCut.mCutCell.assign( tNumCells, false );
            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                Matrix< IndexMat > tNodes = aMesh->get_nodes_connected_to_element_loc_inds( (moris_index)iCell );
                uint               tPos   = 0;
                for ( uint iNode = 0; iNode < tNodes.numel(); iNode++ )
                {
                    tPos += tCut.mNodal( tNodes( iNode ) ) >= 0.0 ? 1 : 0;
                }
                tCut.mCutCell[ iCell ] = tPos > 0 && tPos < tNodes.numel();

                // corner edges (corners are the first four nodes of QUAD4 and QUAD9)
                for ( uint iEdge = 0; iEdge < 4; iEdge++ )
                {
                    moris_index tA  = tNodes( iEdge );
                    moris_index tB  = tNodes( ( iEdge + 1 ) % 4 );
                    real        tVa = tCut.mNodal( tA );
                    real        tVb = tCut.mNodal( tB );
                    if ( ( tVa >= 0.0 ) != ( tVb >= 0.0 ) )
                    {
                        Matrix< DDRMat > tXa = aMesh->get_mtk_vertex( tA ).get_coords();
                        Matrix< DDRMat > tXb = aMesh->get_mtk_vertex( tB ).get_coords();
                        real             tS  = tVa / ( tVa - tVb );
                        tCut.mPoints[ { iCell, iEdge } ] = { tXa( 0 ) + tS * ( tXb( 0 ) - tXa( 0 ) ), tXa( 1 ) + tS * ( tXb( 1 ) - tXa( 1 ) ) };
                    }
                }
            }
            return tCut;
        }

        struct Result
        {
            uint mCutCells       = 0;
            uint mPoints         = 0;
            bool mCutCellsEqual  = true;
            bool mNodalEqual     = true;    // nodal values on every node of every cut cell, bitwise
            bool mPointsEqual    = true;    // every intersection point, bitwise
            real mMaxPointMove   = 0.0;
            real mMaxFarChange   = 0.0;     // max |c - psi| over all coefficients after the redistances
            uint mFrozen         = 0;
            uint mNumCoeffs      = 0;
        };

        Result
        run( uint aOrder, bool aPreserve, uint aApplications )
        {
            Module_Parameter_Lists tParameters( Module_Type::HMR );
            tParameters.set( "number_of_elements_per_dimension", 24, 24 );
            tParameters.set( "domain_dimensions", 2.0, 2.0 );
            tParameters.set( "domain_offset", -1.0, -1.0 );
            tParameters.set( "lagrange_output_meshes", std::string( "0" ) );
            tParameters.set( "lagrange_orders", std::to_string( aOrder ) );
            tParameters.set( "lagrange_pattern", std::string( "0" ) );
            tParameters.set( "bspline_orders", std::to_string( aOrder ) );
            tParameters.set( "bspline_pattern", std::string( "0" ) );
            tParameters.set( "lagrange_to_bspline", "0" );
            tParameters.set( "refinement_buffer", 0 );
            tParameters.set( "staircase_buffer", 1 );
            tParameters.set( "use_number_aura", false );
            tParameters.set( "severity_level", 2 );

            hmr::HMR tHMR( tParameters );
            tHMR.perform_initial_refinement();
            tHMR.finalize();

            hmr::Interpolation_Mesh_HMR* tIPMesh = tHMR.create_interpolation_mesh( 0 );
            mtk::Integration_Mesh*       tIGMesh = mtk::create_integration_mesh_from_interpolation_mesh( mtk::MeshType::HMR, tIPMesh );
            mtk::Mesh_Pair               tMeshPair( tIPMesh, tIGMesh, true );

            uint             tNumNodes = tIPMesh->get_num_nodes();
            Matrix< DDRMat > tNodal( tNumNodes, 1, 0.0 );
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tX = tIPMesh->get_mtk_vertex( iNode ).get_coords();
                tNodal( iNode )     = ellipse( tX( 0 ), tX( 1 ) );
            }
            const Matrix< DDRMat > tPsi    = project( tMeshPair, tNodal );
            const Cut              tCutPsi = cut( tIPMesh, tPsi );

            Vector< Vector< moris_index > > tCellCoefficients = Reinitialize_Performer::cell_coefficients( tIPMesh, 0 );

            Result           tResult;
            Matrix< DDRMat > tCoeffs = tPsi;
            for ( uint iApply = 0; iApply < aApplications; iApply++ )
            {
                // targets: the sign of the current field, the distance to the ellipse (as compute_sdf_from_interface)
                Matrix< DDRMat > tCurrent = Reinitialize_Performer::nodal_values( tIPMesh, 0, tCoeffs );
                Matrix< DDRMat > tTargets( tNumNodes, 1, 0.0 );
                for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                {
                    Matrix< DDRMat > tX = tIPMesh->get_mtk_vertex( iNode ).get_coords();
                    tTargets( iNode )   = ( tCurrent( iNode ) >= 0.0 ? 1.0 : -1.0 ) * distance( tX( 0 ), tX( 1 ) );
                }
                Matrix< DDRMat > tNew = project( tMeshPair, tTargets );

                if ( aPreserve )
                {
                    uint tRings     = 0;
                    tResult.mFrozen = Reinitialize_Performer::freeze_interface_band( tNew, tCoeffs, tCellCoefficients, 1.0e-12, tRings );
                }
                tCoeffs = tNew;
            }

            const Cut tCutPhi  = cut( tIPMesh, tCoeffs );
            tResult.mNumCoeffs = tPsi.numel();
            for ( uint iC = 0; iC < tPsi.numel(); iC++ )
            {
                tResult.mMaxFarChange = std::max( tResult.mMaxFarChange, std::abs( tCoeffs( iC ) - tPsi( iC ) ) );
            }
            for ( uint iCell = 0; iCell < tCutPsi.mCutCell.size(); iCell++ )
            {
                tResult.mCutCells += tCutPsi.mCutCell[ iCell ] ? 1 : 0;
                tResult.mCutCellsEqual = tResult.mCutCellsEqual && tCutPsi.mCutCell[ iCell ] == tCutPhi.mCutCell[ iCell ];
                if ( tCutPsi.mCutCell[ iCell ] )
                {
                    Matrix< IndexMat > tNodes = tIPMesh->get_nodes_connected_to_element_loc_inds( (moris_index)iCell );
                    for ( uint iNode = 0; iNode < tNodes.numel(); iNode++ )
                    {
                        tResult.mNodalEqual = tResult.mNodalEqual && tCutPsi.mNodal( tNodes( iNode ) ) == tCutPhi.mNodal( tNodes( iNode ) );
                    }
                }
            }
            tResult.mPoints      = tCutPsi.mPoints.size();
            tResult.mPointsEqual = tCutPsi.mPoints.size() == tCutPhi.mPoints.size();
            for ( const std::pair< const std::pair< uint, uint >, std::pair< real, real > >& tEntry : tCutPsi.mPoints )
            {
                std::map< std::pair< uint, uint >, std::pair< real, real > >::const_iterator tFound = tCutPhi.mPoints.find( tEntry.first );
                if ( tFound == tCutPhi.mPoints.end() )
                {
                    tResult.mPointsEqual = false;
                    continue;
                }
                tResult.mPointsEqual  = tResult.mPointsEqual && tFound->second == tEntry.second;
                tResult.mMaxPointMove = std::max( tResult.mMaxPointMove,
                        std::hypot( tFound->second.first - tEntry.second.first, tFound->second.second - tEntry.second.second ) );
            }

            // the mesh pair owns (and deletes) both meshes
            return tResult;
        }
    }    // namespace preserve_interface_test

    TEST_CASE( "Interface-preserving SDF redistance keeps the cut bit-identical", "[WRK][sdf_preserve_interface][MB-027]" )
    {
        if ( par_size() == 1 )
        {
            for ( uint tOrder = 1; tOrder <= 2; tOrder++ )
            {
                preserve_interface_test::Result tPreserved = preserve_interface_test::run( tOrder, true, 8 );
                preserve_interface_test::Result tPlain     = preserve_interface_test::run( tOrder, false, 8 );

                MORIS_LOG_INFO( "order %u: %u cut cells, %u edge points; preserved: %u of %u coefficients frozen, far-field change %.3e; "
                                "plain: max point move %.3e",
                        tOrder, tPreserved.mCutCells, tPreserved.mPoints, tPreserved.mFrozen, tPreserved.mNumCoeffs,
                        tPreserved.mMaxFarChange, tPlain.mMaxPointMove );

                REQUIRE( tPreserved.mCutCells > 20 );

                // with the freeze band, eight redistances leave the cut exactly as it was
                CHECK( tPreserved.mCutCellsEqual );
                CHECK( tPreserved.mNodalEqual );
                CHECK( tPreserved.mPointsEqual );
                CHECK( tPreserved.mMaxPointMove == 0.0 );

                // ... while the rest of the field is redistanced
                CHECK( tPreserved.mFrozen < tPreserved.mNumCoeffs );
                CHECK( tPreserved.mMaxFarChange > 1.0e-3 );

                // negative control (MB-027): the plain redistance moves the cut
                CHECK_FALSE( tPlain.mNodalEqual );
                CHECK( tPlain.mMaxPointMove > 0.0 );
            }
        }
    }
}    // namespace moris::wrk

namespace moris::wrk
{
    TEST_CASE( "SDF taper: wean tolerance with patience over accepted iterates", "[WRK][sdf_taper]" )
    {
        SDF_Wean_Monitor tMonitor;
        tMonitor.mTol      = 1.0e-3;
        tMonitor.mPatience = 3;

        // two small changes, then a large one resets the run
        CHECK_FALSE( tMonitor.update( 5.0e-4 ) );
        CHECK_FALSE( tMonitor.update( 5.0e-4 ) );
        CHECK_FALSE( tMonitor.update( 2.0e-3 ) );
        CHECK( tMonitor.mBelow == 0 );

        // three consecutive small changes wean it off, on the third
        CHECK_FALSE( tMonitor.update( 9.0e-4 ) );
        CHECK_FALSE( tMonitor.update( 1.0e-6 ) );
        CHECK( tMonitor.update( 0.0 ) );
        CHECK( tMonitor.mAccepted == 6 );

        // sticky
        CHECK( tMonitor.update( 1.0 ) );
        CHECK( tMonitor.mAccepted == 6 );

        // disabled tolerance never weans
        SDF_Wean_Monitor tOff;
        for ( uint iIter = 0; iIter < 100; iIter++ )
        {
            CHECK_FALSE( tOff.update( 0.0 ) );
        }
    }

    TEST_CASE( "SDF taper: hard stop after sdf_reinit_stop_iter accepted iterates", "[WRK][sdf_taper]" )
    {
        SDF_Wean_Monitor tMonitor;
        tMonitor.mStopIter = 4;
        CHECK_FALSE( tMonitor.update( 1.0 ) );
        CHECK_FALSE( tMonitor.update( 1.0 ) );
        CHECK_FALSE( tMonitor.update( 1.0 ) );
        CHECK( tMonitor.update( 1.0 ) );
        CHECK( tMonitor.mReason.find( "sdf_reinit_stop_iter" ) != std::string::npos );
    }
}    // namespace moris::wrk
