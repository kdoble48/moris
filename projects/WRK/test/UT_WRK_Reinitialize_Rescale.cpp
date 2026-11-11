/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Reinitialize_Rescale.cpp
 *
 * sdf_preserve_interface_slope_tol: slope restoration inside the band that sdf_preserve_interface freezes
 * (Reinitialize_Performer::rescale_interface_band). The band coefficients are multiplied by positive factors whose
 * logarithm is limited to a spread of 4 atanh( tol ) per cell. The tests apply one interface-preserving redistance
 * (sign-preserving distance targets, freeze band, rescale) to analytic level sets and check
 *  - the invariants: no coefficient and no node changes sign, the cut cells are the same, and every intersection point
 *    on a cell edge moves by at most tol of the edge length (coefficients = nodal values, order 1);
 *  - what the slope does: a signed distance is left bit-identical, a uniformly flat or steep field is brought to the end
 *    of the target range with the crossings at round-off, a field whose slope varies along the interface is corrected
 *    only as far as the tolerance allows (a level set of order 1 is fixed by its interface up to one factor);
 *  - the old modes: with the tolerance at zero the frozen band is bit-identical, and the parameter defaults are off.
 */

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <map>
#include <utility>
#include <vector>

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
#include "fn_PRM_OPT_Parameters.hpp"
#include "cl_WRK_Reinitialize_Performer.hpp"

namespace moris::wrk
{
    namespace band_rescale_test
    {
        const real gH = 1.0 / 12.0;    // element size of the 24 x 24 mesh on [-1, 1]^2

        enum class Field
        {
            PLANE,                // signed distance to an oblique plane
            SPHERE,               // signed distance to a circle / sphere
            PLANE_FLAT,           // 0.2 x the plane distance: uniformly too flat
            SPHERE_STEEP,         // 3 x the sphere distance: uniformly too steep
            SPHERE_TANH,          // 0.25 w tanh( d / w ): flat at the interface, flatter away from it
            SPHERE_CUBED,         // d^3 / h^2: zero slope on the interface
            SPHERE_VARYING,       // ( 0.3 + 0.5 ( x + 1 ) ) d: slope 0.3 to 1.3 along the interface
            NODE_ON_SPHERE,       // SPHERE_VARYING with the node at ( 0.5, 0, 0 ) 1e-6 h off the interface
            NODE_ON_SPHERE_SDF    // the signed distance of that sphere
        };

        /** signed distance of the field's interface */
        real
        distance( Field aField, const Matrix< DDRMat >& aX )
        {
            if ( aField == Field::PLANE || aField == Field::PLANE_FLAT )
            {
                return 0.6 * aX( 0 ) + 0.8 * aX( 1 ) - 0.137;
            }
            real tRadiusSquared = 0.0;
            for ( uint iDim = 0; iDim < aX.numel(); iDim++ )
            {
                tRadiusSquared += aX( iDim ) * aX( iDim );
            }
            real tRadius = ( aField == Field::NODE_ON_SPHERE || aField == Field::NODE_ON_SPHERE_SDF ) ? 0.5 - 1.0e-6 * gH : 0.613;
            return std::sqrt( tRadiusSquared ) - tRadius;
        }

        real
        value( Field aField, const Matrix< DDRMat >& aX )
        {
            real tD = distance( aField, aX );
            switch ( aField )
            {
                case Field::PLANE_FLAT:
                    return 0.2 * tD;
                case Field::SPHERE_STEEP:
                    return 3.0 * tD;
                case Field::SPHERE_TANH:
                    return 0.25 * 0.2 * std::tanh( tD / 0.2 );
                case Field::SPHERE_CUBED:
                    return tD * tD * tD / ( gH * gH );
                case Field::SPHERE_VARYING:
                case Field::NODE_ON_SPHERE:
                    return ( 0.3 + 0.5 * ( aX( 0 ) + 1.0 ) ) * tD;
                default:
                    return tD;
            }
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

        /** corner edges of a QUAD or HEX cell as pairs of local corner indices */
        std::vector< std::pair< uint, uint > >
        corner_edges( uint aDim )
        {
            if ( aDim == 2 )
            {
                return { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 } };
            }
            return { { 0, 1 }, { 1, 2 }, { 2, 3 }, { 3, 0 }, { 4, 5 }, { 5, 6 }, { 6, 7 }, { 7, 4 }, { 0, 4 }, { 1, 5 }, { 2, 6 }, { 3, 7 } };
        }

        /** cut of a coefficient vector: nodal values, cut cells, intersection fraction of every cut corner edge */
        struct Cut
        {
            Matrix< DDRMat >                          mNodal;
            std::vector< bool >                       mCutCell;
            std::map< std::pair< uint, uint >, real > mFraction;    // (cell, edge) -> fraction from the edge's first node
        };

        Cut
        cut( mtk::Interpolation_Mesh* aMesh, const Matrix< DDRMat >& aCoeffs )
        {
            uint                                   tDim     = aMesh->get_spatial_dim();
            uint                                   tCorners = tDim == 2 ? 4 : 8;
            std::vector< std::pair< uint, uint > > tEdges   = corner_edges( tDim );

            Cut tCut;
            tCut.mNodal    = Reinitialize_Performer::nodal_values( aMesh, 0, aCoeffs );
            uint tNumCells = aMesh->get_num_elems();
            tCut.mCutCell.assign( tNumCells, false );
            for ( uint iCell = 0; iCell < tNumCells; iCell++ )
            {
                Matrix< IndexMat > tNodes = aMesh->get_nodes_connected_to_element_loc_inds( (moris_index)iCell );
                uint               tPos   = 0;
                for ( uint iNode = 0; iNode < tCorners; iNode++ )
                {
                    tPos += tCut.mNodal( tNodes( iNode ) ) >= 0.0 ? 1 : 0;
                }
                tCut.mCutCell[ iCell ] = tPos > 0 && tPos < tCorners;
                for ( uint iEdge = 0; iEdge < tEdges.size(); iEdge++ )
                {
                    real tVa = tCut.mNodal( tNodes( tEdges[ iEdge ].first ) );
                    real tVb = tCut.mNodal( tNodes( tEdges[ iEdge ].second ) );
                    if ( ( tVa >= 0.0 ) != ( tVb >= 0.0 ) )
                    {
                        tCut.mFraction[ { iCell, iEdge } ] = tVa / ( tVa - tVb );
                    }
                }
            }
            return tCut;
        }

        struct Result
        {
            uint                mCutCells         = 0;
            uint                mCutEdges         = 0;
            uint                mNodeSignFlips    = 0;
            uint                mCoeffSignFlips   = 0;
            bool                mCutCellsEqual    = true;
            bool                mCutEdgesEqual    = true;
            real                mMaxShift         = 0.0;                  // largest move of an intersection point, in edge lengths
            bool                mBandBitIdentical = true;                 // every coefficient equals the frozen-band result
            bool                mFarBitIdentical  = true;                 // every coefficient outside the band cells equals the frozen-band result
            real                mSlopeBefore[ 3 ] = { 0.0, 0.0, 0.0 };    // min, median, max over the cut cells
            real                mSlopeAfter[ 3 ]  = { 0.0, 0.0, 0.0 };
            Band_Rescale_Report mReport;
        };

        void
        slope_statistics( const Matrix< DDRMat >& aSlopes, const std::vector< bool >& aCutCell, real* aStatistics )
        {
            std::vector< real > tValues;
            for ( uint iCell = 0; iCell < aCutCell.size(); iCell++ )
            {
                if ( aCutCell[ iCell ] )
                {
                    tValues.push_back( aSlopes( iCell ) );
                }
            }
            std::sort( tValues.begin(), tValues.end() );
            aStatistics[ 0 ] = tValues.front();
            aStatistics[ 1 ] = tValues[ tValues.size() / 2 ];
            aStatistics[ 2 ] = tValues.back();
        }

        /**
         * one interface-preserving redistance of aField on a 24 x 24 or 12 x 12 x 12 mesh of order aOrder: distance targets with the sign
         * of the field, the freeze band, then the band rescale with crossing tolerance aTol (0: rescale not called with a
         * spread, i.e. the old mode) and target slope range aRange
         */
        Result
        run( Field aField, uint aDim, uint aOrder, real aTol, real aRange )
        {
            uint                   tElements = aDim == 2 ? 24 : 12;
            Module_Parameter_Lists tParameters( Module_Type::HMR );
            if ( aDim == 2 )
            {
                tParameters.set( "number_of_elements_per_dimension", tElements, tElements );
                tParameters.set( "domain_dimensions", 2.0, 2.0 );
                tParameters.set( "domain_offset", -1.0, -1.0 );
            }
            else
            {
                tParameters.set( "number_of_elements_per_dimension", tElements, tElements, tElements );
                tParameters.set( "domain_dimensions", 2.0, 2.0, 2.0 );
                tParameters.set( "domain_offset", -1.0, -1.0, -1.0 );
            }
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

            // field values and distance targets at the nodes
            uint             tNumNodes = tIPMesh->get_num_nodes();
            Matrix< DDRMat > tNodal( tNumNodes, 1, 0.0 );
            Matrix< DDRMat > tTargets( tNumNodes, 1, 0.0 );
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                Matrix< DDRMat > tX = tIPMesh->get_mtk_vertex( iNode ).get_coords();
                tNodal( iNode )     = value( aField, tX );
                tTargets( iNode )   = ( tNodal( iNode ) >= 0.0 ? 1.0 : -1.0 ) * std::abs( distance( aField, tX ) );
            }

            // coefficients: for order 1 the coefficient of a node is its value (checked), otherwise the L2 projection
            Matrix< DDRMat > tPsi;
            Matrix< DDRMat > tNew;
            if ( aOrder == 1 )
            {
                uint tNumCoeffs = tIPMesh->get_num_entities( mtk::EntityRank::BSPLINE );
                tPsi.set_size( tNumCoeffs, 1, 0.0 );
                tNew.set_size( tNumCoeffs, 1, 0.0 );
                for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
                {
                    Matrix< IndexMat > tIndices = tIPMesh->get_coefficient_indices_of_node( iNode, 0 );
                    REQUIRE( tIndices.numel() == 1 );
                    REQUIRE( tIPMesh->get_t_matrix_of_node_loc_ind( iNode, 0 )( 0 ) == 1.0 );
                    tPsi( tIndices( 0 ) ) = tNodal( iNode );
                    tNew( tIndices( 0 ) ) = tTargets( iNode );
                }
            }
            else
            {
                tPsi = project( tMeshPair, tNodal );
                tNew = project( tMeshPair, tTargets );
            }

            Vector< Vector< moris_index > > tCellCoefficients = Reinitialize_Performer::cell_coefficients( tIPMesh, 0 );

            // the old mode: freeze band only
            uint tRings = 0;
            Reinitialize_Performer::freeze_interface_band( tNew, tPsi, tCellCoefficients, 1.0e-12, tRings );
            const Matrix< DDRMat > tFrozen = tNew;

            // the new mode
            Result tResult;
            tResult.mReport = Reinitialize_Performer::rescale_interface_band_on_mesh(
                    tIPMesh, 0, tPsi, tNew, 1.0e-12, 4.0 * std::atanh( aTol ), aRange, -1.0e3, 1.0e3 );

            // invariants against psi
            const Cut tCutPsi = cut( tIPMesh, tPsi );
            const Cut tCutNew = cut( tIPMesh, tNew );
            for ( uint iNode = 0; iNode < tNumNodes; iNode++ )
            {
                tResult.mNodeSignFlips += ( tCutPsi.mNodal( iNode ) >= 0.0 ) != ( tCutNew.mNodal( iNode ) >= 0.0 ) ? 1 : 0;
            }
            std::vector< bool > tBandCoefficient( tPsi.numel(), false );
            for ( uint iCell = 0; iCell < tCutPsi.mCutCell.size(); iCell++ )
            {
                tResult.mCutCells += tCutPsi.mCutCell[ iCell ] ? 1 : 0;
                tResult.mCutCellsEqual = tResult.mCutCellsEqual && tCutPsi.mCutCell[ iCell ] == tCutNew.mCutCell[ iCell ];

                // band cells as the rescale defines them
                bool tAllPos = true;
                bool tAllNeg = true;
                for ( uint iC = 0; iC < tCellCoefficients( iCell ).size(); iC++ )
                {
                    tAllPos = tAllPos && tPsi( tCellCoefficients( iCell )( iC ) ) > 1.0e-12;
                    tAllNeg = tAllNeg && tPsi( tCellCoefficients( iCell )( iC ) ) < -1.0e-12;
                }
                if ( !tAllPos && !tAllNeg )
                {
                    for ( uint iC = 0; iC < tCellCoefficients( iCell ).size(); iC++ )
                    {
                        tBandCoefficient[ tCellCoefficients( iCell )( iC ) ] = true;
                    }
                }
            }
            for ( uint iC = 0; iC < tPsi.numel(); iC++ )
            {
                if ( tBandCoefficient[ iC ] )
                {
                    tResult.mCoeffSignFlips += ( tPsi( iC ) > 0.0 ) != ( tNew( iC ) > 0.0 ) ? 1 : 0;
                }
                else
                {
                    tResult.mFarBitIdentical = tResult.mFarBitIdentical && tNew( iC ) == tFrozen( iC );
                }
                tResult.mBandBitIdentical = tResult.mBandBitIdentical && tNew( iC ) == tFrozen( iC );
            }
            tResult.mCutEdges      = tCutPsi.mFraction.size();
            tResult.mCutEdgesEqual = tCutPsi.mFraction.size() == tCutNew.mFraction.size();
            for ( const std::pair< const std::pair< uint, uint >, real >& tEntry : tCutPsi.mFraction )
            {
                std::map< std::pair< uint, uint >, real >::const_iterator tFound = tCutNew.mFraction.find( tEntry.first );
                if ( tFound == tCutNew.mFraction.end() )
                {
                    tResult.mCutEdgesEqual = false;
                    continue;
                }
                tResult.mMaxShift = std::max( tResult.mMaxShift, std::abs( tFound->second - tEntry.second ) );
            }

            slope_statistics( Reinitialize_Performer::cell_slopes( tIPMesh, 0, tPsi ), tCutPsi.mCutCell, tResult.mSlopeBefore );
            slope_statistics( Reinitialize_Performer::cell_slopes( tIPMesh, 0, tNew ), tCutPsi.mCutCell, tResult.mSlopeAfter );

            // the mesh pair owns (and deletes) both meshes
            return tResult;
        }

        void
        print( const char* aName, const Result& aResult )
        {
            std::printf(
                    "[band rescale] %-22s cut cells %4u edges %4u | slope min/median/max %.3f %.3f %.3f -> %.3f %.3f %.3f | "
                    "max crossing shift %.3e | scaled %u of %u, factors %.3f..%.3f, spread %.4f, sweeps %u, backoffs %u\n",
                    aName,
                    aResult.mCutCells,
                    aResult.mCutEdges,
                    aResult.mSlopeBefore[ 0 ],
                    aResult.mSlopeBefore[ 1 ],
                    aResult.mSlopeBefore[ 2 ],
                    aResult.mSlopeAfter[ 0 ],
                    aResult.mSlopeAfter[ 1 ],
                    aResult.mSlopeAfter[ 2 ],
                    aResult.mMaxShift,
                    aResult.mReport.mScaledCoeffs,
                    aResult.mReport.mBandCoeffs,
                    aResult.mReport.mMinFactor,
                    aResult.mReport.mMaxFactor,
                    aResult.mReport.mMaxLogSpread,
                    aResult.mReport.mLimiterPasses,
                    aResult.mReport.mBackoffs );
        }

        /** the invariants every application must keep; aShiftBound in edge lengths */
        void
        check_invariants( const Result& aResult, real aShiftBound )
        {
            REQUIRE( aResult.mCutCells > 20 );
            CHECK( aResult.mNodeSignFlips == 0 );
            CHECK( aResult.mCoeffSignFlips == 0 );
            CHECK( aResult.mCutCellsEqual );
            CHECK( aResult.mCutEdgesEqual );
            CHECK( aResult.mFarBitIdentical );
            CHECK( aResult.mMaxShift <= aShiftBound );
        }
    }    // namespace band_rescale_test

    //------------------------------------------------------------------------------

    TEST_CASE( "Band rescale: a signed distance is a fixed point", "[WRK][sdf_preserve_interface][sdf_preserve_slope]" )
    {
        if ( par_size() == 1 )
        {
            using band_rescale_test::Field;
            for ( Field tField : { Field::PLANE, Field::SPHERE, Field::NODE_ON_SPHERE_SDF } )
            {
                band_rescale_test::Result tResult = band_rescale_test::run( tField, 2, 1, 0.02, 1.25 );
                band_rescale_test::print( "signed distance", tResult );
                band_rescale_test::check_invariants( tResult, 0.0 );

                // slopes already in [0.8, 1.25]: nothing is asked for, the frozen band is returned bit for bit
                CHECK( tResult.mSlopeBefore[ 0 ] > 0.8 );
                CHECK( tResult.mSlopeBefore[ 2 ] < 1.25 );
                CHECK( tResult.mReport.mScaledCoeffs == 0 );
                CHECK( tResult.mBandBitIdentical );
            }
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "Band rescale: a uniformly flat or steep band gets one factor and keeps its crossings", "[WRK][sdf_preserve_interface][sdf_preserve_slope]" )
    {
        if ( par_size() == 1 )
        {
            // 0.2 x distance to a plane: every band coefficient x 4, slope 0.2 -> 0.8 = 1 / range
            band_rescale_test::Result tFlat = band_rescale_test::run( band_rescale_test::Field::PLANE_FLAT, 2, 1, 0.02, 1.25 );
            band_rescale_test::print( "plane x 0.2", tFlat );
            band_rescale_test::check_invariants( tFlat, 1.0e-12 );
            CHECK( tFlat.mSlopeBefore[ 2 ] == Approx( 0.2 ).epsilon( 1.0e-10 ) );
            CHECK( tFlat.mSlopeAfter[ 0 ] == Approx( 0.8 ).epsilon( 1.0e-10 ) );
            CHECK( tFlat.mSlopeAfter[ 2 ] == Approx( 0.8 ).epsilon( 1.0e-10 ) );
            CHECK( tFlat.mReport.mMinFactor == Approx( 4.0 ).epsilon( 1.0e-10 ) );
            CHECK( tFlat.mReport.mMaxFactor == Approx( 4.0 ).epsilon( 1.0e-10 ) );
            CHECK( tFlat.mReport.mScaledCoeffs == tFlat.mReport.mBandCoeffs );

            // 3 x distance to a circle: slopes about 3 -> about 1.25 = range, crossings within the tolerance
            band_rescale_test::Result tSteep = band_rescale_test::run( band_rescale_test::Field::SPHERE_STEEP, 2, 1, 0.02, 1.25 );
            band_rescale_test::print( "circle x 3", tSteep );
            band_rescale_test::check_invariants( tSteep, 0.02 + 1.0e-12 );
            CHECK( tSteep.mSlopeBefore[ 0 ] > 2.9 );
            CHECK( tSteep.mSlopeAfter[ 0 ] > 1.2 );
            CHECK( tSteep.mSlopeAfter[ 2 ] < 1.3 );
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "Band rescale: flattened fields move towards unit slope within the crossing tolerance", "[WRK][sdf_preserve_interface][sdf_preserve_slope]" )
    {
        if ( par_size() == 1 )
        {
            using band_rescale_test::Field;
            const real tTol = 0.02;

            // 0.25 w tanh( d / w ): flat everywhere near the interface
            band_rescale_test::Result tTanh = band_rescale_test::run( Field::SPHERE_TANH, 2, 1, tTol, 1.25 );
            band_rescale_test::print( "circle tanh", tTanh );
            band_rescale_test::check_invariants( tTanh, tTol + 1.0e-12 );
            CHECK( tTanh.mSlopeBefore[ 2 ] < 0.26 );
            CHECK( tTanh.mSlopeAfter[ 0 ] > 0.7 );
            CHECK( tTanh.mSlopeAfter[ 2 ] < 1.25 );

            // d^3 / h^2: zero slope on the interface; the cubed nodal values put the crossings elsewhere than the circle and
            // the rescale keeps them there
            band_rescale_test::Result tCubed = band_rescale_test::run( Field::SPHERE_CUBED, 2, 1, tTol, 1.25 );
            band_rescale_test::print( "circle cubed", tCubed );
            band_rescale_test::check_invariants( tCubed, tTol + 1.0e-12 );
            CHECK( tCubed.mSlopeAfter[ 1 ] > tCubed.mSlopeBefore[ 1 ] );

            // slope 0.3 to 1.3 along the interface: the limiter is active, the correction is partial
            band_rescale_test::Result tVarying = band_rescale_test::run( Field::SPHERE_VARYING, 2, 1, tTol, 1.25 );
            band_rescale_test::print( "circle varying", tVarying );
            band_rescale_test::check_invariants( tVarying, tTol + 1.0e-12 );
            CHECK( tVarying.mReport.mLimiterPasses > 0 );
            CHECK( tVarying.mReport.mMaxLogSpread <= 4.0 * std::atanh( tTol ) + 1.0e-11 );
            CHECK( tVarying.mSlopeAfter[ 0 ] > tVarying.mSlopeBefore[ 0 ] );

            // the same with a node 1e-6 h off the interface: it keeps its sign and its crossings
            band_rescale_test::Result tNode = band_rescale_test::run( Field::NODE_ON_SPHERE, 2, 1, tTol, 1.25 );
            band_rescale_test::print( "node 1e-6 h off", tNode );
            band_rescale_test::check_invariants( tNode, tTol + 1.0e-12 );
            CHECK( tNode.mSlopeAfter[ 0 ] > tNode.mSlopeBefore[ 0 ] );

            // a larger tolerance buys more slope
            band_rescale_test::Result tLoose = band_rescale_test::run( Field::SPHERE_VARYING, 2, 1, 0.1, 1.25 );
            band_rescale_test::print( "circle varying tol 0.1", tLoose );
            band_rescale_test::check_invariants( tLoose, 0.1 + 1.0e-12 );
            CHECK( tLoose.mSlopeAfter[ 0 ] > tVarying.mSlopeAfter[ 0 ] );

            // 3D, HEX8
            band_rescale_test::Result tHex = band_rescale_test::run( Field::SPHERE_VARYING, 3, 1, tTol, 1.25 );
            band_rescale_test::print( "sphere varying 3D", tHex );
            band_rescale_test::check_invariants( tHex, tTol + 1.0e-12 );
            CHECK( tHex.mSlopeAfter[ 0 ] > tHex.mSlopeBefore[ 0 ] );
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "Band rescale: quadratic B-splines keep every nodal sign", "[WRK][sdf_preserve_interface][sdf_preserve_slope]" )
    {
        if ( par_size() == 1 )
        {
            // coefficients are not nodal values: the coefficient signs and the cut cells are kept, the nodal signs are
            // checked and backed off if needed; the crossing bound of order 1 does not apply and is only reported
            band_rescale_test::Result tResult = band_rescale_test::run( band_rescale_test::Field::SPHERE_VARYING, 2, 2, 0.02, 1.25 );
            band_rescale_test::print( "circle varying, p = 2", tResult );
            REQUIRE( tResult.mCutCells > 20 );
            CHECK( tResult.mNodeSignFlips == 0 );
            CHECK( tResult.mCoeffSignFlips == 0 );
            CHECK( tResult.mCutCellsEqual );
            CHECK( tResult.mCutEdgesEqual );
            CHECK( tResult.mFarBitIdentical );
            CHECK( tResult.mReport.mScaledCoeffs > 0 );
            CHECK( tResult.mMaxShift < 0.1 );
        }
    }

    //------------------------------------------------------------------------------

    TEST_CASE( "Band rescale: off by default, and zero tolerance is the frozen band bit for bit", "[WRK][sdf_preserve_interface][sdf_preserve_slope]" )
    {
        if ( par_size() == 1 )
        {
            // parameter defaults: both the freeze band and its rescale are off
            Parameter_List tProblem = prm::create_opt_problem_parameter_list();
            CHECK( tProblem.get< bool >( "sdf_preserve_interface" ) == false );
            CHECK( tProblem.get< real >( "sdf_preserve_interface_slope_tol" ) == 0.0 );
            CHECK( tProblem.get< real >( "sdf_preserve_interface_slope_range" ) == 1.25 );

            // zero tolerance = zero spread: unequal requests cannot be met, the band stays frozen
            band_rescale_test::Result tZero = band_rescale_test::run( band_rescale_test::Field::SPHERE_VARYING, 2, 1, 0.0, 1.25 );
            band_rescale_test::print( "circle varying tol 0", tZero );
            band_rescale_test::check_invariants( tZero, 0.0 );
            CHECK( tZero.mReport.mScaledCoeffs == 0 );
            CHECK( tZero.mBandBitIdentical );
        }
    }
}    // namespace moris::wrk
