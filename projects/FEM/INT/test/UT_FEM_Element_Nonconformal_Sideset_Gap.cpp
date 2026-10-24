/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Element_Nonconformal_Sideset_Gap.cpp
 *
 * Guards the OUTPUT-ONLY current-configuration follower remap used by nonconformal contact IQIs
 * (e.g. IQI_Gap). Regression target: without the remap the follower field interpolator is evaluated
 * at the frozen reference-configuration VIS pairing, collapsing the contact gap to ~0 everywhere.
 *
 * The test builds a 2D line-facet leader/follower pair with a prescribed, non-uniform deformation
 * (the follower slides tangentially and opens by a spatially varying amount), then exercises the
 * shared current-configuration ray remap fem::remap_nonconformal_ray() - the exact routine the
 * output path (Element_Nonconformal_Sideset::remap_follower_for_output) calls - and asserts that the
 * resulting gap g = (x_follower_current - x_leader_current) . n_leader_current equals the analytic
 * current-configuration ray/segment separation: nonzero, spatially varying, and DIFFERENT from the
 * (buggy) frozen-reference-pairing value.
 */

#include "catch.hpp"

#include "cl_MTK_Interpolation_Rule.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "fn_FEM_Remap_Nonconformal_Ray.hpp"

#include "cl_MSI_Dof_Type_Enums.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "fn_dot.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"

using namespace moris;
using namespace fem;

namespace
{
    // A minimal faithful "side" of a nonconformal contact interface: a 2D bulk QUAD carrying a
    // displacement field interpolator (space dim 2, as the real IG displacement FI has) plus a LINE
    // IG geometry interpolator restricted to one edge (the contact facet), mirroring how the FEM
    // set wires an IG side geometry interpolator to an IP bulk displacement FI.
    struct Contact_Side
    {
        Geometry_Interpolator* mBulkGI = nullptr;    // IP (bulk) geometry interpolator, 2D
        Geometry_Interpolator* mSideGI = nullptr;    // IG (side) geometry interpolator, LINE facet
        Field_Interpolator*    mDispFI = nullptr;    // displacement FI, built on the bulk GI

        ~Contact_Side()
        {
            delete mDispFI;
            delete mSideGI;
            delete mBulkGI;
        }
    };

    // Build one side of the interface.
    //  aBulkPhysCoords    : 4x2 physical coordinates of the bulk QUAD nodes
    //  aBulkDispCoeffs    : 4x2 nodal displacement coefficients (ux, uy) of the bulk QUAD
    //  aFacetPhysCoords   : 2x2 physical coordinates of the two facet (edge) nodes
    //  aFacetParentParams : 2x2 parametric coordinates of those facet nodes within the bulk QUAD
    void build_side(
            Contact_Side&           aSide,
            const Matrix< DDRMat >& aBulkPhysCoords,
            const Matrix< DDRMat >& aBulkDispCoeffs,
            const Matrix< DDRMat >& aFacetPhysCoords,
            const Matrix< DDRMat >& aFacetParentParams )
    {
        Matrix< DDRMat > tTHat  = { { 0.0 }, { 1.0 } };
        Matrix< DDRMat > tTauId = { { -1.0 }, { 1.0 } };

        // reference parametric corners of a LINEAR QUAD
        Matrix< DDRMat > tQuadParam = {
            { -1.0, -1.0 },
            { +1.0, -1.0 },
            { +1.0, +1.0 },
            { -1.0, +1.0 }
        };

        // IP (bulk) geometry interpolator - QUAD, on which the displacement FI lives (space dim 2)
        mtk::Interpolation_Rule tBulkRule(
                mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Geometry_Type::LINE,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        aSide.mBulkGI = new Geometry_Interpolator( tBulkRule );
        aSide.mBulkGI->set_coeff( aBulkPhysCoords, tTHat );
        aSide.mBulkGI->set_space_param_coeff( tQuadParam );

        // displacement field interpolator (2 fields), built on the bulk GI - get_space_dim() == 2
        mtk::Interpolation_Rule tFIRule(
                mtk::Geometry_Type::QUAD,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        Vector< MSI::Dof_Type > tDispDof = { MSI::Dof_Type::UX };
        aSide.mDispFI                    = new Field_Interpolator( 2, tFIRule, aSide.mBulkGI, tDispDof );

        // the FI carries SPACE-TIME coefficients (4 spatial bases x 2 linear-time bases = 8 rows);
        // the field is prescribed time-constant, so both time levels get the same spatial values.
        uint const       tNumSpaceBases = aBulkDispCoeffs.n_rows();
        uint const       tNumFields     = aBulkDispCoeffs.n_cols();
        Matrix< DDRMat > tSpaceTimeDisp( 2 * tNumSpaceBases, tNumFields, 0.0 );
        for ( uint iNode = 0; iNode < tNumSpaceBases; ++iNode )
        {
            for ( uint iField = 0; iField < tNumFields; ++iField )
            {
                tSpaceTimeDisp( iNode, iField )                  = aBulkDispCoeffs( iNode, iField );
                tSpaceTimeDisp( iNode + tNumSpaceBases, iField ) = aBulkDispCoeffs( iNode, iField );
            }
        }
        aSide.mDispFI->set_coeff( tSpaceTimeDisp );

        // initialize the FI evaluation point (bulk param + time) so save/restore during the current-
        // configuration node evaluation (get_space_coeff_current) operates on a valid point
        Matrix< DDRMat > tInitFIPoint = { { 0.0 }, { 0.0 }, { -1.0 } };
        aSide.mDispFI->set_space_time( tInitFIPoint );

        // IG (side) geometry interpolator - LINE facet embedded in 2D, exactly as a sideset IG GI:
        // its space_param_coeff are the facet nodes' parametric coordinates within the parent QUAD so
        // that the bulk displacement FI can be evaluated at the facet nodes (get_space_coeff_current).
        mtk::Interpolation_Rule tSideRule(
                mtk::Geometry_Type::LINE,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR,
                mtk::Geometry_Type::LINE,
                mtk::Interpolation_Type::LAGRANGE,
                mtk::Interpolation_Order::LINEAR );

        aSide.mSideGI = new Geometry_Interpolator( tSideRule, tBulkRule, mtk::CellShape::GENERAL, true, false );
        aSide.mSideGI->set_space_coeff( aFacetPhysCoords );
        aSide.mSideGI->set_time_coeff( tTHat );
        aSide.mSideGI->set_space_param_coeff( aFacetParentParams );
        aSide.mSideGI->set_time_param_coeff( tTauId );
    }

    // analytic ray / line-segment intersection (2D), mirroring mtk::Ray_Line_Intersection.
    // returns the physical intersection point; aHasHit reports whether v in [0,1].
    Matrix< DDRMat > ray_segment_intersection(
            const Matrix< DDRMat >& aOrigin,       // 2x1
            const Matrix< DDRMat >& aDirection,    // 2x1 (unit)
            const Matrix< DDRMat >& aSegStart,     // 2x1
            const Matrix< DDRMat >& aSegEnd,       // 2x1
            bool&                   aHasHit )
    {
        Matrix< DDRMat > tSpan     = aSegEnd - aSegStart;
        Matrix< DDRMat > tdOrigins = aSegStart - aOrigin;

        real tDetD = tSpan( 0 ) * aDirection( 1 ) - tSpan( 1 ) * aDirection( 0 );
        real tV    = ( tdOrigins( 1 ) * aDirection( 0 ) - tdOrigins( 0 ) * aDirection( 1 ) ) / tDetD;

        aHasHit = ( tV >= 0.0 && tV <= 1.0 );
        return aSegStart + tV * tSpan;
    }
}    // namespace

TEST_CASE( "FEM Nonconformal Sideset output gap remap",
        "[moris],[fem],[FEM_NCSS_Gap]" )
{
    real const tEps = 1.0e-9;

    // ---------------------------------------------------------------------------------------------
    // Leader: bulk QUAD BELOW the interface; contact facet is the reference segment (2,0)-(0,0).
    // Node order of the facet (A=(2,0), B=(0,0)) is chosen so the outward LINE normal points +y,
    // i.e. toward the follower. Leader is undeformed (its current config == reference).
    // ---------------------------------------------------------------------------------------------
    Contact_Side tLeader;
    {
        Matrix< DDRMat > tBulkPhys = {
            { 0.0, 0.0 },
            { 2.0, 0.0 },
            { 2.0, -1.0 },
            { 0.0, -1.0 }
        };
        Matrix< DDRMat > tBulkDisp( 4, 2, 0.0 );    // leader fixed
        Matrix< DDRMat > tFacetPhys = {
            { 2.0, 0.0 },
            { 0.0, 0.0 }
        };
        Matrix< DDRMat > tFacetParam = {
            { 1.0, -1.0 },
            { -1.0, -1.0 }
        };
        build_side( tLeader, tBulkPhys, tBulkDisp, tFacetPhys, tFacetParam );
    }

    // ---------------------------------------------------------------------------------------------
    // Follower: bulk QUAD ABOVE the interface; contact facet is the reference segment (0,1)-(2,1).
    // Prescribed deformation slides the follower tangentially and opens it non-uniformly:
    //   node (0,1) -> (0.30, 1.20),   node (2,1) -> (1.90, 1.50).
    // ---------------------------------------------------------------------------------------------
    Contact_Side tFollower;
    {
        Matrix< DDRMat > tBulkPhys = {
            { 0.0, 1.0 },
            { 2.0, 1.0 },
            { 2.0, 2.0 },
            { 0.0, 2.0 }
        };
        Matrix< DDRMat > tBulkDisp = {
            { 0.30, 0.20 },
            { -0.10, 0.50 },
            { 0.0, 0.0 },
            { 0.0, 0.0 }
        };
        Matrix< DDRMat > tFacetPhys = {
            { 0.0, 1.0 },
            { 2.0, 1.0 }
        };
        Matrix< DDRMat > tFacetParam = {
            { -1.0, -1.0 },
            { 1.0, -1.0 }
        };
        build_side( tFollower, tBulkPhys, tBulkDisp, tFacetPhys, tFacetParam );
    }

    // follower current facet nodes (for the analytic reference)
    Matrix< DDRMat > tF0c = { { 0.30 }, { 1.20 } };
    Matrix< DDRMat > tF1c = { { 1.90 }, { 1.50 } };

    // leader integration/eval points on the facet (LINE parametric)
    Vector< real > tEtaPoints = { -0.6, 0.0, 0.6 };

    Vector< real > tGaps;
    Vector< real > tGapsFrozenReference;

    for ( real tEta : tEtaPoints )
    {
        // place the leader interpolators at the current leader point (side param -> mapped bulk param)
        Matrix< DDRMat > tLeaderParam = { { tEta }, { 0.0 } };
        tLeader.mSideGI->set_space_time( tLeaderParam );
        tLeader.mSideGI->map_integration_point();

        // seed the follower parametric point (as the output path does before remapping)
        Matrix< DDRMat > tFollowerSeed = { { 0.0 }, { 0.0 } };
        tFollower.mSideGI->set_space_time( tFollowerSeed );
        tFollower.mSideGI->map_integration_point();

        // current leader point and normal (from the code under test)
        Matrix< DDRMat > tLeaderCurrent = tLeader.mSideGI->valx_current( tLeader.mDispFI );
        Matrix< DDRMat > tNormalCurrent = tLeader.mSideGI->get_normal_current( tLeader.mDispFI );

        // the normal must point toward the follower (+y) for the ray to hit
        REQUIRE( tNormalCurrent( 1 ) > 0.0 );

        // ---- current-configuration remap (the fix) --------------------------------------------
        Matrix< DDRMat > tRemapped = remap_nonconformal_ray(
                tLeader.mSideGI, tFollower.mSideGI, tLeader.mDispFI, tFollower.mDispFI );

        tFollower.mSideGI->set_space_time( tRemapped );
        tFollower.mSideGI->map_integration_point();
        Matrix< DDRMat > tFollowerCurrent = tFollower.mSideGI->valx_current( tFollower.mDispFI );

        real tGap = dot( tFollowerCurrent - tLeaderCurrent, tNormalCurrent );
        tGaps.push_back( tGap );

        // ---- analytic current-config ray/segment intersection ---------------------------------
        Matrix< DDRMat > tOrigin = trans( tLeaderCurrent );
        bool             tHit    = false;
        Matrix< DDRMat > tPStar  = ray_segment_intersection( tOrigin, tNormalCurrent, tF0c, tF1c, tHit );

        REQUIRE( tHit );

        real tGapAnalytic = dot( tPStar - tOrigin, tNormalCurrent );

        // the code's remapped follower must land on the analytic current-config intersection
        CHECK( std::abs( tFollowerCurrent( 0 ) - tPStar( 0 ) ) < tEps );
        CHECK( std::abs( tFollowerCurrent( 1 ) - tPStar( 1 ) ) < tEps );
        CHECK( std::abs( tGap - tGapAnalytic ) < tEps );

        // ---- control: gap WITHOUT remap (frozen reference-config pairing) ----------------------
        // reference ray from the (undeformed) leader point hits the follower reference facet at the
        // same x; that frozen parametric location is xi_ref = x0 - 1 on this facet.
        real             tXiFrozen    = tLeaderCurrent( 0 ) - 1.0;
        Matrix< DDRMat > tFrozenParam = { { tXiFrozen }, { 0.0 } };
        tFollower.mSideGI->set_space_time( tFrozenParam );
        tFollower.mSideGI->map_integration_point();
        Matrix< DDRMat > tFollowerFrozen = tFollower.mSideGI->valx_current( tFollower.mDispFI );
        real             tGapFrozen      = dot( tFollowerFrozen - tLeaderCurrent, tNormalCurrent );
        tGapsFrozenReference.push_back( tGapFrozen );
    }

    // gaps must be strictly positive (open, separated interface) - the bug produced ~0
    for ( real tGap : tGaps )
    {
        CHECK( tGap > 0.5 );
    }

    // gaps must be spatially varying (not a constant field)
    real tMin = tGaps( 0 );
    real tMax = tGaps( 0 );
    for ( real tGap : tGaps )
    {
        tMin = std::min( tMin, tGap );
        tMax = std::max( tMax, tGap );
    }
    CHECK( ( tMax - tMin ) > 0.1 );

    // the remap must matter: at least one point differs measurably from the frozen-reference value
    real tMaxRemapEffect = 0.0;
    for ( uint i = 0; i < tGaps.size(); ++i )
    {
        tMaxRemapEffect = std::max( tMaxRemapEffect, std::abs( tGaps( i ) - tGapsFrozenReference( i ) ) );
    }
    CHECK( tMaxRemapEffect > 1.0e-3 );
}
