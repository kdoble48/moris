/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Remap_Nonconformal_Ray.cpp
 *
 */

#include "fn_FEM_Remap_Nonconformal_Ray.hpp"

#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "fn_trans.hpp"

#include <cl_MTK_Ray_Line_Intersection.hpp>

namespace moris::fem
{
    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    remap_nonconformal_ray(
            Geometry_Interpolator* aLeaderIGGeometryInterpolator,
            Geometry_Interpolator* aFollowerIGGeometryInterpolator,
            Field_Interpolator*    aLeaderDisplacementFI,
            Field_Interpolator*    aFollowerDisplacementFI )
    {
        uint const tDim = aLeaderDisplacementFI->get_space_dim();

        // Since the leader field interpolator might carry a different (perturbed / current)
        // displacement field, force the geometry interpolator to use the updated coordinates for
        // all computations (i.e. not the memoized reference ones).
        aLeaderIGGeometryInterpolator->reset_eval_flags_coordinates();

        // deformed coordinates of the leader element in the current configuration
        Matrix< DDRMat > const tLeaderCoordsCurrent = aLeaderIGGeometryInterpolator->get_space_coeff_current( aLeaderDisplacementFI );

        // current location of the (already-set) leader integration/nodal point
        Matrix< DDRMat > const tIGPointCurrent = aLeaderIGGeometryInterpolator->valx_current( aLeaderDisplacementFI );

        // current leader normal (direction of the contact ray)
        Matrix< DDRMat > const tNormalCurrent = aLeaderIGGeometryInterpolator->get_normal_current( aLeaderDisplacementFI );

        MORIS_ASSERT( tLeaderCoordsCurrent.n_rows() == 2 && tLeaderCoordsCurrent.n_cols() == 2,
                "fem::remap_nonconformal_ray - Nonconformal ray remap is only implemented for line elements." );

        // reset internal flags to make sure the geometry interpolator re-evaluates the current
        // follower element coordinates
        aFollowerIGGeometryInterpolator->reset_eval_flags_coordinates();
        Matrix< DDRMat > const tFollowerCoordinates = aFollowerIGGeometryInterpolator->get_space_coeff_current( aFollowerDisplacementFI );

        // cast the ray from the current leader point along the current leader normal onto the
        // current follower facet
        mtk::Ray_Line_Intersection tRLI( tDim );
        tRLI.set_ray_origin( trans( tIGPointCurrent ) );
        tRLI.set_ray_direction( tNormalCurrent );
        tRLI.set_target_origin( trans( tFollowerCoordinates.get_row( 0 ) ) );
        tRLI.set_target_span( trans( tFollowerCoordinates.get_row( 1 ) - tFollowerCoordinates.get_row( 0 ) ) );
        tRLI.perform_raytracing();

        // start from the follower's current parametric space-time; only the spatial parametric
        // coordinate is updated on a successful hit
        Matrix< DDRMat > tFollowerSpaceTime;
        aFollowerIGGeometryInterpolator->get_space_time( tFollowerSpaceTime );

        if ( tRLI.has_intersection() )
        {
            tFollowerSpaceTime( 0 ) = tRLI.get_intersection_parametric()( 0 );
        }

        return tFollowerSpaceTime;
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
