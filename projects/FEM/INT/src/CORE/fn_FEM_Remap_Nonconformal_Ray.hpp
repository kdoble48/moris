/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Remap_Nonconformal_Ray.hpp
 *
 */

#ifndef SRC_FEM_FN_FEM_REMAP_NONCONFORMAL_RAY_HPP_
#define SRC_FEM_FN_FEM_REMAP_NONCONFORMAL_RAY_HPP_

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    class Geometry_Interpolator;
    class Field_Interpolator;

    //------------------------------------------------------------------------------
    // Re-fires the nonconformal-contact ray in the CURRENT (deformed) configuration.
    //
    // Given the leader IG geometry interpolator - already positioned at the current leader
    // integration/nodal point - the ray is cast from the current leader point x_leader_current
    // along the current leader normal n_leader_current onto the current follower facet. The
    // returned value is the follower parametric space-time coordinate of the intersection, ready
    // to be handed to Field_Interpolator_Manager::set_space_time_from_local_IG_point() on the
    // follower side.
    //
    // This is the same current-configuration remap that the IWG applies during its
    // residual/Jacobian finite-difference sweep (cl_FEM_IWG::remap_nonconformal_rays delegates
    // here). It is required whenever a nonconformal quantity (e.g. the contact gap) must be
    // evaluated against the true deformed follower surface rather than a frozen reference-config
    // pairing.
    //
    // TODO: currently limited to 2D line facets (mtk::Ray_Line_Intersection). The 3D case needs a
    // ray-triangle generalization.
    //
    // @param[in] aLeaderIGGeometryInterpolator    leader IG geometry interpolator (positioned at
    //                                             the current leader point)
    // @param[in] aFollowerIGGeometryInterpolator  follower IG geometry interpolator
    // @param[in] aLeaderDisplacementFI            leader displacement field interpolator
    // @param[in] aFollowerDisplacementFI          follower displacement field interpolator
    // @return    follower parametric space-time coordinate of the ray hit (unchanged if no hit)
    Matrix< DDRMat >
    remap_nonconformal_ray(
            Geometry_Interpolator* aLeaderIGGeometryInterpolator,
            Geometry_Interpolator* aFollowerIGGeometryInterpolator,
            Field_Interpolator*    aLeaderDisplacementFI,
            Field_Interpolator*    aFollowerDisplacementFI );

    //------------------------------------------------------------------------------

}    // namespace moris::fem

#endif /* SRC_FEM_FN_FEM_REMAP_NONCONFORMAL_RAY_HPP_ */
