/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_FEM_Element_Nonconformal_Sideset.cpp
 *
 */

#include "cl_FEM_Element_Nonconformal_Sideset.hpp"
#include "cl_FEM_Geometry_Interpolator.hpp"
#include "cl_FEM_Element_Double_Sideset.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_Remap_Nonconformal_Ray.hpp"
#include <iostream>
#include "cl_MTK_Enums.hpp"
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Vector.hpp"
#include "linalg_typedefs.hpp"

#include <cl_FEM_Model.hpp>

namespace moris::fem
{
    Matrix< DDRMat > Element_Nonconformal_Sideset::get_leader_integration_point( uint const aGPIndex ) const
    {
        return mIntegrationPointPairs.get_leader_coordinates().get_column( aGPIndex );
    }

    //----------------------------------------------------------------

    Matrix< DDRMat > Element_Nonconformal_Sideset::get_follower_integration_point( uint const aGPIndex ) const
    {
        return mIntegrationPointPairs.get_follower_coordinates().get_column( aGPIndex );
    }

    //----------------------------------------------------------------

    real Element_Nonconformal_Sideset::get_integration_weight( uint aGPIndex ) const
    {
        return mIntegrationPointPairs.get_integration_weights()( aGPIndex );
    }

    //----------------------------------------------------------------

    uint Element_Nonconformal_Sideset::get_number_of_integration_points() const
    {
        return mIntegrationPointPairs.get_integration_weights().size();
    }

    //----------------------------------------------------------------

    void Element_Nonconformal_Sideset::compute_jacobian_and_residual()
    {
        Element_Double_Sideset::compute_jacobian_and_residual();

        // undo the increment from the double sided sideset
        mSet->mFemModel->mDoubleSidedSideSetsGaussPoints -= get_number_of_integration_points();

        mSet->mFemModel->mNonconformalSideSetsGaussPoints += get_number_of_integration_points();
    }

    //----------------------------------------------------------------

    moris_index Element_Nonconformal_Sideset::get_follower_local_cell_index() const
    {
        return mFollowerCellIndexInCluster;
    }

    //----------------------------------------------------------------

    void Element_Nonconformal_Sideset::remap_follower_for_output() const
    {
        // OUTPUT-ONLY current-configuration follower remap.
        //
        // The VIS clusters carry integration/nodal point pairs that were copied from the FEM
        // clusters when the VIS mesh was built (near the reference configuration) and are never
        // refreshed by the solver's per-iteration remap. Evaluating a nonconformal quantity such as
        // the contact gap against that frozen follower pairing collapses the gap to ~0 everywhere.
        //
        // Here we re-fire the contact ray in the current (deformed) configuration - from the current
        // leader point (already set on the leader interpolators) along the current leader normal onto
        // the current follower facet - and reposition the follower field interpolator at the hit.
        // This mirrors the remap the IWG applies during its residual/Jacobian finite-difference sweep.
        //
        // TODO: 2D line facets only (see fn_FEM_Remap_Nonconformal_Ray); 3D needs a ray-triangle map.

        // the remap requires a displacement field on both sides; if none is present, do nothing
        Field_Interpolator* tLeaderDisplacementFI =
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )
                        ->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        Field_Interpolator* tFollowerDisplacementFI =
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
                        ->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        if ( tLeaderDisplacementFI == nullptr || tFollowerDisplacementFI == nullptr )
        {
            return;
        }

        Geometry_Interpolator* tLeaderIGGI =
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::LEADER )->get_IG_geometry_interpolator();
        Geometry_Interpolator* tFollowerIGGI =
                mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )->get_IG_geometry_interpolator();

        Matrix< DDRMat > const tRemappedFollowerCoords =
                remap_nonconformal_ray( tLeaderIGGI, tFollowerIGGI, tLeaderDisplacementFI, tFollowerDisplacementFI );

        mSet->get_field_interpolator_manager( mtk::Leader_Follower::FOLLOWER )
                ->set_space_time_from_local_IG_point( tRemappedFollowerCoords );
    }
}    // namespace moris::fem
