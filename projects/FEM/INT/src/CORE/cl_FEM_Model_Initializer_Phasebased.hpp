/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 * ------------------------------------------------------------------------------------
 *
 * cl_FEM_Model_Initializer_Phasebased.hpp
 *
 */

#ifndef MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP
#define MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP

#include "cl_FEM_Model_Initializer.hpp"
#include "cl_FEM_Phase_User_Info.hpp"
#include <map>
#include <utility>
#include "cl_FEM_Set_User_Info.hpp"
namespace moris::fem
{
    class Model_Initializer_Phasebased : public Model_Initializer
    {
      public:
        Model_Initializer_Phasebased(
                const Module_Parameter_Lists        &aParameterList,
                std::shared_ptr< Library_IO >                    aLibrary,
                mtk::Mesh_Pair const                            *aMeshPair,
                uint                                             aSpatialDimension,
                bool                                             aUseNewGhostSets,
                std::unordered_map< MSI::Dof_Type, moris_index > aDofTypeToBsplineMeshIndex )
                : Model_Initializer(
                          aParameterList,
                          aMeshPair,
                          std::move( aLibrary ),
                          aSpatialDimension,
                          aUseNewGhostSets,
                          std::move( aDofTypeToBsplineMeshIndex ) ) {};

      protected:

      public:
        void initialize() override;

      protected:
        void create_material_models() override;
        void create_constitutive_models() override;
        void create_follower_side_models() override;
        void create_stabilization_parameters() override;
        void create_iwgs() override;
        void create_iqis() override;
        void create_set_info() override;
        void print_physics_model() override;

      private:
        void create_phases();

        void get_mesh_set_names(
                fem::Element_Type      aBulkType,
                const std::string     &aLeaderPhaseName,
                const std::string     &aFollowerPhaseName,
                const std::string     &aNeighborPhaseString,
                const std::string     &aSideOrdinalsString,
                bool                   aIsGhost,
                Vector< std::string > &aMeshSetNames );

        // property of the given index, and the phase info holding the MMs and CMs, for the leader or the follower side
        std::shared_ptr< Property > const &get_property( uint aIndex, mtk::Leader_Follower aLeaderFollowerType ) const;
        fem::Phase_User_Info               &get_phase_info( uint aPhaseIndex, mtk::Leader_Follower aLeaderFollowerType );

        Vector< fem::Phase_User_Info > mPhaseInfo;
        std::map< std::string, uint >  mPhaseMap;

        // follower-side instances (MB-047): properties, and per phase the MMs and CMs, built from the same parameter
        // lists as mProperties and mPhaseInfo. The DOF and DV type lists of a phase stay in mPhaseInfo.
        PointerCell< fem::Property >   mFollowerProperties;
        Vector< fem::Phase_User_Info > mFollowerPhaseInfo;
    };
}    // namespace moris::fem

#endif    // MORIS_CL_FEM_MODEL_INITIALIZER_PHASEBASED_HPP
