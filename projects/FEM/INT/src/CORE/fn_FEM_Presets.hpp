/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Presets.hpp
 *
 * Physics presets for input decks (Deck API v2, see doc/internal/DECK_API_RFC.md):
 * one call emits the standard property -> constitutive model -> stabilization ->
 * IWG -> IQI block that is otherwise copied (~150 lines) between decks of the same
 * physics family. Presets emit through the ordinary Module_Parameter_Lists API, so
 * they are usable from BOTH deck styles (a legacy FEMParameterList function or a
 * MORISInputDeck entry point), and everything they create can be tweaked afterwards
 * through normal parameter-list access — a preset is a macro, not a black box. The
 * returned name bundle lists every created entity for follow-up edits or IQI wiring.
 *
 * Header-only on purpose: deck .so files resolve symbols against the moris binary at
 * dlopen, and a preset translation unit referenced only by decks would be dropped
 * from the statically-linked executable.
 */

#pragma once

#include <optional>
#include <string>

#include "cl_Module_Parameter_Lists.hpp"
#include "parameters.hpp"
#include "cl_FEM_Enums.hpp"

namespace moris::fem::presets
{
    //------------------------------------------------------------------------------------------------------------------

    /**
     * Configuration for the linear-elastic preset. Mesh-set strings use the usual
     * XTK/STK grammar; empty strings disable the corresponding piece.
     */
    struct Linear_Elastic_Config
    {
        // mesh sets
        std::string mBulkSets;              // required: bulk sets the material occupies
        std::string mDirichletSets;         // "" = no Dirichlet BC
        std::string mNeumannSets;           // "" = no traction BC
        std::string mGhostSets;             // "" = no ghost stabilization

        // material and BC values (function_parameters strings)
        std::string mDofs           = "UX,UY";
        std::string mYoungs         = "1.0";
        std::string mPoisson        = "0.0";
        std::string mDensity        = "1.0";
        std::string mBedding        = "";              // "" = no bedding property/term
        std::string mDirichletValue = "0.0;0.0";
        std::string mTraction       = "0.0;1.0";
        std::string mTractionFunction = "";            // named value_function for the traction
                                                       // property ("" = constant mTraction)

        // stabilization
        real mNitschePenalty = 100.0;
        real mGhostPenalty   = 0.01;

        // standard IQIs
        bool mAddStrainEnergyIQI = true;
        bool mAddVolumeIQI       = true;

        // name prefix, so several preset instances (e.g. per phase) can coexist
        std::string mPrefix = "";
    };

    /**
     * Names of every entity created by linear_elastic(), for follow-up tweaks
     * (via normal parameter-list access) and for wiring IQIs into GEN/VIS.
     */
    struct Linear_Elastic_Names
    {
        std::string mDensityProp;
        std::string mYoungsProp;
        std::string mPoissonProp;
        std::string mBeddingProp;      // "" if not created
        std::string mDirichletProp;    // "" if not created
        std::string mTractionProp;     // "" if not created
        std::string mCM;
        std::string mNitscheSP;        // "" if not created
        std::string mGhostSP;          // "" if not created
        std::string mBulkIWG;
        std::string mDirichletIWG;     // "" if not created
        std::string mNeumannIWG;       // "" if not created
        std::string mGhostIWG;         // "" if not created
        std::string mStrainEnergyIQI;    // "" if not created
        std::string mVolumeIQI;          // "" if not created
    };

    //------------------------------------------------------------------------------------------------------------------

    /**
     * Emits the standard isotropic linear-elastic block: constant properties (using
     * the built-in constant value function — no Func_Const needed in the deck), the
     * STRUC_LIN_ISO constitutive model, Nitsche Dirichlet + traction Neumann BCs,
     * optional ghost stabilization, and the standard strain-energy / volume IQIs.
     *
     * @param aFemParameterLists The FEM module parameter lists to append to
     * @param aConfig Preset configuration
     * @return Names of all created entities
     */
    inline Linear_Elastic_Names
    linear_elastic(
            Module_Parameter_Lists&       aFemParameterLists,
            const Linear_Elastic_Config& aConfig )
    {
        MORIS_ERROR( !aConfig.mBulkSets.empty(),
                "fem::presets::linear_elastic - mBulkSets is required." );

        Linear_Elastic_Names tNames;
        const std::string&   tP = aConfig.mPrefix;

        // ---- properties (constant: the built-in value function applies) ----------
        auto tAddConstantProperty = [ &aFemParameterLists ]( const std::string& aName, const std::string& aValue ) {
            aFemParameterLists( FEM::PROPERTIES ).add_parameter_list();
            aFemParameterLists.set( "property_name", aName );
            aFemParameterLists.set( "function_parameters", aValue );
        };

        tNames.mDensityProp = tP + "PropDensity";
        tAddConstantProperty( tNames.mDensityProp, aConfig.mDensity );

        tNames.mYoungsProp = tP + "PropYoungs";
        tAddConstantProperty( tNames.mYoungsProp, aConfig.mYoungs );

        tNames.mPoissonProp = tP + "PropPoisson";
        tAddConstantProperty( tNames.mPoissonProp, aConfig.mPoisson );

        if ( !aConfig.mBedding.empty() )
        {
            tNames.mBeddingProp = tP + "PropBedding";
            tAddConstantProperty( tNames.mBeddingProp, aConfig.mBedding );
        }

        if ( !aConfig.mDirichletSets.empty() )
        {
            tNames.mDirichletProp = tP + "PropDirichletU";
            tAddConstantProperty( tNames.mDirichletProp, aConfig.mDirichletValue );
        }

        if ( !aConfig.mNeumannSets.empty() )
        {
            tNames.mTractionProp = tP + "PropTraction";
            tAddConstantProperty( tNames.mTractionProp, aConfig.mTraction );
            if ( !aConfig.mTractionFunction.empty() )
            {
                aFemParameterLists.set( "value_function", aConfig.mTractionFunction );
            }
        }

        // ---- constitutive model ---------------------------------------------------
        tNames.mCM = tP + "CMStrucLinIso";
        aFemParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aFemParameterLists.set( "constitutive_name", tNames.mCM );
        aFemParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aFemParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( aConfig.mDofs, "Displacement" ) );
        aFemParameterLists.set( "properties",
                tNames.mYoungsProp + ",YoungsModulus;" + tNames.mPoissonProp + ",PoissonRatio" );

        // ---- stabilization --------------------------------------------------------
        if ( !aConfig.mDirichletSets.empty() )
        {
            tNames.mNitscheSP = tP + "SPNitscheDirichlet";
            aFemParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aFemParameterLists.set( "stabilization_name", tNames.mNitscheSP );
            aFemParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
            aFemParameterLists.set( "function_parameters", std::to_string( aConfig.mNitschePenalty ) );
            aFemParameterLists.set( "leader_properties", tNames.mYoungsProp + ",Material" );
        }

        if ( !aConfig.mGhostSets.empty() )
        {
            tNames.mGhostSP = tP + "SPGhost";
            aFemParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aFemParameterLists.set( "stabilization_name", tNames.mGhostSP );
            aFemParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aFemParameterLists.set( "function_parameters", std::to_string( aConfig.mGhostPenalty ) );
            aFemParameterLists.set( "leader_properties", tNames.mYoungsProp + ",Material" );
        }

        // ---- IWGs -----------------------------------------------------------------
        tNames.mBulkIWG = tP + "IWGBulkU";
        aFemParameterLists( FEM::IWG ).add_parameter_list();
        aFemParameterLists.set( "IWG_name", tNames.mBulkIWG );
        aFemParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aFemParameterLists.set( "dof_residual", aConfig.mDofs );
        aFemParameterLists.set( "leader_dof_dependencies", aConfig.mDofs );
        aFemParameterLists.set( "leader_constitutive_models", tNames.mCM + ",ElastLinIso" );
        if ( !tNames.mBeddingProp.empty() )
        {
            aFemParameterLists.set( "leader_properties", tNames.mBeddingProp + ",Bedding" );
        }
        aFemParameterLists.set( "mesh_set_names", aConfig.mBulkSets );

        if ( !aConfig.mDirichletSets.empty() )
        {
            tNames.mDirichletIWG = tP + "IWGDirichletU";
            aFemParameterLists( FEM::IWG ).add_parameter_list();
            aFemParameterLists.set( "IWG_name", tNames.mDirichletIWG );
            aFemParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            aFemParameterLists.set( "dof_residual", aConfig.mDofs );
            aFemParameterLists.set( "leader_dof_dependencies", aConfig.mDofs );
            aFemParameterLists.set( "leader_properties", tNames.mDirichletProp + ",Dirichlet" );
            aFemParameterLists.set( "leader_constitutive_models", tNames.mCM + ",ElastLinIso" );
            aFemParameterLists.set( "stabilization_parameters", tNames.mNitscheSP + ",DirichletNitsche" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mDirichletSets );
        }

        if ( !aConfig.mNeumannSets.empty() )
        {
            tNames.mNeumannIWG = tP + "IWGTraction";
            aFemParameterLists( FEM::IWG ).add_parameter_list();
            aFemParameterLists.set( "IWG_name", tNames.mNeumannIWG );
            aFemParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
            aFemParameterLists.set( "dof_residual", aConfig.mDofs );
            aFemParameterLists.set( "leader_dof_dependencies", aConfig.mDofs );
            aFemParameterLists.set( "leader_properties", tNames.mTractionProp + ",Traction" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mNeumannSets );
        }

        if ( !aConfig.mGhostSets.empty() )
        {
            tNames.mGhostIWG = tP + "IWGGhost";
            aFemParameterLists( FEM::IWG ).add_parameter_list();
            aFemParameterLists.set( "IWG_name", tNames.mGhostIWG );
            aFemParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aFemParameterLists.set( "dof_residual", aConfig.mDofs );
            aFemParameterLists.set( "leader_dof_dependencies", aConfig.mDofs );
            aFemParameterLists.set( "follower_dof_dependencies", aConfig.mDofs );
            aFemParameterLists.set( "stabilization_parameters", tNames.mGhostSP + ",GhostSP" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mGhostSets );
        }

        // ---- standard IQIs --------------------------------------------------------
        if ( aConfig.mAddStrainEnergyIQI )
        {
            tNames.mStrainEnergyIQI = tP + "IQIBulkStrainEnergy";
            aFemParameterLists( FEM::IQI ).add_parameter_list();
            aFemParameterLists.set( "IQI_name", tNames.mStrainEnergyIQI );
            aFemParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
            aFemParameterLists.set( "leader_dof_dependencies", aConfig.mDofs );
            aFemParameterLists.set( "leader_constitutive_models", tNames.mCM + ",Elast" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mBulkSets );
        }

        if ( aConfig.mAddVolumeIQI )
        {
            tNames.mVolumeIQI = tP + "IQIBulkVolume";
            aFemParameterLists( FEM::IQI ).add_parameter_list();
            aFemParameterLists.set( "IQI_name", tNames.mVolumeIQI );
            aFemParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
            aFemParameterLists.set( "leader_properties", tNames.mDensityProp + ",Density" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mBulkSets );
        }

        return tNames;
    }

    //------------------------------------------------------------------------------------------------------------------

    /**
     * Configuration for the linear-diffusion (heat conduction) preset.
     */
    struct Diffusion_Config
    {
        std::string mBulkSets;         // required
        std::string mDirichletSets;    // "" = none
        std::string mNeumannSets;      // "" = none (heat flux)

        std::string mConductivity   = "1.0";
        std::string mDensity        = "1.0";
        std::string mHeatCapacity   = "0.0";
        std::string mDirichletValue = "0.0";
        std::string mFlux           = "1.0";

        real mNitschePenalty = 100.0;

        bool mAddVolumeIQI = true;

        std::string mPrefix = "";
    };

    struct Diffusion_Names
    {
        std::string mConductivityProp;
        std::string mDensityProp;
        std::string mHeatCapacityProp;
        std::string mDirichletProp;    // "" if not created
        std::string mFluxProp;         // "" if not created
        std::string mCM;
        std::string mNitscheSP;        // "" if not created
        std::string mBulkIWG;
        std::string mDirichletIWG;    // "" if not created
        std::string mNeumannIWG;      // "" if not created
        std::string mVolumeIQI;       // "" if not created
    };

    /**
     * Emits the standard isotropic linear-diffusion block (TEMP dof).
     *
     * @param aFemParameterLists The FEM module parameter lists to append to
     * @param aConfig Preset configuration
     * @return Names of all created entities
     */
    inline Diffusion_Names
    diffusion(
            Module_Parameter_Lists&  aFemParameterLists,
            const Diffusion_Config& aConfig )
    {
        MORIS_ERROR( !aConfig.mBulkSets.empty(),
                "fem::presets::diffusion - mBulkSets is required." );

        Diffusion_Names    tNames;
        const std::string& tP = aConfig.mPrefix;

        auto tAddConstantProperty = [ &aFemParameterLists ]( const std::string& aName, const std::string& aValue ) {
            aFemParameterLists( FEM::PROPERTIES ).add_parameter_list();
            aFemParameterLists.set( "property_name", aName );
            aFemParameterLists.set( "function_parameters", aValue );
        };

        tNames.mConductivityProp = tP + "PropConductivity";
        tAddConstantProperty( tNames.mConductivityProp, aConfig.mConductivity );

        tNames.mDensityProp = tP + "PropDensity";
        tAddConstantProperty( tNames.mDensityProp, aConfig.mDensity );

        tNames.mHeatCapacityProp = tP + "PropHeatCapacity";
        tAddConstantProperty( tNames.mHeatCapacityProp, aConfig.mHeatCapacity );

        if ( !aConfig.mDirichletSets.empty() )
        {
            tNames.mDirichletProp = tP + "PropDirichletTemp";
            tAddConstantProperty( tNames.mDirichletProp, aConfig.mDirichletValue );
        }

        if ( !aConfig.mNeumannSets.empty() )
        {
            tNames.mFluxProp = tP + "PropFlux";
            tAddConstantProperty( tNames.mFluxProp, aConfig.mFlux );
        }

        tNames.mCM = tP + "CMDiffusion";
        aFemParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aFemParameterLists.set( "constitutive_name", tNames.mCM );
        aFemParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aFemParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aFemParameterLists.set( "properties",
                tNames.mConductivityProp + ",Conductivity;"
                        + tNames.mDensityProp + ",Density;"
                        + tNames.mHeatCapacityProp + ",HeatCapacity" );

        if ( !aConfig.mDirichletSets.empty() )
        {
            tNames.mNitscheSP = tP + "SPNitscheTemp";
            aFemParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aFemParameterLists.set( "stabilization_name", tNames.mNitscheSP );
            aFemParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
            aFemParameterLists.set( "function_parameters", std::to_string( aConfig.mNitschePenalty ) );
            aFemParameterLists.set( "leader_properties", tNames.mConductivityProp + ",Material" );
        }

        tNames.mBulkIWG = tP + "IWGBulkTemp";
        aFemParameterLists( FEM::IWG ).add_parameter_list();
        aFemParameterLists.set( "IWG_name", tNames.mBulkIWG );
        aFemParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aFemParameterLists.set( "dof_residual", "TEMP" );
        aFemParameterLists.set( "leader_dof_dependencies", "TEMP" );
        aFemParameterLists.set( "leader_constitutive_models", tNames.mCM + ",Diffusion" );
        aFemParameterLists.set( "mesh_set_names", aConfig.mBulkSets );

        if ( !aConfig.mDirichletSets.empty() )
        {
            tNames.mDirichletIWG = tP + "IWGDirichletTemp";
            aFemParameterLists( FEM::IWG ).add_parameter_list();
            aFemParameterLists.set( "IWG_name", tNames.mDirichletIWG );
            aFemParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE );
            aFemParameterLists.set( "dof_residual", "TEMP" );
            aFemParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aFemParameterLists.set( "leader_properties", tNames.mDirichletProp + ",Dirichlet" );
            aFemParameterLists.set( "leader_constitutive_models", tNames.mCM + ",Diffusion" );
            aFemParameterLists.set( "stabilization_parameters", tNames.mNitscheSP + ",DirichletNitsche" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mDirichletSets );
        }

        if ( !aConfig.mNeumannSets.empty() )
        {
            tNames.mNeumannIWG = tP + "IWGFlux";
            aFemParameterLists( FEM::IWG ).add_parameter_list();
            aFemParameterLists.set( "IWG_name", tNames.mNeumannIWG );
            aFemParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_NEUMANN );
            aFemParameterLists.set( "dof_residual", "TEMP" );
            aFemParameterLists.set( "leader_dof_dependencies", "TEMP" );
            aFemParameterLists.set( "leader_properties", tNames.mFluxProp + ",Neumann" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mNeumannSets );
        }

        if ( aConfig.mAddVolumeIQI )
        {
            tNames.mVolumeIQI = tP + "IQIBulkVolume";
            aFemParameterLists( FEM::IQI ).add_parameter_list();
            aFemParameterLists.set( "IQI_name", tNames.mVolumeIQI );
            aFemParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
            aFemParameterLists.set( "leader_properties", tNames.mDensityProp + ",Density" );
            aFemParameterLists.set( "mesh_set_names", aConfig.mBulkSets );
        }

        return tNames;
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem::presets
