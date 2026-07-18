/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Presets.cpp
 *
 * Structural tests for the input-deck physics presets (fn_FEM_Presets.hpp): the
 * preset must emit the same property -> CM -> SP -> IWG -> IQI block a hand-written
 * deck carries, with consistent name wiring. (Enum-valued parameters are exercised
 * end-to-end by the EXA twin decks; here the structure and string links are pinned.)
 */

#include <catch.hpp>

#include "cl_Module_Parameter_Lists.hpp"
#include "fn_FEM_Presets.hpp"
#include "parameters.hpp"

namespace moris
{
    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "fem::presets::linear_elastic emits the standard block", "[moris],[fem],[presets]" )
    {
        Module_Parameter_Lists tFem( Module_Type::FEM );

        fem::presets::Linear_Elastic_Config tConfig;
        tConfig.mBulkSets       = "HMR_dummy_n_p1,HMR_dummy_c_p1";
        tConfig.mDirichletSets  = "SideSet_4_n_p1,SideSet_4_c_p1";
        tConfig.mNeumannSets    = "SideSet_2_n_p1,SideSet_2_c_p1";
        tConfig.mGhostSets      = "ghost_p1";
        tConfig.mBedding        = "1.0e-6";
        tConfig.mDirichletValue = "0.0;0.0";
        tConfig.mTraction       = "0.0;1.0";

        fem::presets::Linear_Elastic_Names tNames = fem::presets::linear_elastic( tFem, tConfig );

        SECTION( "counts: 6 properties, 1 CM, 2 SPs, 4 IWGs, 2 IQIs" )
        {
            CHECK( tFem( FEM::PROPERTIES ).size() == 6 );
            CHECK( tFem( FEM::CONSTITUTIVE_MODELS ).size() == 1 );
            CHECK( tFem( FEM::STABILIZATION ).size() == 2 );
            CHECK( tFem( FEM::IWG ).size() == 4 );
            CHECK( tFem( FEM::IQI ).size() == 2 );
        }

        SECTION( "constant properties omit value_function (builtin constant applies)" )
        {
            REQUIRE( tFem( FEM::PROPERTIES ).size() >= 2 );
            CHECK( tFem( FEM::PROPERTIES )( 0 ).get< std::string >( "property_name" ) == tNames.mDensityProp );
            CHECK( tFem( FEM::PROPERTIES )( 1 ).get< std::string >( "property_name" ) == tNames.mYoungsProp );
            CHECK( tFem( FEM::PROPERTIES )( 1 ).get< std::string >( "function_parameters" ) == "1.0" );
            CHECK( tFem( FEM::PROPERTIES )( 1 ).get< std::string >( "value_function" ).empty() );
        }

        SECTION( "cross-references are wired by name" )
        {
            // CM references the created properties
            CHECK( tFem( FEM::CONSTITUTIVE_MODELS )( 0 ).get< std::string >( "constitutive_name" ) == tNames.mCM );
            CHECK( tFem( FEM::CONSTITUTIVE_MODELS )( 0 ).get< std::string >( "properties" )
                    == tNames.mYoungsProp + ",YoungsModulus;" + tNames.mPoissonProp + ",PoissonRatio" );

            // bulk IWG references the CM and the bedding property
            CHECK( tFem( FEM::IWG )( 0 ).get< std::string >( "IWG_name" ) == tNames.mBulkIWG );
            CHECK( tFem( FEM::IWG )( 0 ).get< std::string >( "leader_constitutive_models" ) == tNames.mCM + ",ElastLinIso" );
            CHECK( tFem( FEM::IWG )( 0 ).get< std::string >( "leader_properties" ) == tNames.mBeddingProp + ",Bedding" );
            CHECK( tFem( FEM::IWG )( 0 ).get< std::string >( "mesh_set_names" ) == tConfig.mBulkSets );

            // Dirichlet IWG references the Nitsche SP
            CHECK( tFem( FEM::IWG )( 1 ).get< std::string >( "stabilization_parameters" ) == tNames.mNitscheSP + ",DirichletNitsche" );

            // ghost IWG references the ghost SP and has follower dofs
            CHECK( tFem( FEM::IWG )( 3 ).get< std::string >( "stabilization_parameters" ) == tNames.mGhostSP + ",GhostSP" );
            CHECK( tFem( FEM::IWG )( 3 ).get< std::string >( "follower_dof_dependencies" ) == "UX,UY" );

            // IQIs reference the CM / density property on the bulk sets
            CHECK( tFem( FEM::IQI )( 0 ).get< std::string >( "IQI_name" ) == tNames.mStrainEnergyIQI );
            CHECK( tFem( FEM::IQI )( 0 ).get< std::string >( "leader_constitutive_models" ) == tNames.mCM + ",Elast" );
            CHECK( tFem( FEM::IQI )( 1 ).get< std::string >( "leader_properties" ) == tNames.mDensityProp + ",Density" );
        }

        SECTION( "optional pieces are dropped when not configured" )
        {
            Module_Parameter_Lists tFemMinimal( Module_Type::FEM );

            fem::presets::Linear_Elastic_Config tMinimalConfig;
            tMinimalConfig.mBulkSets          = "HMR_dummy_n_p1";
            tMinimalConfig.mAddVolumeIQI      = false;
            tMinimalConfig.mAddStrainEnergyIQI = false;

            fem::presets::Linear_Elastic_Names tMinimalNames = fem::presets::linear_elastic( tFemMinimal, tMinimalConfig );

            CHECK( tFemMinimal( FEM::PROPERTIES ).size() == 3 );      // density, youngs, poisson
            CHECK( tFemMinimal( FEM::STABILIZATION ).size() == 0 );
            CHECK( tFemMinimal( FEM::IWG ).size() == 1 );             // bulk only
            CHECK( tFemMinimal( FEM::IQI ).size() == 0 );
            CHECK( tMinimalNames.mDirichletIWG.empty() );
            CHECK( tMinimalNames.mGhostSP.empty() );
        }

        SECTION( "prefix keeps two preset instances collision-free" )
        {
            Module_Parameter_Lists tFemTwoPhase( Module_Type::FEM );

            fem::presets::Linear_Elastic_Config tPhaseOne;
            tPhaseOne.mBulkSets = "HMR_dummy_n_p1";
            tPhaseOne.mPrefix   = "P1";

            fem::presets::Linear_Elastic_Config tPhaseTwo;
            tPhaseTwo.mBulkSets = "HMR_dummy_n_p2";
            tPhaseTwo.mPrefix   = "P2";

            auto tNamesOne = fem::presets::linear_elastic( tFemTwoPhase, tPhaseOne );
            auto tNamesTwo = fem::presets::linear_elastic( tFemTwoPhase, tPhaseTwo );

            CHECK( tNamesOne.mCM == "P1CMStrucLinIso" );
            CHECK( tNamesTwo.mCM == "P2CMStrucLinIso" );
            CHECK( tFemTwoPhase( FEM::CONSTITUTIVE_MODELS ).size() == 2 );
        }
    }

    //------------------------------------------------------------------------------------------------------------------

    TEST_CASE( "fem::presets::diffusion emits the standard block", "[moris],[fem],[presets]" )
    {
        Module_Parameter_Lists tFem( Module_Type::FEM );

        fem::presets::Diffusion_Config tConfig;
        tConfig.mBulkSets      = "HMR_dummy_n_p0,HMR_dummy_c_p0";
        tConfig.mDirichletSets = "SideSet_4_n_p0";
        tConfig.mNeumannSets   = "SideSet_2_n_p0";

        fem::presets::Diffusion_Names tNames = fem::presets::diffusion( tFem, tConfig );

        CHECK( tFem( FEM::PROPERTIES ).size() == 5 );
        CHECK( tFem( FEM::CONSTITUTIVE_MODELS ).size() == 1 );
        CHECK( tFem( FEM::STABILIZATION ).size() == 1 );
        CHECK( tFem( FEM::IWG ).size() == 3 );
        CHECK( tFem( FEM::IQI ).size() == 1 );

        CHECK( tFem( FEM::CONSTITUTIVE_MODELS )( 0 ).get< std::string >( "properties" )
                == tNames.mConductivityProp + ",Conductivity;"
                        + tNames.mDensityProp + ",Density;"
                        + tNames.mHeatCapacityProp + ",HeatCapacity" );

        CHECK( tFem( FEM::IWG )( 2 ).get< std::string >( "leader_properties" ) == tNames.mFluxProp + ",Neumann" );
    }

    //------------------------------------------------------------------------------------------------------------------

}    // namespace moris
