/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * WRK_Input_TD_Field.cpp
 *
 * Input for UT_WRK_TD_Field (FEM FIELDS key field_elemental_normalization, IQI TOPOLOGICAL_DERIVATIVE).
 * Two-phase XFEM bar on a 4 x 2 background mesh, h = 1, interface x = 1.3 (phase 0, x < 1.3: E0 = 1;
 * phase 1: E1 = 0.25), nu = 0, plane strain. Left side clamped, unit traction in x on the right side, so the exact
 * state is piecewise constant (sigma_xx = 1, eps_xx = 1/E_p) and is reproduced by the enriched linear basis.
 * ELEMENTAL fields (written to hdf5 in the working directory):
 *   CutRatio  VOLUME on the phase-0 sets, ip_cell        -> rho_e = |e cap phase 0| / |e|
 *   TDP0      TD on the phase-0 sets (r = 0.25, delta = -1), ip_cell
 *   TDP1      TD on the phase-1 sets (r = 4, delta = +1), ip_cell
 *   TDP0Avg   TD on the phase-0 sets, default normalization (cluster)
 * (l = 0.6, chi override 2 in both TD IQIs)
 *
 */

#include <string>
#include <iostream>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "parameters.hpp"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    std::string tTDPhase0 = "HMR_dummy_n_p0,HMR_dummy_c_p0";
    std::string tTDPhase1 = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tTDAll    = "HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1";

    void
    Func_Const_TD(
            Matrix< DDRMat >&           aPropMatrix,
            Vector< Matrix< DDRMat > >& aParameters,
            fem::Field_Interpolator_Manager* )
    {
        aPropMatrix = aParameters( 0 );
    }

    real
    Plane_TD(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >& )
    {
        return aCoordinates( 0 ) - 1.3;
    }

    bool
    Output_Criterion( tsa::Time_Solver* )
    {
        return true;
    }

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 4, 2 );
        aParameterLists.set( "domain_dimensions", 4.0, 2.0 );
        aParameterLists.set( "domain_offset", 0.0, 0.0 );
        aParameterLists.set( "lagrange_output_meshes", "0" );
        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );
        aParameterLists.set( "lagrange_to_bspline", "0" );
        aParameterLists.set( "refinement_buffer", 0 );
        aParameterLists.set( "staircase_buffer", 0 );
    }

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", false );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Plane_TD" );
    }

    void
    add_property_TD(
            Module_Parameter_Lists& aParameterLists,
            const std::string&      aName,
            const std::string&      aValue )
    {
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", aName );
        aParameterLists.set( "function_parameters", aValue );
        aParameterLists.set( "value_function", "Func_Const_TD" );
    }

    void
    add_elemental_field_TD(
            Module_Parameter_Lists& aParameterLists,
            const std::string&      aFieldName,
            const std::string&      aIQIName,
            const std::string&      aNormalization )
    {
        aParameterLists( FEM::FIELDS ).add_parameter_list();
        aParameterLists.set( "field_name", aFieldName );
        aParameterLists.set( "field_entity_type", "ELEMENTAL" );
        aParameterLists.set( "field_type", "FIELD_1" );
        aParameterLists.set( "IQI_Name", aIQIName );
        aParameterLists.set( "field_output_to_file", "WRK_TD_Field_" + aFieldName + ".hdf5" );
        if ( !aNormalization.empty() )
        {
            aParameterLists.set( "field_elemental_normalization", aNormalization );
        }
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();

        add_property_TD( aParameterLists, "PropDirichlet", "0.0;0.0" );
        add_property_TD( aParameterLists, "PropTraction", "1.0;0.0" );
        add_property_TD( aParameterLists, "PropE0", "1.0" );
        add_property_TD( aParameterLists, "PropE1", "0.25" );
        add_property_TD( aParameterLists, "PropNu", "0.0" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CM0" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "model_type", fem::Model_Type::PLANE_STRAIN );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropE0,YoungsModulus;PropNu,PoissonRatio" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CM1" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "model_type", fem::Model_Type::PLANE_STRAIN );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropE1,YoungsModulus;PropNu,PoissonRatio" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPDirichlet" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropE0,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterface" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropE1,Material" );
        aParameterLists.set( "follower_properties", "PropE0,Material" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulk0" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CM0,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tTDPhase0 );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulk1" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CM1,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tTDPhase1 );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichlet" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CM0,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPDirichlet,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", "SideSet_4_n_p0" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", "SideSet_2_n_p1" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterface" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CM1,ElastLinIso" );
        aParameterLists.set( "follower_constitutive_models", "CM0,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPInterface,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", "dbl_iside_p0_1_p1_0" );

        // cut ratio
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQICutRatio" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "mesh_set_names", tTDPhase0 );

        // topological derivative, material set ( r = E1/E0, delta = -1, l, chi = 2 )
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQITDP0" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TOPOLOGICAL_DERIVATIVE );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CM0,Elast" );
        aParameterLists.set( "function_parameters", "0.25/-1.0/0.6/2.0" );
        aParameterLists.set( "mesh_set_names", tTDPhase0 );

        // topological derivative, void set ( r = E0/E1, delta = +1, l, chi = 2 )
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQITDP1" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::TOPOLOGICAL_DERIVATIVE );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CM1,Elast" );
        aParameterLists.set( "function_parameters", "4.0/1.0/0.6/2.0" );
        aParameterLists.set( "mesh_set_names", tTDPhase1 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tTDAll );

        add_elemental_field_TD( aParameterLists, "CutRatio", "IQICutRatio", "ip_cell" );
        add_elemental_field_TD( aParameterLists, "TDP0", "IQITDP0", "ip_cell" );
        add_elemental_field_TD( aParameterLists, "TDP1", "IQITDP1", "ip_cell" );
        add_elemental_field_TD( aParameterLists, "TDP0Avg", "IQITDP0", "" );

        aParameterLists( FEM::COMPUTATION );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aParameterLists.set( "Solver_Type", "Amesos_Klu" );
        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_max_iter", 1 );
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& )
    {
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "WRK_TD_Field.exo" ) );
        aParameterLists.set( "Set_Names", tTDAll );
        aParameterLists.set( "Field_Names", "UX" );
        aParameterLists.set( "Field_Type", "NODAL" );
        aParameterLists.set( "IQI_Names", "IQIUX" );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
