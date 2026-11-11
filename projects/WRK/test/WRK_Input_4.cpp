/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * WRK_Input_4.cpp
 *
 * Test deck of UT_WRK_Parallel_Redistance.cpp (MB-113): a 2D heat-conduction problem on [0,1]^2, 16 x 16 elements. The
 * design is a level set on linear B-splines (one ADV per coefficient), started from an ellipse inclusion whose field
 * is deliberately not a distance field. OPT asks for the interface redistance in place at every evaluation
 * (sdf_reinit_mode = "interface", sdf_reinit_every_step) with the footprint diagnostics on. Read when the parameter
 * lists are built:
 *  - WRK_INPUT_4_PRESERVE=1: sdf_preserve_interface;
 *  - WRK_INPUT_4_SMOOTHING=<n>: n Laplacian smoothing iterations of the distance values.
 *
 */

#include <string>
#include <iostream>
#include <cmath>
#include <cstdlib>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Bitset.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    bool
    env_flag( const char* aName )
    {
        const char* tFlag = std::getenv( aName );
        return tFlag != nullptr && std::string( tFlag ) == "1";
    }

    sint
    env_int( const char* aName )
    {
        const char* tValue = std::getenv( aName );
        return tValue != nullptr ? std::atoi( tValue ) : 0;
    }

    //------------------------------------------------------------------------------

    // ellipse inclusion (positive inside), bounded and not a distance field: its slope varies along the interface
    real
    Ellipse_Field(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        real tX = ( aCoordinates( 0 ) - 0.47 ) / 0.31;
        real tY = ( aCoordinates( 1 ) - 0.52 ) / 0.23;
        return 0.25 * std::tanh( 2.0 * ( 1.0 - tX * tX - tY * tY ) );
    }

    //------------------------------------------------------------------------------

    void
    Func_Const(
            Matrix< DDRMat >&                aPropMatrix,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    bool
    Output_Criterion( tsa::Time_Solver* aTimeSolver )
    {
        return false;
    }

    //------------------------------------------------------------------------------

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );
        return tConstraintTypes;
    }

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives( 1, 1 );
        tObjectives( 0 ) = aCriteria( 0 );
        return tObjectives;
    }

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        return Matrix< DDRMat >( 1, aADVs.size(), 0.0 );
    }

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tD( 1, aCriteria.size(), 0.0 );
        tD( 0, 0 ) = 1.0;
        return tD;
    }

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints( 1, 1 );
        tConstraints( 0 ) = aCriteria( 1 ) / 0.3 - 1.0;
        return tConstraints;
    }

    Matrix< DDRMat >
    compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        return Matrix< DDRMat >( 1, aADVs.size(), 0.0 );
    }

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tD( 1, aCriteria.size(), 0.0 );
        tD( 0, 1 ) = 1.0 / 0.3;
        return tD;
    }

    //------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );

        // interface redistance in place at every evaluation
        aParameterLists.set( "sdf_reinit_mode", "interface" );
        aParameterLists.set( "sdf_reinit_every_step", true );
        aParameterLists.set( "sdf_reinit_diagnostics", true );
        aParameterLists.set( "sdf_preserve_interface", env_flag( "WRK_INPUT_4_PRESERVE" ) );
        aParameterLists.set( "sdf_smoothing_iterations", env_int( "WRK_INPUT_4_SMOOTHING" ) );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::MMA );
        aParameterLists.set( "max_its", 1 );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 16, 16 );
        aParameterLists.set( "domain_dimensions", 1.0, 1.0 );
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );

        aParameterLists( HMR::LAGRANGE_MESHES ).add_parameter_list();
        aParameterLists.set( "order", 1 );

        aParameterLists( HMR::BSPLINE_MESHES ).add_parameter_list();
        aParameterLists.set( "orders", 1 );
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
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIStoredSolid", "IQIVolumeDisk" );
        aParameterLists.set( "number_of_phases", 2 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "name", "Level_Set_Field" );
        aParameterLists.set( "field_function_name", "Ellipse_Field" );
        aParameterLists.set( "discretization_mesh_index", 0 );
        aParameterLists.set( "discretization_lower_bound", -0.3 );
        aParameterLists.set( "discretization_upper_bound", 0.3 );
    }

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "phase_indices", "0" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseDisk" );
        aParameterLists.set( "phase_indices", "1" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivityDisk" );
        aParameterLists.set( "function_parameters", "10.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropHeatLoad" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichlet" );
        aParameterLists.set( "function_parameters", "0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropOne" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion" );
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties", "PropConductivity,Conductivity" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusionDisk" );
        aParameterLists.set( "phase_name", "PhaseDisk" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties", "PropConductivityDisk,Conductivity" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitscheT" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "follower_phase_name", "PhaseDisk" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );
        aParameterLists.set( "follower_properties", "PropConductivityDisk,Material" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkSolid" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropHeatLoad,Load" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkDisk" );
        aParameterLists.set( "leader_phase_name", "PhaseDisk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropHeatLoad,Load" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusionDisk,Diffusion" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterface" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::DOUBLE_SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "follower_phase_name", "PhaseDisk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_INTERFACE_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "follower_constitutive_models", "CMDiffusionDisk,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichlet" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "4" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_UNSYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropDirichlet,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitscheT,DirichletNitsche" );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIStoredSolid" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIVolumeDisk" );
        aParameterLists.set( "leader_phase_name", "PhaseDisk" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropOne,Density" );

        aParameterLists( FEM::COMPUTATION );
        aParameterLists.set( "print_physics_model", false );
    }

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aParameterLists.set( "Solver_Type", "Amesos_Klu" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-12 );
        aParameterLists.set( "NLA_max_iter", 3 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );
        aParameterLists.set( "NLA_Solver_Implementation", NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_DofTypes", "TEMP" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 0 );
        aParameterLists.set( "TSA_Nonlinear_Sensitivity_Solver", 0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "TEMP" );
        aParameterLists.set( "TSA_Initialize_Sol_Vec", "TEMP,0.0" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );
        aParameterLists.set( "TSA_time_level_per_type", "TEMP,1" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "TEMP", 0 );
    }

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "WRK_Input_4.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1" );
        aParameterLists.set( "Field_Names", "TEMP" );
        aParameterLists.set( "Field_Type", "NODAL" );
        aParameterLists.set( "IQI_Names", "IQIStoredSolid" );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // no remeshing (0) and no refinement (1); the mapping sublist (2) creates the Reinitialize_Performer
        aParameterLists( 2 ).add_parameter_list( Parameter_List( "" ) );
        prm::create_mapping_parameterlist( aParameterLists( 2 )( 0 ) );
        aParameterLists.set( "adv_field", "Level_Set_Field" );
        aParameterLists.set( "reinitialization_frequency", 1 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
