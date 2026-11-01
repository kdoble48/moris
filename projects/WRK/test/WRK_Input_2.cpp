/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * WRK_Input_2.cpp
 *
 * Test deck of UT_WRK_Reinitialize_User_Function.cpp (sdf_reinit_mode = "user_function"): a 2D heat-conduction
 * problem on [0,1]^2 with a disk inclusion whose centre x and radius are the two ADVs (user-defined field), optimized by
 * the native MMA. The OPT reinitialization schedule (reinitialize_interface_iter = 3) calls Reinit_Known_Transform(),
 * which moves the disk to a known design. compute_objectives() records every evaluated ADV vector, so that the test can
 * check that the optimizer continues from exactly the function's result.
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
    // the known transform: the design the reinitialization moves to (centre x, radius)
    static const real gReinitCentreX = 0.43;
    static const real gReinitRadius  = 0.27;

    // record of the reinitialization calls and of the evaluated designs (read by the test through the getters below)
    static uint                     gNumReinitCalls      = 0;
    static uint                     gLastReinitIteration = 0;
    static Vector< real >           gLastReinitInput;
    static Vector< Vector< real > > gEvaluatedADVs;
    static Vector< uint >           gReinitAfterEvaluation;    // number of recorded evaluations when each call happened

    // WRK_INPUT_2_GEN_REFINE=1 (read when the parameter lists are built): GEN refines the Lagrange mesh around the disk, so
    // a restart rebuilds the mesh on the current design (Workflow_HMR_XTK::refine_hmr_on_current_design(), the GS-LS decks'
    // path); unset: the mesh does not depend on the design
    bool
    gen_refines()
    {
        const char* tFlag = std::getenv( "WRK_INPUT_2_GEN_REFINE" );
        return tFlag != nullptr && std::string( tFlag ) == "1";
    }

    //------------------------------------------------------------------------------

    real
    Disk_Field(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        real tDx = aCoordinates( 0 ) - aGeometryParameters( 0 );
        real tDy = aCoordinates( 1 ) - 0.5;
        return aGeometryParameters( 1 ) - std::sqrt( tDx * tDx + tDy * tDy );
    }

    void
    Disk_Field_Grad(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters,
            Matrix< DDRMat >&       aSensitivities )
    {
        real tDx   = aCoordinates( 0 ) - aGeometryParameters( 0 );
        real tDy   = aCoordinates( 1 ) - 0.5;
        real tDist = std::max( std::sqrt( tDx * tDx + tDy * tDy ), 1.0e-12 );
        aSensitivities.set_size( 1, 2 );
        aSensitivities( 0 ) = tDx / tDist;
        aSensitivities( 1 ) = 1.0;
    }

    //------------------------------------------------------------------------------

    // sdf_reinit_mode = "user_function": the known transform
    Vector< real >
    Reinit_Known_Transform( const Vector< real >& aADVs, uint aIteration )
    {
        gNumReinitCalls++;
        gLastReinitIteration = aIteration;
        gLastReinitInput     = aADVs;
        gReinitAfterEvaluation.push_back( gEvaluatedADVs.size() );

        Vector< real > tNew( aADVs.size() );
        tNew( 0 ) = gReinitCentreX;
        tNew( 1 ) = gReinitRadius;
        return tNew;
    }

    // a failing function: returns no ADVs (the core must stop the run)
    Vector< real >
    Reinit_Fails( const Vector< real >& aADVs, uint aIteration )
    {
        return Vector< real >( 0 );
    }

    // getters for the test (the library may stay loaded between test sections: reset first)
    void
    reset_record()
    {
        gNumReinitCalls        = 0;
        gLastReinitIteration   = 0;
        gLastReinitInput       = Vector< real >( 0 );
        gEvaluatedADVs         = Vector< Vector< real > >( 0 );
        gReinitAfterEvaluation = Vector< uint >( 0 );
    }

    uint
    get_num_reinit_calls()
    {
        return gNumReinitCalls;
    }

    uint
    get_last_reinit_iteration()
    {
        return gLastReinitIteration;
    }

    Vector< real >
    get_last_reinit_input()
    {
        return gLastReinitInput;
    }

    uint
    get_num_evaluations()
    {
        return gEvaluatedADVs.size();
    }

    Vector< real >
    get_evaluated_advs( uint aIndex )
    {
        return gEvaluatedADVs( aIndex );
    }

    uint
    get_reinit_after_evaluation( uint aCall )
    {
        return gReinitAfterEvaluation( aCall );
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

    // objective: stored thermal measure of the solid (criterion 0); constraint: disk volume below 0.3 (criterion 1)
    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        bool tFinite = true;
        for ( uint i = 0; i < aCriteria.size(); i++ )
        {
            tFinite = tFinite && std::isfinite( aCriteria( i ) );
        }
        // record the evaluated designs only (the NaN call that requests the restart is not an evaluation)
        if ( tFinite )
        {
            gEvaluatedADVs.push_back( aADVs );
        }
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
        aParameterLists.set( "sdf_reinit_mode", "user_function" );
        aParameterLists.set( "reinit_user_function", "Reinit_Known_Transform" );
        aParameterLists.set( "reinitialize_interface_iter", 3 );

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::MMA );
        aParameterLists.set( "max_its", 6 );
        aParameterLists.set( "step_size", 0.05 );
        aParameterLists.set( "max_inner_its", 0 );
        aParameterLists.set( "norm_drop", 1.0e-12 );
    }

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", 8, 8 );
        aParameterLists.set( "domain_dimensions", 1.0, 1.0 );
        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );
        if ( gen_refines() )
        {
            // the GS-LS decks' layout: the Lagrange mesh (pattern 0) is refined around the design, the B-spline mesh
            // (pattern 1) never
            aParameterLists.set( "pattern_initial_refinement", 0, 0 );
        }

        aParameterLists( HMR::LAGRANGE_MESHES ).add_parameter_list();
        aParameterLists.set( "order", 1 );

        aParameterLists( HMR::BSPLINE_MESHES ).add_parameter_list();
        aParameterLists.set( "orders", 1 );
        if ( gen_refines() )
        {
            aParameterLists.set( "pattern_index", 1 );
        }
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
        aParameterLists.set( "name", "Disk" );
        aParameterLists.set( "field_function_name", "Disk_Field" );
        aParameterLists.set( "sensitivity_function_name", "Disk_Field_Grad" );
        aParameterLists.set( "discretization_mesh_index", -2 );
        if ( gen_refines() )
        {
            aParameterLists.set( "number_of_refinements", 1 );
            aParameterLists.set( "refinement_mesh_index", 0 );
        }
        aParameterLists( GEN::GEOMETRIES ).insert( "variable_00000", Design_Variable( 0.3, 0.5, 0.7 ) );
        aParameterLists( GEN::GEOMETRIES ).insert( "variable_00001", Design_Variable( 0.15, 0.31, 0.4 ) );
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
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "WRK_Input_2.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1" );
        aParameterLists.set( "Field_Names", "TEMP" );
        aParameterLists.set( "Field_Type", "NODAL" );
        aParameterLists.set( "IQI_Names", "IQIStoredSolid" );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // mapping sublist (index 2): the design the user_function reinitialization replaces
        aParameterLists( 2 ).add_parameter_list( Parameter_List( "" ) );
        prm::create_mapping_parameterlist( aParameterLists( 2 )( 0 ) );
        aParameterLists.set( "adv_field", "Disk" );
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
