/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * WRK_Input_3.cpp
 *
 * Test deck of UT_WRK_Parallel_Refinement_Redo.cpp (MB-097): the 2D heat-conduction problem of WRK_Input_2.cpp on
 * [0,1]^2, 8 x 8 elements, with the union of two disk inclusions as the design (user-defined field, 3 ADVs):
 *  - a large disk at (x_c, 0.5) with radius r (ADVs 0 and 1);
 *  - a small disk at (0.5, 0.9375), the midpoint of a coarse background edge, with radius r_s (ADV 2; r_s < 0: absent).
 * With WRK_INPUT_3_GEN_REFINE=1 (read when the parameter lists are built) GEN refines the Lagrange mesh once around the
 * interface, with the MB-084 sub-cell feature refinement and the MB-089 re-detection (feature_refinement_redetect_edges
 * = 0); unset, the mesh does not depend on the design.
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
    // centre of the small disk: the midpoint of the vertical background edge x = 0.5, 0.875 <= y <= 1 (a coarse cell
    // edge that no refinement for the large disk reaches), so a small disk crosses that edge twice and no other
    static const real gSmallDiskX = 0.5;
    static const real gSmallDiskY = 0.9375;

    bool
    gen_refines()
    {
        const char* tFlag = std::getenv( "WRK_INPUT_3_GEN_REFINE" );
        return tFlag != nullptr && std::string( tFlag ) == "1";
    }

    //------------------------------------------------------------------------------

    // union of the two disks (positive inside): max of the signed distances
    real
    Two_Disk_Field(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        real tLargeDx = aCoordinates( 0 ) - aGeometryParameters( 0 );
        real tLargeDy = aCoordinates( 1 ) - 0.5;
        real tSmallDx = aCoordinates( 0 ) - gSmallDiskX;
        real tSmallDy = aCoordinates( 1 ) - gSmallDiskY;
        real tLarge   = aGeometryParameters( 1 ) - std::sqrt( tLargeDx * tLargeDx + tLargeDy * tLargeDy );
        real tSmall   = aGeometryParameters( 2 ) - std::sqrt( tSmallDx * tSmallDx + tSmallDy * tSmallDy );
        return std::max( tLarge, tSmall );
    }

    void
    Two_Disk_Field_Grad(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters,
            Matrix< DDRMat >&       aSensitivities )
    {
        real tLargeDx   = aCoordinates( 0 ) - aGeometryParameters( 0 );
        real tLargeDy   = aCoordinates( 1 ) - 0.5;
        real tSmallDx   = aCoordinates( 0 ) - gSmallDiskX;
        real tSmallDy   = aCoordinates( 1 ) - gSmallDiskY;
        real tLargeDist = std::max( std::sqrt( tLargeDx * tLargeDx + tLargeDy * tLargeDy ), 1.0e-12 );
        real tSmallDist = std::sqrt( tSmallDx * tSmallDx + tSmallDy * tSmallDy );
        aSensitivities.set_size( 1, 3, 0.0 );
        if ( aGeometryParameters( 1 ) - tLargeDist >= aGeometryParameters( 2 ) - tSmallDist )
        {
            aSensitivities( 0 ) = tLargeDx / tLargeDist;
            aSensitivities( 1 ) = 1.0;
        }
        else
        {
            aSensitivities( 2 ) = 1.0;
        }
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

        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::MMA );
        aParameterLists.set( "max_its", 1 );
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
        aParameterLists.set( "geometry_resolution_samples", 2u );
    }

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIStoredSolid", "IQIVolumeDisk" );
        aParameterLists.set( "number_of_phases", 2 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "name", "TwoDisks" );
        aParameterLists.set( "field_function_name", "Two_Disk_Field" );
        aParameterLists.set( "sensitivity_function_name", "Two_Disk_Field_Grad" );
        aParameterLists.set( "discretization_mesh_index", -2 );
        if ( gen_refines() )
        {
            aParameterLists.set( "number_of_refinements", 1 );
            aParameterLists.set( "refinement_mesh_index", 0 );
            aParameterLists.set( "feature_refinement_samples", 8u );
            aParameterLists.set( "feature_refinement_redetect_edges", 0 );
        }
        aParameterLists( GEN::GEOMETRIES ).insert( "variable_00000", Design_Variable( 0.2, 0.31, 0.8 ) );
        aParameterLists( GEN::GEOMETRIES ).insert( "variable_00001", Design_Variable( 0.1, 0.18, 0.3 ) );
        aParameterLists( GEN::GEOMETRIES ).insert( "variable_00002", Design_Variable( -0.1, 0.03, 0.1 ) );
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
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "WRK_Input_3.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", "HMR_dummy_n_p0,HMR_dummy_c_p0,HMR_dummy_n_p1,HMR_dummy_c_p1" );
        aParameterLists.set( "Field_Names", "TEMP" );
        aParameterLists.set( "Field_Type", "NODAL" );
        aParameterLists.set( "IQI_Names", "IQIStoredSolid" );
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //------------------------------------------------------------------------------
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
