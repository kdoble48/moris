/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Shape_Sensitivity_Sweep_NewIO.cpp
 *
 * Single-entry-point (Deck API v2) twin of Shape_Sensitivity_Sweep: the same
 * problem, gated by the same analytical-vs-finite-difference sensitivity checks,
 * written in the MORISInputDeck style (see doc/internal/DECK_API_RFC.md).
 *
 * Compared to the legacy twin (416 lines): no extern "C" module functions, no
 * Func_Const / Output_Criterion (builtins), 4 OPT callbacks instead of 7 (the
 * constraint types come from the 'constraint_types' parameter and the explicit
 * ADV gradients from the built-in zero defaults) registered in-process with
 * internal linkage, and the FEM block emitted by the linear-elastic preset.
 */

#include <string>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "parameters.hpp"
#include "cl_Input_Deck.hpp"
#include "fn_FEM_Presets.hpp"

namespace
{
    using namespace moris;

    const std::string tMeshSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    const std::string tDBCSets  = "iside_b0_1_b1_3";
    const std::string tNBCSets  = "iside_b0_1_b1_0";

    // FD in adjoint
    const real tFEMFdEpsilon = 1.0e-7;

    // FD step size in sweep
    const std::string tFDsweep = "1.0e-7";

    //------------------------------------------------------------------------------
    // OPT callbacks — internal linkage; resolved through the in-process registry

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        return { { aCriteria( 0 ) } };
    }

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        return { { aCriteria( 1 ) } };
    }

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        return { { 1.0, 0.0 } };
    }

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        return { { 0.0, 1.0 } };
    }
}    // namespace

//------------------------------------------------------------------------------

MORIS_DECK( aDeck )
{
    using namespace moris;

    // ---- OPT: sweep-based FD verification of the shape sensitivities ------------
    auto& tOpt = aDeck.opt();
    tOpt.set( "is_optimization_problem", true );
    tOpt.set( "problem", "user_defined" );
    tOpt.set( "constraint_types", "1" );    // replaces get_constraint_types()
    // NOTE: no 'library' parameter — the loaded deck is handed through to OPT

    tOpt( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
    tOpt.set( "hdf5_path", "shape_opt_test.hdf5" );
    tOpt.set( "evaluate_objective_gradients", true );
    tOpt.set( "evaluate_constraint_gradients", true );
    tOpt.set( "num_evaluations_per_adv", "1" );
    tOpt.set( "include_bounds", false );
    tOpt.set( "finite_difference_type", "all" );
    tOpt.set( "finite_difference_epsilons", tFDsweep );

    aDeck.register_function( "compute_objectives", &compute_objectives );
    aDeck.register_function( "compute_constraints", &compute_constraints );
    aDeck.register_function( "compute_dobjective_dcriteria", &compute_dobjective_dcriteria );
    aDeck.register_function( "compute_dconstraint_dcriteria", &compute_dconstraint_dcriteria );

    // ---- HMR: background mesh ----------------------------------------------------
    auto& tHmr = aDeck.hmr();
    tHmr.set( "number_of_elements_per_dimension", 2, 2 );
    tHmr.set( "domain_dimensions", 2.0, 2.0 );
    tHmr.set( "domain_offset", -1.0, -1.0 );
    tHmr.set( "lagrange_output_meshes", "0" );
    tHmr.set( "lagrange_orders", "1" );
    tHmr.set( "lagrange_pattern", "0" );
    tHmr.set( "bspline_orders", "1" );
    tHmr.set( "bspline_pattern", "0" );
    tHmr.set( "lagrange_to_bspline", "0" );
    tHmr.set( "refinement_buffer", 3 );
    tHmr.set( "staircase_buffer", 3 );
    tHmr.set( "adaptive_refinement_level", 1 );

    // ---- XTK: conformal cut ------------------------------------------------------
    auto& tXtk = aDeck.xtk();
    tXtk.set( "decompose", true );
    tXtk.set( "decomposition_type", "conformal" );
    tXtk.set( "enrich_mesh_indices", "0" );
    tXtk.set( "ghost_stab", true );
    tXtk.set( "multigrid", false );
    tXtk.set( "verbose", false );
    tXtk.set( "print_enriched_ig_mesh", false );
    tXtk.set( "exodus_output_XTK_ig_mesh", true );
    tXtk.set( "high_to_low_dbl_side_sets", true );

    // ---- GEN: two design lines; criteria = (strain energy, volume) --------------
    auto& tGen = aDeck.gen();
    tGen.set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );

    // vertical line (x-position and normal are design variables)
    tGen( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
    tGen.set( "center_x", 0.8, 0.8, 0.8 );
    tGen.set( "center_y", 0.3 );
    tGen.set( "normal_x", 1.0 );
    tGen.set( "normal_y", 0.0, 0.0, 0.0 );

    // oblique line
    tGen( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
    tGen.set( "center_x", -0.6, -0.6, -0.6 );
    tGen.set( "center_y", -0.3, -0.3, -0.3 );
    tGen.set( "normal_x", 1.0, 1.0, 1.0 );
    tGen.set( "normal_y", 1.0, 1.0, 1.0 );

    // ---- FEM: linear elasticity via the preset -----------------------------------
    auto& tFem = aDeck.fem();
    tFem.hack_for_legacy_fem();

    fem::presets::Linear_Elastic_Config tElasticity;
    tElasticity.mBulkSets      = tMeshSets;
    tElasticity.mDirichletSets = tDBCSets;
    tElasticity.mNeumannSets   = tNBCSets;
    tElasticity.mTraction      = "0.0;0.1";
    fem::presets::linear_elastic( tFem, tElasticity );

    // displacement IQIs for the output mesh (raw parameter-list escape hatch)
    tFem( FEM::IQI ).add_parameter_list();
    tFem.set( "IQI_name", "IQIDispX" );
    tFem.set( "IQI_type", fem::IQI_Type::DOF );
    tFem.set( "dof_quantity", "UX,UY" );
    tFem.set( "leader_dof_dependencies", "UX,UY" );
    tFem.set( "vectorial_field_index", 0 );
    tFem.set( "mesh_set_names", tMeshSets );

    tFem( FEM::IQI ).add_parameter_list();
    tFem.set( "IQI_name", "IQIDispY" );
    tFem.set( "IQI_type", fem::IQI_Type::DOF );
    tFem.set( "dof_quantity", "UX,UY" );
    tFem.set( "leader_dof_dependencies", "UX,UY" );
    tFem.set( "vectorial_field_index", 1 );
    tFem.set( "mesh_set_names", tMeshSets );

    // FD scheme for the analytical sensitivity checks in FEM
    tFem( FEM::COMPUTATION );
    tFem.set( "finite_difference_scheme", fem::FDScheme_Type::POINT_3_CENTRAL );
    tFem.set( "finite_difference_perturbation_size", tFEMFdEpsilon );

    // ---- SOL: Belos + ILU, single static step, adjoint sensitivities ------------
    auto& tSol = aDeck.sol();
    tSol( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::BELOS_IMPL );
    tSol.set( "preconditioners", "0" );
    tSol( SOL::LINEAR_SOLVERS ).add_parameter_list();
    tSol( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
    tSol( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
    tSol.set( "NLA_DofTypes", "UX,UY" );
    tSol( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
    tSol.set( "TSA_Num_Time_Steps", 1 );
    tSol.set( "TSA_Time_Frame", 1.0 );
    tSol( SOL::TIME_SOLVERS ).add_parameter_list();
    tSol.set( "TSA_DofTypes", "UX,UY" );
    tSol.set( "TSA_Output_Indices", "0" );
    tSol.set( "TSA_Output_Criteria", "Output_Criterion" );    // built-in always-true criterion
    tSol( SOL::SOLVER_WAREHOUSE ).set( "Sensitivity_Analysis_Type", sol::SensitivityAnalysisType::ADJOINT );
    tSol( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::IFPACK );
    tSol.set( "ifpack_prec_type", "ILU" );

    // ---- MSI: defaults -----------------------------------------------------------
    aDeck.msi();

    // ---- VIS: displacement output ------------------------------------------------
    auto& tVis = aDeck.vis();
    tVis.set( "File_Name", std::pair< std::string, std::string >( "./", "shape_sensitivities_newio.exo" ) );
    tVis.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
    tVis.set( "Set_Names", tMeshSets );
    tVis.set( "Field_Names", "UX,UY" );
    tVis.set( "Field_Type", "NODAL,NODAL" );
    tVis.set( "IQI_Names", "IQIDispX,IQIDispY" );
    tVis.set( "Save_Frequency", 1 );

    // ---- MORISGENERAL: defaults --------------------------------------------------
    aDeck.morisgeneral();
}
