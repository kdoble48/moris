/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * LevelSet_Boxbeam_NewIO.cpp
 *
 * Single-entry-point (Deck API v2) twin of Levelset_Boxbeam: the same two-phase
 * level-set topology-optimization problem (superellipse frame + hole-seeded
 * interior, GCMMA), written in the MORISInputDeck style (see
 * doc/internal/DECK_API_RFC.md).
 *
 * Compared to the legacy twin (616 lines): no extern "C" module functions, no
 * Func_Const / Output_Criterion / get_constraint_types (builtins + parameter),
 * ZERO OPT callback functions — the objective and constraint are criteria
 * expressions (criteria gradients by reverse-mode differentiation, GEN IQI_types
 * and constraint types derived automatically), the mesh-set strings are generated
 * by the typed vocabulary, and the per-phase FEM blocks are emitted by the
 * linear-elastic preset. The one remaining callback (the spatially limited
 * traction) is registered in-process — no dlsym, no extern "C".
 *
 * BOXBEAM_MAX_ITS overrides the GCMMA iteration count (default 2, the CI gate);
 * smoke runs set e.g. BOXBEAM_MAX_ITS=10.
 */

#include <cmath>
#include <cstdlib>
#include <string>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "parameters.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_stringify_matrix.hpp"
#include "cl_Input_Deck.hpp"
#include "cl_Input_Deck_Vocabulary.hpp"
#include "fn_FEM_Presets.hpp"

namespace
{
    using namespace moris;

    // ---- main problem parameters (identical to the legacy twin) -----------------
    const real tInitialStrainEnergy = 3.26439;
    const real tInitialPerimeter    = 29.0669;
    const real tPerimeterPenalty    = 0.2;

    const real tMaxMass = 1.0;

    const real tMMAPenalty  = 5.0;
    const real tMMAStepSize = 0.01;

    const int tInitialRef = 2;

    const real tElementEdgeLength = 1.0 / 15.0 / std::pow( 2, tInitialRef );
    const real tLoadLimitY        = std::floor( 0.2 / tElementEdgeLength ) * tElementEdgeLength;

    const int  tNumHoleX     = 15;
    const int  tNumHoleY     = 5;
    const real tHoleDiameter = 0.065;

    const real tWallThickness = 0.05;

    const real tBsplineLimit = 0.06;

    const bool tUseGhost = true;

    // GCMMA iterations: CI gate default, overridable for smoke/production runs
    int
    max_iterations()
    {
        const char* tEnvMaxIts = std::getenv( "BOXBEAM_MAX_ITS" );
        return tEnvMaxIts != nullptr ? std::atoi( tEnvMaxIts ) : 2;
    }

    // GCMMA step size: legacy default, overridable for smoke tuning
    real
    step_size()
    {
        const char* tEnvStepSize = std::getenv( "BOXBEAM_STEP_SIZE" );
        return tEnvStepSize != nullptr ? std::atof( tEnvStepSize ) : tMMAStepSize;
    }

    // ---- phases and generated set names -----------------------------------------
    // geometry sign pattern -> bulk phase via the phase table [0,1,2,2]:
    // 0 = void, 1 = interior (design), 2 = frame
    const deck::Phase tVoid( 0 );
    const deck::Phase tInterior( 1 );
    const deck::Phase tFrame( 2 );

    const std::string tFrameSets       = tFrame.bulk();                                   // HMR_dummy_n_p2,HMR_dummy_c_p2
    const std::string tInteriorSets    = tInterior.bulk();                                // HMR_dummy_n_p1,HMR_dummy_c_p1
    const std::string tTotalDomainSets = deck::join( tFrameSets, tInteriorSets );

    const std::string tFrameLoadSSets    = tFrame.side( deck::Side::Right );              // SideSet_2_*_p2
    const std::string tFrameSupportSSets = tFrame.side( deck::Side::Left );               // SideSet_4_*_p2

    const std::string tInterfaceVoidSSets =
            deck::join( deck::interface( tFrame, tVoid ), deck::interface( tInterior, tVoid ) );

    const std::string tFrameInteriorDSets = deck::between( tInterior, tFrame );           // dbl_iside_p0_1_p1_2

    //------------------------------------------------------------------------------
    // Traction: load only the lower part (y < tLoadLimitY) of the right edge.
    // Registered in-process by the deck — the one genuinely custom callback.
    void
    Func_Neumann_U(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        if ( aFIManager->get_IG_geometry_interpolator()->valx()( 1 ) < tLoadLimitY )
        {
            aPropMatrix = { { 0.0 }, { aParameters( 0 )( 0 ) } };
        }
        else
        {
            aPropMatrix = { { 0.0 }, { 0.0 } };
        }
    }

}    // namespace

//------------------------------------------------------------------------------

MORIS_DECK( aDeck )
{
    using namespace moris;

    aDeck.register_function( "Func_Neumann_U", &Func_Neumann_U );

    // ---- HMR: background mesh ----------------------------------------------------
    Module_Parameter_Lists& tHmr = aDeck.hmr();
    tHmr.set( "number_of_elements_per_dimension", 45, 15 );
    tHmr.set( "domain_dimensions", 3.0, 1.0 );
    tHmr.set( "refinement_buffer", 3 );
    tHmr.set( "staircase_buffer", 3 );
    tHmr.set( "pattern_initial_refinement", tInitialRef );

    tHmr( HMR::LAGRANGE_MESHES ).add_parameter_list();
    tHmr.set( "order", 1 );

    tHmr( HMR::BSPLINE_MESHES ).add_parameter_list();
    tHmr.set( "orders", 1 );

    tHmr( HMR::BSPLINE_MESHES ).add_parameter_list();
    tHmr.set( "orders", 1 );

    // ---- XTK: conformal cut, ghost on both material phases -----------------------
    Module_Parameter_Lists& tXtk = aDeck.xtk();
    tXtk.set( "decompose", true );
    tXtk.set( "decomposition_type", "conformal" );
    tXtk.set( "enrich_mesh_indices", "0,1" );
    tXtk.set( "ghost_stab", tUseGhost );
    tXtk.set( "multigrid", false );
    tXtk.set( "verbose", false );
    tXtk.set( "print_enriched_ig_mesh", false );
    tXtk.set( "exodus_output_XTK_ig_mesh", true );
    tXtk.set( "high_to_low_dbl_side_sets", true );

    // ---- GEN: superellipse frame + hole-seeded design field ----------------------
    // (IQI_types is written automatically from the objective/constraint expressions)
    Module_Parameter_Lists& tGen = aDeck.gen();
    tGen.set( "output_mesh_file", "GEN_LevelSet_Boxbeam_NewIO.exo" );

    Matrix< DDUMat > tPhaseMap( 4, 1, 0 );
    tPhaseMap( 1 ) = 1;
    tPhaseMap( 2 ) = 2;
    tPhaseMap( 3 ) = 2;
    tGen.set( "phase_table", moris::ios::stringify( tPhaseMap ) );
    tGen.set( "print_phase_table", true );

    // outer frame
    tGen( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
    tGen.set( "center_x", 1.5 );
    tGen.set( "center_y", 0.5 );
    tGen.set( "semidiameter_x", 1.5 - tWallThickness );
    tGen.set( "semidiameter_y", 0.5 - tWallThickness );
    tGen.set( "exponent", 24.0 );

    // hole seeding array; its B-spline coefficients are the ADVs
    tGen( GEN::GEOMETRIES ).add_parameter_list( prm::create_field_array_parameter_list( gen::Field_Type::SUPERELLIPSE ) );
    tGen.set( "semidiameter_x", tHoleDiameter );
    tGen.set( "semidiameter_y", tHoleDiameter );
    tGen.set( "exponent", 8.0 );
    tGen.set( "lower_bound_x", 0.12 );
    tGen.set( "upper_bound_x", 2.88 );
    tGen.set( "lower_bound_y", 0.12 );
    tGen.set( "upper_bound_y", 0.88 );
    tGen.set( "number_of_fields_x", tNumHoleX );
    tGen.set( "number_of_fields_y", tNumHoleY );
    tGen.set( "discretization_mesh_index", 0 );
    tGen.set( "discretization_lower_bound", -tBsplineLimit );
    tGen.set( "discretization_upper_bound", tBsplineLimit );

    // ---- FEM: one linear-elastic preset per material phase ------------------------
    Module_Parameter_Lists& tFem = aDeck.fem();
    tFem.hack_for_legacy_fem();

    fem::presets::Linear_Elastic_Config tFrameConfig;
    tFrameConfig.mPrefix           = "Frame_";
    tFrameConfig.mBulkSets         = tFrameSets;
    tFrameConfig.mDirichletSets    = tFrameSupportSSets;
    tFrameConfig.mNeumannSets      = tFrameLoadSSets;
    tFrameConfig.mTraction         = "1.0";
    tFrameConfig.mTractionFunction = "Func_Neumann_U";
    tFrameConfig.mGhostSets        = tUseGhost ? tFrame.ghost() : "";
    tFrameConfig.mBedding          = "1.0e-6";
    fem::presets::Linear_Elastic_Names tFrameNames = fem::presets::linear_elastic( tFem, tFrameConfig );

    fem::presets::Linear_Elastic_Config tInteriorConfig;
    tInteriorConfig.mPrefix    = "Interior_";
    tInteriorConfig.mBulkSets  = tInteriorSets;
    tInteriorConfig.mGhostSets = tUseGhost ? tInterior.ghost() : "";
    tInteriorConfig.mBedding   = "1.0e-6";
    fem::presets::Linear_Elastic_Names tInteriorNames = fem::presets::linear_elastic( tFem, tInteriorConfig );

    // frame-interior tie: Nitsche interface (raw parameter-list escape hatch)
    tFem( FEM::STABILIZATION ).add_parameter_list();
    tFem.set( "stabilization_name", "SPNitscheFrameInteriorInterface" );
    tFem.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
    tFem.set( "function_parameters", "100.0" );
    tFem.set( "leader_properties", tFrameNames.mYoungsProp + ",Material" );
    tFem.set( "follower_properties", tInteriorNames.mYoungsProp + ",Material" );

    tFem( FEM::IWG ).add_parameter_list();
    tFem.set( "IWG_name", "IWGFrameInteriorInterface" );
    tFem.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE );
    tFem.set( "dof_residual", deck::Dofs::Displacement2D );
    tFem.set( "leader_dof_dependencies", deck::Dofs::Displacement2D );
    tFem.set( "follower_dof_dependencies", deck::Dofs::Displacement2D );
    tFem.set( "leader_constitutive_models", tFrameNames.mCM + ",ElastLinIso" );
    tFem.set( "follower_constitutive_models", tInteriorNames.mCM + ",ElastLinIso" );
    tFem.set( "stabilization_parameters", "SPNitscheFrameInteriorInterface,NitscheInterface" );
    tFem.set( "mesh_set_names", tFrameInteriorDSets );

    // material-void interface perimeter (regularizer criterion)
    tFem( FEM::IQI ).add_parameter_list();
    tFem.set( "IQI_name", "IQIPerimeter_InterfaceVoid" );
    tFem.set( "IQI_type", fem::IQI_Type::VOLUME );
    tFem.set( "leader_dof_dependencies", deck::Dofs::Displacement2D );
    tFem.set( "mesh_set_names", tInterfaceVoidSSets );

    // displacement IQIs for the output mesh
    tFem( FEM::IQI ).add_parameter_list();
    tFem.set( "IQI_name", "IQIBulkUX" );
    tFem.set( "IQI_type", fem::IQI_Type::DOF );
    tFem.set( "dof_quantity", deck::Dofs::Displacement2D );
    tFem.set( "leader_dof_dependencies", deck::Dofs::Displacement2D );
    tFem.set( "vectorial_field_index", 0 );
    tFem.set( "mesh_set_names", tTotalDomainSets );

    tFem( FEM::IQI ).add_parameter_list();
    tFem.set( "IQI_name", "IQIBulkUY" );
    tFem.set( "IQI_type", fem::IQI_Type::DOF );
    tFem.set( "dof_quantity", deck::Dofs::Displacement2D );
    tFem.set( "leader_dof_dependencies", deck::Dofs::Displacement2D );
    tFem.set( "vectorial_field_index", 1 );
    tFem.set( "mesh_set_names", tTotalDomainSets );

    tFem( FEM::COMPUTATION );

    // ---- OPT: GCMMA on the criteria expressions -----------------------------------
    // Replaces the legacy 7-callback set: values evaluated on the expression trees,
    // criteria gradients by reverse-mode differentiation, GEN IQI_types in
    // first-appearance order, constraint type from the comparison.
    Module_Parameter_Lists& tOpt = aDeck.opt();
    tOpt.set( "is_optimization_problem", true );
    tOpt.set( "problem", "user_defined" );
    tOpt.set( "restart_file", "" );

    tOpt( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::GCMMA );
    tOpt.set( "step_size", step_size() );
    tOpt.set( "penalty", tMMAPenalty );
    tOpt.set( "max_its", max_iterations() );

    aDeck.objective(
            deck::criterion( tFrameNames.mStrainEnergyIQI ) / tInitialStrainEnergy
            + deck::criterion( tInteriorNames.mStrainEnergyIQI ) / tInitialStrainEnergy
            + tPerimeterPenalty * deck::criterion( "IQIPerimeter_InterfaceVoid" ) / tInitialPerimeter );

    aDeck.constraint(
            deck::criterion( tInteriorNames.mVolumeIQI ) / tMaxMass - 1.0 <= 0.0 );

    // ---- SOL: direct solve, single static step ------------------------------------
    Module_Parameter_Lists& tSol = aDeck.sol();
    tSol( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
    // KLU explicitly: the Amesos default (Pardiso) is not built in every Trilinos
    tSol.set( "Solver_Type", "Amesos_Klu" );

    tSol( SOL::LINEAR_SOLVERS ).add_parameter_list();

    tSol( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
    tSol.set( "NLA_combined_res_jac_assembly", true );
    tSol.set( "NLA_rel_res_norm_drop", 1.00 );
    tSol.set( "NLA_relaxation_parameter", 1.00 );
    tSol.set( "NLA_max_iter", 1 );

    tSol( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
    tSol.set( "NLA_DofTypes", deck::Dofs::Displacement2D );

    tSol( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
    tSol.set( "TSA_Num_Time_Steps", 1 );
    tSol.set( "TSA_Time_Frame", 1.0 );

    tSol( SOL::TIME_SOLVERS ).add_parameter_list();
    tSol.set( "TSA_DofTypes", deck::Dofs::Displacement2D );
    tSol.set( "TSA_Output_Indices", "0" );
    tSol.set( "TSA_Output_Criteria", "Output_Criterion" );    // built-in always-true criterion

    tSol( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );

    // ---- MSI: displacements on B-spline mesh 1 ------------------------------------
    Module_Parameter_Lists& tMsi = aDeck.msi();
    tMsi.set( "UX", 1 );
    tMsi.set( "UY", 1 );

    // ---- VIS: displacement output on the material sets ----------------------------
    Module_Parameter_Lists& tVis = aDeck.vis();
    tVis.set( "File_Name", std::pair< std::string, std::string >( "./", "LevelSet_Boxbeam_NewIO.exo" ) );
    tVis.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
    tVis.set( "Set_Names", tTotalDomainSets );
    tVis.set( "Field_Names", "UX,UY" );
    tVis.set( "Field_Type", "NODAL,NODAL" );
    tVis.set( "IQI_Names", "IQIBulkUX,IQIBulkUY" );
    tVis.set( "Save_Frequency", 1 );
    tVis.set( "Time_Offset", 10.0 );

    // ---- MORISGENERAL: defaults ---------------------------------------------------
    aDeck.morisgeneral();
}
