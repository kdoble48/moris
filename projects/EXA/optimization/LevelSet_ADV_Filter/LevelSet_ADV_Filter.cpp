/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * LevelSet_ADV_Filter.cpp
 *
 * B-spline level-set designs with and without the Kreissl-Maute linear filter of the level-set B-spline coefficients
 * (GEN level-set parameter "adv_filter_radius"). With the filter the design field is F psi, with
 * F_ij ~ max( 0, r - |x_i - x_j| ) (rows normalized) on the Greville points of the coefficients; the optimizer keeps
 * psi and the sensitivities are F^T times the coefficient sensitivities.
 *
 * Common: 2 x 1 domain on [0, 2] x [0, 1], quadratic Lagrange and B-splines, linear elasticity (E = 1, nu = 0.3),
 * clamped at x = 0 (Nitsche), material where the level set is positive, ghost stabilization, KLU direct solver.
 *
 * Cases 0 and 1: sensitivity sweep (OPT SWEEP) at the initial design, analytic vs finite differences for every ADV
 *   12 x 6 cells (h = 1/6); one circular hole, center (1.0, 0.5), radius 0.3 (112 ADVs); traction (0, -1) on x = 2
 *   criteria: strain energy (objective) and volume (constraint)
 *   case 0: no filter; case 1: filter radius 1.5 h
 *
 * Cases 2 and 3: cantilever optimization, min strain energy s.t. volume <= 50 % of the domain
 *   40 x 20 cells (h = 1/20); seed: staggered grid of 13 circular holes of radius 0.12 (924 ADVs, bounds +- 3 h)
 *   load: traction (0, -10) on x = 2, y in [0.45, 0.55]
 *   cut on the quadratic interpolant of the level set (use_cell_order_intersection)
 *   native MMA (GCMMA inner loop), 100 iterations, no level-set reinitialization
 *   case 2: no filter; case 3: filter radius 1.5 h
 */

#include <string>
#include <iostream>
#include <algorithm>
#include <cmath>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "parameters.hpp"

// case index set by the test (EXA_Globals)
extern moris::uint gCaseIndex;

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    // cases 0, 1: sweep on 12 x 6 cells; cases 2, 3: optimization on 40 x 20 cells
    bool
    is_optimization()
    {
        return gCaseIndex >= 2;
    }

    bool
    use_filter()
    {
        return gCaseIndex == 1 or gCaseIndex == 3;
    }

    real
    cell_size()
    {
        return is_optimization() ? 1.0 / 20.0 : 1.0 / 6.0;
    }

    // optimization cases: volume limit (50 % of the domain) and the strain energy of the initial design (normalization)
    real tMaxVolume  = 0.5 * 2.0;
    real tInitialSE  = -1.0;
    real tSeedRadius = 0.12;
    real tLoadY0     = 0.45;
    real tLoadY1     = 0.55;

    // material = outside the hole (level set > 0, phase 1)
    std::string tMaterialSets = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tClampSets    = "SideSet_4_n_p1,SideSet_4_c_p1";
    std::string tLoadSets     = "SideSet_2_n_p1,SideSet_2_c_p1";
    std::string tGhostSets    = "ghost_p1";

    //--------------------------------------------------------------------------------------------------------------

    void
    Func_Const(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    // traction (0, t) on x = 2: everywhere (sweep cases) or on y in [0.45, 0.55] (optimization cases)
    void
    Func_Traction(
            moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        real tY       = aFIManager->get_IG_geometry_interpolator()->valx()( 1 );
        bool tOnPatch = not is_optimization() or ( tY >= tLoadY0 - 1e-12 and tY <= tLoadY1 + 1e-12 );
        aPropMatrix   = { { 0.0 }, { tOnPatch ? aParameters( 0 )( 0 ) : 0.0 } };
    }

    // optimization seed: staggered grid of circular holes, phi = min_k ( |x - c_k| - r ), clipped to the ADV bounds
    real
    Func_Hole_Grid(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        real tPhi = 1.0e10;
        for ( uint iColumn = 0; iColumn < 5; iColumn++ )
        {
            real           tX    = 0.2 + 0.4 * iColumn;
            Vector< real > tRows = iColumn % 2 == 0 ? Vector< real >{ 0.25, 0.75 } : Vector< real >{ 0.0, 0.5, 1.0 };
            for ( real tY : tRows )
            {
                tPhi = std::min( tPhi, std::hypot( aCoordinates( 0 ) - tX, aCoordinates( 1 ) - tY ) - tSeedRadius );
            }
        }
        real tBound = 3.0 * cell_size();
        return std::max( -tBound, std::min( tBound, tPhi ) );
    }

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------
    // OPT callbacks: criterion 0 = strain energy (objective), criterion 1 = volume (constraint)

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );
        return tConstraintTypes;
    }

    // sweep cases: the criteria as they are; optimization cases: SE / SE(initial design) and V / V_max - 1
    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        if ( is_optimization() and tInitialSE < 0.0 )
        {
            tInitialSE = aCriteria( 0 );
        }
        Matrix< DDRMat > tObjectives = { { is_optimization() ? aCriteria( 0 ) / tInitialSE : aCriteria( 0 ) } };
        return tObjectives;
    }

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints = { { is_optimization() ? aCriteria( 1 ) / tMaxVolume - 1.0 : aCriteria( 1 ) } };
        return tConstraints;
    }

    Matrix< DDRMat >
    compute_dobjective_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDADV( 1, aADVs.size(), 0.0 );
        return tDObjectiveDADV;
    }

    Matrix< DDRMat >
    compute_dobjective_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDObjectiveDCriteria = { { is_optimization() ? 1.0 / tInitialSE : 1.0, 0.0 } };
        return tDObjectiveDCriteria;
    }

    Matrix< DDRMat >
    compute_dconstraint_dadv( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDADV( 1, aADVs.size(), 0.0 );
        return tDConstraintDADV;
    }

    Matrix< DDRMat >
    compute_dconstraint_dcriteria( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tDConstraintDCriteria = { { 0.0, is_optimization() ? 1.0 / tMaxVolume : 1.0 } };
        return tDConstraintDCriteria;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        if ( is_optimization() )
        {
            aParameterLists.set( "number_of_elements_per_dimension", 40, 20 );
        }
        else
        {
            aParameterLists.set( "number_of_elements_per_dimension", 12, 6 );
        }
        aParameterLists.set( "domain_dimensions", 2.0, 1.0 );
        aParameterLists.set( "domain_offset", 0.0, 0.0 );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "2" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "2" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", 3 );
        aParameterLists.set( "staircase_buffer", 3 );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "print_enriched_ig_mesh", false );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", false );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIBulkStrainEnergy", "IQIBulkVolume" );
        aParameterLists.set( "output_mesh_file", "gen_levelset_adv_filter.exo" );

        if ( is_optimization() )
        {
            // seed: staggered grid of holes (material where phi > 0), bounds +- 3 h
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Func_Hole_Grid" );
            aParameterLists.set( "discretization_lower_bound", -3.0 * cell_size() );
            aParameterLists.set( "discretization_upper_bound", 3.0 * cell_size() );

            // cut on the cell's own-order (quadratic) interpolant of the level set, so the analysed boundary follows the
            // design field itself rather than its corner values
            aParameterLists.set( "use_cell_order_intersection", true );
        }
        else
        {
            // hole: phi = distance to the center - radius (material where phi > 0)
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::CIRCLE ) );
            aParameterLists.set( "center_x", 1.0 );
            aParameterLists.set( "center_y", 0.5 );
            aParameterLists.set( "radius", 0.3 );
            aParameterLists.set( "discretization_lower_bound", -2.0 );
            aParameterLists.set( "discretization_upper_bound", 2.0 );
        }
        aParameterLists.set( "discretization_mesh_index", 0 );

        // cases 1 and 3: Kreissl-Maute filter of the B-spline coefficients, radius 1.5 cells
        if ( use_filter() )
        {
            aParameterLists.set( "adv_filter_radius", 1.5 * cell_size() );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoungs" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.3" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDirichletU" );
        aParameterLists.set( "function_parameters", "0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropTraction" );
        aParameterLists.set( "function_parameters", is_optimization() ? "-10.0" : "-1.0" );
        aParameterLists.set( "value_function", "Func_Traction" );

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPGhost" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
        aParameterLists.set( "function_parameters", "0.01" );
        aParameterLists.set( "leader_properties", "PropYoungs,Material" );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulkU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichletU" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropDirichletU,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tClampSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGTraction" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_NEUMANN );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_properties", "PropTraction,Traction" );
        aParameterLists.set( "mesh_set_names", tLoadSets );

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGGhost" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
        aParameterLists.set( "dof_residual", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists.set( "stabilization_parameters", "SPGhost,GhostSP" );
        aParameterLists.set( "mesh_set_names", tGhostSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkStrainEnergy" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,Elast" );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkVolume" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
        aParameterLists.set( "leader_properties", "PropDensity,Density" );
        aParameterLists.set( "mesh_set_names", tMaterialSets );

        aParameterLists( FEM::COMPUTATION );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // direct sparse solve with KLU (part of every Trilinos Amesos build)
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aParameterLists.set( "Solver_Type", "Amesos_Klu" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );
        aParameterLists.set( "NLA_max_iter", 1 );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Num_Time_Steps", 1 );
        aParameterLists.set( "TSA_Time_Frame", 1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_DofTypes", "UX,UY" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "levelset_adv_filter.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tMaterialSets );
        aParameterLists.set( "Field_Names", "UX,StrainEnergy,Volume" );
        aParameterLists.set( "Field_Type", "NODAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkUX,IQIBulkStrainEnergy,IQIBulkVolume" );
        aParameterLists.set( "Save_Frequency", 1 );
        if ( is_optimization() )
        {
            // one file per evaluation (levelset_adv_filter.exo.e-s.NNNN): the design history
            aParameterLists.set( "Time_Offset", 10.0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", true );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", "LevelSet_ADV_Filter.so" );

        if ( is_optimization() )
        {
            // native MMA with the GCMMA conservativeness inner loop; no level-set reinitialization
            aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::MMA );
            aParameterLists.set( "step_size", 0.05 );
            aParameterLists.set( "penalty", 5.0 );
            aParameterLists.set( "max_its", 100 );
            aParameterLists.set( "max_inner_its", 5 );
            aParameterLists.set( "norm_drop", 1e-9 );
            return;
        }

        // sensitivity sweep at the initial design: analytic and finite differences (forward, backward, central) for every ADV
        // ("all": a single finite-difference mode skips the analytic gradients)
        aParameterLists( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
        aParameterLists.set( "hdf5_path", "levelset_adv_filter.hdf5" );
        aParameterLists.set( "evaluate_objective_gradients", true );
        aParameterLists.set( "evaluate_constraint_gradients", true );
        aParameterLists.set( "num_evaluations_per_adv", "1" );
        aParameterLists.set( "include_bounds", false );
        aParameterLists.set( "finite_difference_type", "all" );
        aParameterLists.set( "finite_difference_epsilons", "1E-6" );
    }

    //--------------------------------------------------------------------------------------------------------------

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
