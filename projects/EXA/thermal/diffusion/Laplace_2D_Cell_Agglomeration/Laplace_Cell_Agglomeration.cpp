/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Laplace_Cell_Agglomeration.cpp
 *
 * Cell agglomeration on a sliver cut: when to turn "activate_cell_agglomeration" on, and what it does.
 *
 * Problem
 *   Laplace equation on the solid x < x_c of the unit square, background mesh 20 x 20 with element size h = 1/20.
 *   The immersed face sits at x_c = ( 10 + alpha ) * h with alpha = 1e-3, so every cell of column 10 is cut and keeps
 *   only the fraction alpha = 0.1 % of its area as material (a "sliver").
 *
 *   Exact solution (harmonic): T = sin( pi y ) * cosh( pi x ) / cosh( pi x_c )
 *     left face   x = 0     : zero flux (natural)
 *     bottom/top  y = 0, 1  : T = 0, symmetric Nitsche
 *     immersed    x = x_c   : T = sin( pi y ), symmetric Nitsche
 *
 * Why cell agglomeration
 *   The basis functions that live only on the sliver have almost no material under them. Their diagonal entries
 *   scale with alpha, while the symmetric Nitsche terms on the immersed face do not. With the Nitsche penalty of this
 *   deck ( 10 k / h ) and no stabilization, the system matrix is indefinite: 21 negative eigenvalues for linear and
 *   22 for quadratic B-splines (252 and 286 unknowns). The discrete problem is then no longer coercive, and the
 *   condition number of the quadratic system is about 3e10.
 *
 *   Cell agglomeration ( XTK: "activate_cell_agglomeration" = true ) attaches each cut cell to an uncut neighbour, its
 *   root cell. Every basis function that is supported only on cut cells is removed from the list of unknowns and
 *   replaced by the extension of the root cell's polynomial. Here each of the 20 sliver cells is attached to its left
 *   neighbour, the unknowns drop from 252 to 231 (linear) and from 286 to 264 (quadratic), the matrix is positive
 *   definite, and its condition number is about 370 without any ghost penalty.
 *
 *   The face-oriented ghost penalty ( "ghost_stab" ) is the alternative. It keeps all unknowns and adds a penalty
 *   term with a user-chosen weight; cell agglomeration has no weight to tune. This deck leaves the ghost penalty off
 *   so that the result is due to the agglomeration alone.
 *
 *   Note: the cell-agglomeration path agglomerates every cut cell. The XTK parameter "volume_fraction" is used only
 *   by the basis agglomeration ( "activate_basis_agglomeration" ).
 *
 * Variants (global gTestIndex, set by the test case)
 *   0: linear Lagrange elements on linear B-splines
 *   1: quadratic Lagrange elements on quadratic B-splines
 *
 * Output
 *   <name>.exo                : temperature, exact temperature, integrals of the squared L2 and H1 errors, solid volume
 *   xtk_temp_ip.exo           : interpolation mesh with the subphase-group index and the root-cell index of every cell
 *   <name>.0.0.{jac,res}.dat  : system matrix and residual (the length of the residual is the number of unknowns)
 *   with <name> = Laplace_Cell_Agglomeration_Linear or Laplace_Cell_Agglomeration_Quadratic
 */

#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <cmath>

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "cl_Bitset.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

// variant: 0 linear on linear, 1 quadratic on quadratic
extern uint gTestIndex;

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // helpers

    std::string
    moris_to_string( real tValue )
    {
        std::ostringstream streamObj;
        streamObj << std::scientific;
        streamObj << std::setprecision( 15 );
        streamObj << tValue;
        return streamObj.str();
    }

    // Lagrange and B-spline order of the variant; evaluated when the parameter lists are built
    uint
    get_interpolation_order()
    {
        return gTestIndex == 0 ? 1 : 2;
    }

    // base name of the output files of the variant
    std::string
    get_output_name()
    {
        return gTestIndex == 0 ? "Laplace_Cell_Agglomeration_Linear" : "Laplace_Cell_Agglomeration_Quadratic";
    }

    /* ------------------------------------------------------------------------ */
    // stabilization of the cut cells

    // cell agglomeration: each cut cell is attached to an uncut root cell (see the header comment)
    bool tUseCellAgglomeration = true;

    // face-oriented ghost penalty: off, the agglomeration alone stabilizes the sliver cells
    bool tUseGhost = false;

    // material fraction of the cut cells
    real tAlpha = 1.0e-3;

    // background elements per direction
    uint tNumElem = 20;

    /* ------------------------------------------------------------------------ */
    // geometry

    real tDimX = 1.0;
    real tDimY = 1.0;

    real tElementSize = tDimX / (real)tNumElem;

    // position of the immersed face
    real tPlanePos = ( (real)( tNumElem / 2 ) + tAlpha ) * tElementSize;

    /* ------------------------------------------------------------------------ */
    // material and stabilization parameters

    std::string tDensity      = moris_to_string( 1.0 );
    std::string tCapacity     = moris_to_string( 1.0 );
    std::string tConductivity = moris_to_string( 1.0 );

    // Nitsche penalty: gamma * k / h with gamma = 10
    std::string tNitschePenThermal = moris_to_string( 10.0 );

    /* ------------------------------------------------------------------------ */
    // solver parameters

    int  tNLA_max_iter             = 1;
    real tNLA_rel_res_norm_drop    = 1e-12;
    real tNLA_relaxation_parameter = 1.0;

    /* ------------------------------------------------------------------------ */
    // file names

    std::string tSoFile = "Laplace_Cell_Agglomeration.so";

    /* ------------------------------------------------------------------------ */
    // mesh sets

    std::string tSolidPhase = "HMR_dummy_n_p0,HMR_dummy_c_p0";

    /* ------------------------------------------------------------------------ */
    // level set: negative in the solid

    real
    Func_Plane(
            const Matrix< DDRMat >& aCoordinates,
            const Vector< real >&   aGeometryParameters )
    {
        return aCoordinates( 0 ) - tPlanePos;
    }

    /* ------------------------------------------------------------------------ */
    // exact solution and its gradient

    void
    Func_Exact_Temp(
            Matrix< DDRMat >&                aPropMatrix,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        aPropMatrix.set_size( 1, 1 );
        aPropMatrix( 0, 0 ) = std::sin( M_PI * tY ) * std::cosh( M_PI * tX ) / std::cosh( M_PI * tPlanePos );
    }

    void
    Func_Exact_Grad(
            Matrix< DDRMat >&                aPropMatrix,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        real tX = aFIManager->get_IP_geometry_interpolator()->valx()( 0 );
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        aPropMatrix.set_size( 2, 1 );
        aPropMatrix( 0, 0 ) = M_PI * std::sin( M_PI * tY ) * std::sinh( M_PI * tX ) / std::cosh( M_PI * tPlanePos );
        aPropMatrix( 1, 0 ) = M_PI * std::cos( M_PI * tY ) * std::cosh( M_PI * tX ) / std::cosh( M_PI * tPlanePos );
    }

    // temperature prescribed on the immersed face
    void
    Func_Immersed_Temp(
            Matrix< DDRMat >&                aPropMatrix,
            Vector< Matrix< DDRMat > >&      aParameters,
            fem::Field_Interpolator_Manager* aFIManager )
    {
        real tY = aFIManager->get_IP_geometry_interpolator()->valx()( 1 );

        aPropMatrix.set_size( 1, 1 );
        aPropMatrix( 0, 0 ) = std::sin( M_PI * tY );
    }

    /* ------------------------------------------------------------------------ */
    // phase assignment: 0 solid (level set negative), 1 void

    uint
    get_phase_index( const Bitset< 1 >& aGeometrySigns )
    {
        if ( aGeometrySigns.test( 0 ) )
        {
            return 1;
        }
        return 0;
    }

    /* ------------------------------------------------------------------------ */
    // optimization interface (unused: forward analysis only)

    Matrix< DDSMat >
    get_constraint_types()
    {
        Matrix< DDSMat > tConstraintTypes( 1, 1, 1 );
        return tConstraintTypes;
    }

    Matrix< DDRMat >
    compute_objectives( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tObjectives = { { aCriteria( 0 ) } };
        return tObjectives;
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
        Matrix< DDRMat > tDObjectiveDCriteria( 1, aCriteria.size(), 0.0 );
        tDObjectiveDCriteria( 0 ) = 1.0;
        return tDObjectiveDCriteria;
    }

    Matrix< DDRMat >
    compute_constraints( const Vector< real >& aADVs, const Vector< real >& aCriteria )
    {
        Matrix< DDRMat > tConstraints = { { aCriteria( 0 ) } };
        return tConstraints;
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
        Matrix< DDRMat > tDConstraintDCriteria( 1, aCriteria.size(), 0.0 );
        tDConstraintDCriteria( 0, 0 ) = 1.0;
        return tDConstraintDCriteria;
    }

    bool
    Output_Criterion( tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
        aParameterLists.set( "problem", "user_defined" );
        aParameterLists.set( "library", tSoFile );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "number_of_elements_per_dimension", tNumElem, tNumElem );
        aParameterLists.set( "domain_dimensions", tDimX, tDimY );
        aParameterLists.set( "domain_offset", 0.0, 0.0 );

        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", std::to_string( get_interpolation_order() ) );
        aParameterLists.set( "lagrange_pattern", "0" );

        aParameterLists.set( "bspline_orders", std::to_string( get_interpolation_order() ) );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "refinement_buffer", 1 );
        aParameterLists.set( "staircase_buffer", 1 );

        aParameterLists.set( "pattern_initial_refinement", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "use_SPG_based_enrichment", true );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ip_mesh", true );

        // attach every cut cell to an uncut root cell and replace the basis functions that are
        // supported only on cut cells by the extension of the root cell's polynomial
        aParameterLists.set( "activate_cell_agglomeration", tUseCellAgglomeration );

        // write the subphase-group index and the root-cell index of every cell to xtk_temp_ip.exo
        aParameterLists.set( "visualize_cell_association", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "IQI_types", "IQIBulkL2Error" );

        aParameterLists.set( "number_of_phases", 2 );
        aParameterLists.set( "phase_function_name", "get_phase_index" );
        aParameterLists.set( "output_mesh_file", get_output_name() + "_GEN.exo" );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "field_function_name", "Func_Plane" );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );
        aParameterLists.set( "use_multilinear_interpolation", true );
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // phases
        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "phase_indices", "0" );

        aParameterLists( FEM::PHASES ).add_parameter_list();
        aParameterLists.set( "phase_name", "PhaseVoid" );
        aParameterLists.set( "phase_indices", "1" );

        //------------------------------------------------------------------------------
        // properties

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropDensity" );
        aParameterLists.set( "function_parameters", tDensity );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropCapacity" );
        aParameterLists.set( "function_parameters", tCapacity );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropConductivity" );
        aParameterLists.set( "function_parameters", tConductivity );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropWallTemp" );
        aParameterLists.set( "function_parameters", "0.0" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropImmersedTemp" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Immersed_Temp" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropExactTemp" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Exact_Temp" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropExactGrad" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Exact_Grad" );

        //------------------------------------------------------------------------------
        // constitutive model

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMDiffusion" );
        aParameterLists.set( "phase_name", "PhaseSolid" );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::DIFF_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "TEMP", "Temperature" ) );
        aParameterLists.set( "properties",
                "PropConductivity,Conductivity;"
                "PropDensity     ,Density;"
                "PropCapacity    ,HeatCapacity" );

        //------------------------------------------------------------------------------
        // stabilization parameters

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", tNitschePenThermal );
        aParameterLists.set( "leader_properties", "PropConductivity,Material" );

        //------------------------------------------------------------------------------
        // IWGs

        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDiffusionBulk" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_BULK );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );

        // T = 0 on the bottom and top faces of the background mesh
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGWallTemp" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "side_ordinals", "1,3" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropWallTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        // T = sin( pi y ) on the immersed face
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGImmersedTemp" );
        aParameterLists.set( "IWG_bulk_type", fem::Element_Type::SIDESET );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "neighbor_phases", "PhaseVoid" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::SPATIALDIFF_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "TEMP" );
        aParameterLists.set( "leader_properties", "PropImmersedTemp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMDiffusion,Diffusion" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );

        //------------------------------------------------------------------------------
        // IQIs

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkTEMP" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "vectorial_field_index", 0 );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIExactTemp" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::PROPERTY );
        aParameterLists.set( "leader_properties", "PropExactTemp,Property" );
        aParameterLists.set( "vectorial_field_index", 0 );

        // integral of ( T_h - T )^2
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkL2Error" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_properties", "PropExactTemp,L2Check" );

        // integral of | grad T_h - grad T |^2
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIBulkH1Error" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::H1_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "TEMP" );
        aParameterLists.set( "leader_properties", "PropExactGrad,H1Check" );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQISolidVolume" );
        aParameterLists.set( "leader_phase_name", "PhaseSolid" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );

        //------------------------------------------------------------------------------
        aParameterLists( FEM::COMPUTATION );
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // GMRES with ILU(0), both from Aztec; no external preconditioner
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AZTEC_IMPL );
        aParameterLists.set( "AZ_solver", AZ_gmres );
        aParameterLists.set( "AZ_max_iter", 1000 );
        aParameterLists.set( "rel_residual", 1e-10 );
        aParameterLists.set( "AZ_conv", AZ_r0 );
        aParameterLists.set( "AZ_diagnostics", AZ_all );
        aParameterLists.set( "AZ_output", AZ_all );
        aParameterLists.set( "AZ_kspace", 1000 );
        aParameterLists.set( "AZ_orthog", AZ_modified );
        aParameterLists.set( "AZ_precond", AZ_dom_decomp );
        aParameterLists.set( "AZ_subdomain_solve", AZ_ilu );
        aParameterLists.set( "AZ_graph_fill", 0 );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        // print the condition number of the system matrix to the log
        aParameterLists.set( "DLA_operator_condition_number_with_moris", "dense" );

        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Solver_Implementation", NLA::NonlinearSolverType::NEWTON_SOLVER );
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_rel_res_norm_drop", tNLA_rel_res_norm_drop );
        aParameterLists.set( "NLA_relaxation_parameter", tNLA_relaxation_parameter );
        aParameterLists.set( "NLA_max_iter", tNLA_max_iter );

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

        // writes <name>.0.0.jac.dat and <name>.0.0.res.dat
        aParameterLists( SOL::SOLVER_WAREHOUSE ).set( "SOL_save_operator_to_matlab", get_output_name() );
    }

    /* ------------------------------------------------------------------------ */

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "TEMP", 0 );
    }

    /* ------------------------------------------------------------------------ */

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", get_output_name() + ".exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tSolidPhase );
        aParameterLists.set( "Field_Names", "TEMP,TEMP_EXACT,L2ErrSq,H1ErrSq,Volume" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,GLOBAL,GLOBAL,GLOBAL" );
        aParameterLists.set( "IQI_Names", "IQIBulkTEMP,IQIExactTemp,IQIBulkL2Error,IQIBulkH1Error,IQISolidVolume" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    /* ------------------------------------------------------------------------ */

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
