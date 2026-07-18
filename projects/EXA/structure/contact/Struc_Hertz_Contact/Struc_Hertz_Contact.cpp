/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Struc_Hertz_Contact.cpp
 *
 * 2D frictionless Hertzian-contact benchmark (small-patch regime).
 *
 * A rigid-ish parabolic indenter (top body, curved lower boundary) is pressed a small
 * DISPLACEMENT-CONTROLLED penetration into a flat elastic body (bottom). The parabola
 * surface is  y_s(x) = tTopYShift + tFactor*(x - 0.5)^2  (apex DOWN at (0.5, ~0.5035),
 * local radius R = 1/(2*tFactor)). With a small penetration the contact half-width a is
 * small relative to R (a/R ~ 0.1), so the response sits in the Hertz small-patch regime.
 *
 * Loading is displacement-controlled (we MEASURE the load P and the contact pressure in
 * post, we do NOT prescribe P). The KEY output is the FULL 2D stress tensor per body on
 * the unzipped (FULL_DISCONTINUOUS) VIS mesh, so that the interface normal traction
 *     sigma_nn = nx^2 * sigma_xx + 2*nx*ny * sigma_xy + ny^2 * sigma_yy
 * can be recovered in post using the ANALYTIC parabola normal (see below), then compared
 * to the Hertz pressure distribution p(x) = p0*sqrt(1 - (x/a)^2).
 *
 * Analytic parabola normal (for post-processing sigma_nn):
 *   level set (top geometry) phi_top = (y - tTopYShift) - tFactor*(x - tTopXShift)^2
 *   surface  y_s(x)  = tTopYShift + tFactor*(x - tTopXShift)^2
 *   grad(phi_top)    = ( -2*tFactor*(x - tTopXShift), 1 )   (points UP, +phi side = top body)
 *   OUTWARD normal of the top (indenter) body, pointing DOWN into contact:
 *     n = ( 2*tFactor*(x - tTopXShift), -1 ) / sqrt( 1 + (2*tFactor*(x - tTopXShift))^2 )
 *   at the apex x = tTopXShift:  n = (0, -1)  ->  sigma_nn = sigma_yy.
 *
 * Contact formulation: neutral Nitsche, unbiased, theta = 0
 *   (STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED) -- the robust variant from the
 *   Struc_Contact_Patch_Test study -- with ghost stabilization ON and the combined
 *   residual/jacobian assembly path (the Nitsche IWG implements it).
 *
 * Adapted from ../Parabolic_Indenter/Parabolic_Indenter_Linear.cpp (geometry) and
 * ../Struc_Contact_Patch_Test/Struc_Contact_Patch_Test.cpp (phase-based contact wiring,
 * one-workflow-per-TEST_CASE harness, FULL_DISCONTINUOUS bulk-only VIS).
 *
 * Case matrix: gCaseIndex = load*3 + mesh  (0..8), a 3 (penetration) x 3 (mesh) sweep.
 *   load = gCaseIndex / 3 :  0 -> 0.0025, 1 -> 0.005, 2 -> 0.01   (tPenetrationLevels)
 *   mesh = gCaseIndex % 3 :  0 -> 40,     1 -> 80,    2 -> 160     (tMeshLevels, uniform NxN)
 */

#include "cl_Bitset.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_Matrix.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "parameters.hpp"
#include "fn_equal_to.hpp"
#include "fn_norm.hpp"
#include "linalg_typedefs.hpp"
#include "moris_typedefs.hpp"
#include <iostream>
#include <string>

#include "AztecOO.h"
//---------------------------------------------------------------

#define F2STR( func ) #func

using moris::fem::IWG_Type;

//---------------------------------------------------------------

// global variable for interpolation order
extern uint gInterpolationOrder;

// global variable for case index (gCaseIndex = load*3 + mesh, see header comment)
extern uint gCaseIndex;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
namespace moris
{
    /* ============================================================================================ */
    /*        Benchmark control constants -- RETUNE these after the first run measures a/R.          */
    /* ============================================================================================ */
    // Parabola curvature. local radius R = 1/(2*tFactor). tFactor = 1 -> R = 0.5.
    real tFactor = 1.0;

    // Penetration ladder (small-patch regime). geometric a ~ sqrt(2*R*delta), a/R ~ 2*sqrt(delta):
    //   0.0025 -> a/R ~ 0.10,  0.005 -> a/R ~ 0.14,  0.01 -> a/R ~ 0.20.  Retune after measuring a.
    real tPenetrationLevels[ 3 ] = { 0.0025, 0.005, 0.01 };

    // Uniform mesh ladder (NxN). Need >= 8-10 cells across the (small) contact half-width a;
    //   the finest level (160 -> h = 1/160 = 0.00625) gives ~8 cells across a ~ 0.05. If still
    //   under-resolved after measuring a, raise these or add HMR interface refinement.
    uint tMeshLevels[ 3 ] = { 40, 80, 160 };

    // Neutral-Nitsche interface penalty scaling (gamma).
    real tContactGamma = 100.0;

    /* -------------------------------- Case decode (gCaseIndex) ------------------------------------ */
    uint tLoadIndex = gCaseIndex / 3;    // 0..2 : penetration level
    uint tMeshIndex = gCaseIndex % 3;    // 0..2 : mesh level

    /* -------------------------------------- Geometry Setting -------------------------------------- */
    real tHeight = 1.0;
    real tWidth  = 1.0;

    uint tNumElems  = tMeshLevels[ tMeshIndex ];
    uint tNumElemsX = tNumElems;
    uint tNumElemsY = tNumElems;

    real tInitialGap           = 0.001;
    real tInitialGapMiddle     = 0.503 * tHeight;
    real tPenetration          = tPenetrationLevels[ tLoadIndex ];
    real tVerticalDisplacement = -tInitialGap - tPenetration;

    /* --------------------------------------- Top Geometry ----------------------------------------- */
    // HERTZ: curved parabolic indenter (apex down). Keep the source geometry (tTopIsPlane = false).
    bool tTopIsPlane = false;

    real tInterfaceTopAngle = 0 * M_PI / 180.0;
    real tTopYShift         = tInitialGapMiddle + tInitialGap / 2;
    real tTopXShift         = 0.5;
    real tTopParabolaFactor = tFactor;

    /* --------------------------------------- Bottom Geometry -------------------------------------- */
    // HERTZ: flat half-space body (flat horizontal interface).
    bool tBottomIsPlane = true;

    real tInterfaceBottomAngle = 0 * M_PI / 180.0;
    real tBottomYShift         = tInitialGapMiddle - tInitialGap / 2;
    real tBottomXShift         = 0.5;
    real tBottomParabolaFactor = -0.7;

    /* ------------------------------------- Material Parameters ------------------------------------ */
    real tTopEmod    = 1000.0;
    real tBottomEmod = 1000.0;
    real tTopPois    = 0.0;
    real tBottomPois = 0.0;

    /* --------------------------------------- Domain Setting --------------------------------------- */
    std::string tInterpolationOrder     = std::to_string( gInterpolationOrder );
    std::string tDomainTop              = "HMR_dummy_n_p1,HMR_dummy_c_p1";
    std::string tDomainBottom           = "HMR_dummy_n_p2,HMR_dummy_c_p2";
    std::string tContactInterfaceTop    = "ncss|iside_b0_1_b1_0|iside_b0_2_b1_0";
    std::string tContactInterfaceBottom = "ncss|iside_b0_2_b1_0|iside_b0_1_b1_0";
    std::string tContactInterface       = tContactInterfaceBottom + "," + tContactInterfaceTop;
    std::string tDomain                 = tDomainTop + "," + tDomainBottom;

    /* ------------------------------------------ File I/O ------------------------------------------ */
    std::string tOutputFileName =
            "Struc_Hertz_Contact_Case_" + std::to_string( gCaseIndex );

    /* ------------------------------------------- Solver ------------------------------------------- */
    int tMaxIterations = 200;

    real tRelResNormDrop     = 1e-8;
    real tRelaxation         = 1.0;
    real tFDPerturbationSize = 1e-8;

    fem::Perturbation_Type tFDPerturbationStrategy = fem::Perturbation_Type::ABSOLUTE;
    fem::FDScheme_Type     tFDScheme               = fem::FDScheme_Type::POINT_3_CENTRAL;

    int  tLoadControlSteps    = 10;
    real tLoadControlFactor   = 0.1;
    real tLoadControlRelRes   = 1e-2;
    real tLoadControlExponent = 1.0;

    real tRemapResidualChange     = -1e-2;
    int  tRemapLoadStepFrequency  = 1;
    int  tRemapIterationFrequency = 10;
    auto tRaytracingStrategy      = sol::SolverRaytracingStrategy::MixedNthLoadStepAndResidualChange;
    real tMaxNegativeRayLength    = -0.005;
    real tMaxPositiveRayLength    = 0.01;

    mtk::Integration_Order tNonconformalIntegrationOrder = mtk::Integration_Order::BAR_6;

    /* ------------------------------------------- Contact ------------------------------------------ */
    // Neutral Nitsche, unbiased (theta = 0) -- the robust variant from the patch-test study.
    real tContactStabilization = tContactGamma;    // Nitsche interface penalty scaling ("STAB")

    /* -------------------------------------- Control Variables ------------------------------------- */
    bool tOnlyGenerateMesh   = false;
    bool tUseGhost           = true;
    bool tTopHasDirichlet    = true;
    bool tBottomHasDirichlet = true;

    /* ---------------------------------------------------------------------------------------------- */
    /*                                         Field Functions                                        */
    /* ---------------------------------------------------------------------------------------------- */

    /* --------------------------------------- Phase Indexing --------------------------------------- */
    // geometry order: 0 = top parabola, 1 = bottom line.  bit(0) = top sign, bit(1) = bottom sign.
    //   above bottom line & above parabola -> 1 (top indenter body)
    //   above bottom line & below parabola -> 0 (void gap band between the bodies)
    //   below bottom line                  -> 2 (bottom half-space body)
    uint Phase_Index_Split( const Bitset< 2 > &aGeometrySigns )
    {
        if ( aGeometrySigns.test( 1 ) )
        {
            if ( aGeometrySigns.test( 0 ) )
            {
                return 1;    // upper (indenter) domain
            }
            return 0;    // interface (thin void band)
        }
        return 2;    // bottom (half-space)
    }

    /* ----------------------------------- Geometry Field Function ---------------------------------- */

    real tMinimumLevelSetValue = 1.0e-8;

    // Level set for the parabola splitting the domain into top (indenter) and bottom (half-space).
    //   phi = (y - tYShift) - tFactor*(x - tXShift)^2
    real Parabola( const Matrix< DDRMat > &aCoordinates,
            const Vector< real >          &aGeometryParameters )
    {
        real tXShift = aGeometryParameters( 0 );
        real tYShift = aGeometryParameters( 1 );
        real tFactorLocal = aGeometryParameters( 2 );

        real x = aCoordinates( 0 ) - tXShift;
        real y = aCoordinates( 1 ) - tYShift;

        // Compute Signed-Distance field
        real tVal = y - ( tFactorLocal * x * x );

        return ( std::abs( tVal ) < tMinimumLevelSetValue ) ? tMinimumLevelSetValue : tVal;
    }

    /* ----------------------------------- Property Field Function ---------------------------------- */

    void Func_Const( Matrix< DDRMat >       &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    void Func_Dirichlet( Matrix< DDRMat >   &aPropMatrix,
            Vector< Matrix< DDRMat > >      &aParameters,
            fem::Field_Interpolator_Manager *aFIManager )
    {
        real tLoadFactor = gLogger.get_action_data( "NonLinearAlgorithm", "Newton", "Solve", "LoadFactor" );
        aPropMatrix      = aParameters( 0 ) * tLoadFactor;
    }

    bool Output_Criterion( tsa::Time_Solver *aTimeSolver ) { return true; }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### OPT ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void OPTParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### HMR ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void HMRParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* --------------------------------------- Domain Settings ------------------------------------ */
        aParameterLists.set( "number_of_elements_per_dimension", tNumElemsX, tNumElemsY );
        aParameterLists.set( "domain_dimensions", tWidth, tHeight );

        /* ---------------------------------------- Lagrange Mesh ------------------------------------- */
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "lagrange_orders", tInterpolationOrder );
        aParameterLists.set( "lagrange_output_meshes", "0" );

        /* -------------------------------------- B-Spline Mesh --------------------------------------- */
        aParameterLists.set( "bspline_pattern", "0" );
        aParameterLists.set( "bspline_orders", tInterpolationOrder );
        aParameterLists.set( "lagrange_to_bspline", "0" );

        /* ----------------------------------------- Refinement --------------------------------------- */
        aParameterLists.set( "refinement_buffer", 0 );
        aParameterLists.set( "staircase_buffer", 0 );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### XTK ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void XTKParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", true );
        aParameterLists.set( "print_enriched_ig_mesh", true );
        aParameterLists.set( "exodus_output_XTK_ig_mesh", true );
        aParameterLists.set( "high_to_low_dbl_side_sets", true );
        aParameterLists.set( "only_generate_xtk_temp", tOnlyGenerateMesh );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### GEN ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void GENParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists.set( "number_of_phases", 3 );
        aParameterLists.set( "phase_function_name", F2STR( Phase_Index_Split ) );

        /* ----------------------------------------- Top (parabola) ----------------------------------- */
        if ( tTopIsPlane )
        {
            real tNormalX = std::cos( tInterfaceTopAngle + M_PI_2 );
            real tNormalY = std::sin( tInterfaceTopAngle + M_PI_2 );

            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            aParameterLists.set( "center_x", tTopXShift );
            aParameterLists.set( "center_y", tTopYShift );
            aParameterLists.set( "normal_x", tNormalX );
            aParameterLists.set( "normal_y", tNormalY );
            aParameterLists.set( "number_of_refinements", 0 );
            aParameterLists.set( "refinement_mesh_index", 0 );
        }
        else
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Parabola" );
            aParameterLists.insert( "variable_1", tTopXShift );
            aParameterLists.insert( "variable_2", tTopYShift );
            aParameterLists.insert( "variable_3", tTopParabolaFactor );
        }

        /* --------------------------------------- Bottom (flat) -------------------------------------- */
        if ( tBottomIsPlane )
        {
            real tNormalX = std::cos( tInterfaceBottomAngle + M_PI_2 );
            real tNormalY = std::sin( tInterfaceBottomAngle + M_PI_2 );

            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::LINE ) );
            aParameterLists.set( "center_x", tBottomXShift );
            aParameterLists.set( "center_y", tBottomYShift );
            aParameterLists.set( "normal_x", tNormalX );
            aParameterLists.set( "normal_y", tNormalY );
            aParameterLists.set( "number_of_refinements", 0 );
            aParameterLists.set( "refinement_mesh_index", 0 );
        }
        else
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "field_function_name", "Parabola" );
            aParameterLists.insert( "variable_1", tBottomXShift );
            aParameterLists.insert( "variable_2", tBottomYShift );
            aParameterLists.insert( "variable_3", tBottomParabolaFactor );
        }
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### FEM ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void FEMParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* -------------------------------------------- Phases ---------------------------------------- */
        Parameter_List pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseVoid" );
        pl.set( "phase_indices", "0" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "phase_indices", "1" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "phase_indices", "2" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        pl = prm::create_phase_parameter_list();
        pl.set( "phase_name", "PhaseAll" );
        pl.set( "phase_indices", "1,2" );
        aParameterLists( FEM::PHASES ).add_parameter_list( pl );

        /* ------------------------------------------ Properties -------------------------------------- */

        /* -------------------------------------- Youngs Modulus -------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsTop" );
        pl.set( "function_parameters", std::to_string( tTopEmod ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropYoungsBottom" );
        pl.set( "function_parameters", std::to_string( tBottomEmod ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* ------------------------------------------ Poisson ----------------------------------------- */
        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonTop" );
        pl.set( "function_parameters", std::to_string( tTopPois ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        pl = prm::create_property_parameter_list();
        pl.set( "property_name", "PropPoissonBottom" );
        pl.set( "function_parameters", std::to_string( tBottomPois ) );
        pl.set( "value_function", F2STR( Func_Const ) );
        aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );

        /* --------------------------------------- Dirichlet Top -------------------------------------- */
        // displacement-controlled vertical penetration on the top face (grows with load factor)
        if ( tTopHasDirichlet )
        {
            pl = prm::create_property_parameter_list();
            pl.set( "property_name", "PropDirichletTop" );
            pl.set( "function_parameters", "0.0;" + std::to_string( tVerticalDisplacement ) );
            pl.set( "value_function", F2STR( Func_Dirichlet ) );
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );
        }

        /* ------------------------------------- Dirichlet Bottom ------------------------------------- */
        // bottom face fully fixed
        if ( tBottomHasDirichlet )
        {
            pl = prm::create_property_parameter_list();
            pl.set( "property_name", "PropDirichletBottom" );
            pl.set( "function_parameters", "0.0;0.0" );
            pl.set( "value_function", F2STR( Func_Const ) );
            aParameterLists( FEM::PROPERTIES ).add_parameter_list( pl );
        }

        /* ---------------------------------------- Constitutive -------------------------------------- */

        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialTop" );
        pl.set( "phase_name", "PhaseTop" );
        pl.set( "constitutive_type", (uint)fem::Constitutive_Type::STRUC_LIN_ISO );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsTop,YoungsModulus;"
                "PropPoissonTop,PoissonRatio" );
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        pl = prm::create_constitutive_model_parameter_list();
        pl.set( "constitutive_name", "MaterialBottom" );
        pl.set( "phase_name", "PhaseBottom" );
        pl.set( "constitutive_type", (uint)fem::Constitutive_Type::STRUC_LIN_ISO );
        pl.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY", "Displacement" ) );
        pl.set( "properties",
                "PropYoungsBottom,YoungsModulus;"
                "PropPoissonBottom,PoissonRatio" );
        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list( pl );

        /* ------------------------------------- Stabilization ---------------------------------------- */

        // Dirichlet Nitsche SPs (top loading face / bottom fixed face)
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPDirichletNitscheBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", "100.0" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        /* ----------------------------- Contact Nitsche interface SPs ---------------------------- */
        // Wired into the contact IWG's "NitscheInterface" slot. DIRICHLET_NITSCHE-typed SP,
        // params "gamma;1.0" (interface penalty scaling), leader Material -> Youngs.
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", std::to_string( tContactStabilization ) + ";1.0" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPContactInterfaceBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::DIRICHLET_NITSCHE );
        pl.set( "function_parameters", std::to_string( tContactStabilization ) + ";1.0" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        // ghost penalty displacement (top)
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplTop" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_properties", "PropYoungsTop,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        // ghost penalty displacement (bottom)
        pl = prm::create_stabilization_parameter_parameter_list();
        pl.set( "stabilization_name", "SPGPDisplBottom" );
        pl.set( "stabilization_type", (uint)fem::Stabilization_Type::GHOST_DISPL );
        pl.set( "function_parameters", "0.05" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_properties", "PropYoungsBottom,Material" );
        aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

        /* ------------------------------------------- IWG -------------------------------------------- */

        /* -------------------------------------- Body Material Top ----------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructTop" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_BULK );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* ------------------------------------ Body Material Bottom ---------------------------------- */
        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGBulkStructBottom" );
        pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_BULK );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::BULK );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        /* ---------------------------------------- Dirichlet Top ------------------------------------- */
        if ( tTopHasDirichlet )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGDirichletTop" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "side_ordinals", "3" );    // north / top face
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropDirichletTop,Dirichlet" );
            pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheTop,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* -------------------------------------- Dirichlet Bottom ------------------------------------ */
        if ( tBottomHasDirichlet )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGDirichletBottom" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "side_ordinals", "1" );    // south / bottom face
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_properties", "PropDirichletBottom,Dirichlet" );
            pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
            pl.set( "stabilization_parameters", "SPDirichletNitscheBottom,DirichletNitsche" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* ------------------------------------- Contact Interface ------------------------------------ */
        // Neutral Nitsche, unbiased (theta = 0). FD jacobian, combined res/jac assembly (see SOL).
        IWG_Type tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED;

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGContactInterfaceTop" );
        pl.set( "IWG_type", (uint)tContactIWGType );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseTop" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseBottom" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "follower_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialTop,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPContactInterfaceTop,NitscheInterface" );
        pl.set( "analytical_jacobian", false );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        pl = prm::create_IWG_parameter_list();
        pl.set( "IWG_name", "IWGContactInterfaceBottom" );
        pl.set( "IWG_type", (uint)tContactIWGType );
        pl.set( "IWG_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
        pl.set( "dof_residual", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "follower_dof_dependencies", "UX,UY" );
        pl.set( "leader_constitutive_models", "MaterialBottom,ElastLinIso" );
        pl.set( "follower_constitutive_models", "MaterialTop,ElastLinIso" );
        pl.set( "stabilization_parameters", "SPContactInterfaceBottom,NitscheInterface" );
        pl.set( "analytical_jacobian", false );
        aParameterLists( FEM::IWG ).add_parameter_list( pl );

        if ( tUseGhost )
        {
            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGGPDisplTop" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "follower_phase_name", "PhaseTop" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "follower_dof_dependencies", "UX,UY" );
            pl.set( "stabilization_parameters", "SPGPDisplTop,GhostSP" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );

            pl = prm::create_IWG_parameter_list();
            pl.set( "IWG_name", "IWGGPDisplBottom" );
            pl.set( "IWG_type", (uint)fem::IWG_Type::GHOST_NORMAL_FIELD );
            pl.set( "IWG_bulk_type", (uint)fem::Element_Type::DOUBLE_SIDESET );
            pl.set( "dof_residual", "UX,UY" );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "follower_phase_name", "PhaseBottom" );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "follower_dof_dependencies", "UX,UY" );
            pl.set( "stabilization_parameters", "SPGPDisplBottom,GhostSP" );
            aParameterLists( FEM::IWG ).add_parameter_list( pl );
        }

        /* ------------------------------------------- IQI -------------------------------------------- */

        pl = prm::create_IQI_parameter_list();
        pl.set( "IQI_name", "IQIDispX" );
        pl.set( "IQI_type", (uint)fem::IQI_Type::DOF );
        pl.set( "dof_quantity", "UX,UY" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseAll" );
        pl.set( "vectorial_field_index", 0 );
        aParameterLists( FEM::IQI ).add_parameter_list( pl );

        pl = prm::create_IQI_parameter_list();
        pl.set( "IQI_name", "IQIDispY" );
        pl.set( "IQI_type", (uint)fem::IQI_Type::DOF );
        pl.set( "dof_quantity", "UX,UY" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "leader_phase_name", "PhaseAll" );
        pl.set( "vectorial_field_index", 1 );
        aParameterLists( FEM::IQI ).add_parameter_list( pl );

        // FULL 2D stress tensor per body (for sigma_nn recovery in post via the parabola normal).
        //   NORMAL_STRESS vectorial index 0 -> sigma_xx  ("XX")
        //   NORMAL_STRESS vectorial index 1 -> sigma_yy  ("YY")
        //   SHEAR_STRESS  vectorial index 2 -> sigma_xy  ("XY")
        // NOTE: SHEAR_STRESS reads the standardized 6-vector at slot (index+3); in 2D the in-plane
        //   shear sigma_xy lives at slot 5, so vectorial_field_index MUST be 2 (not 0). Both stress
        //   IQI types require the ElastLinIso CM slot (they pull the CM flux vector).
        Vector< std::string > tSides = { "Top", "Bottom" };
        struct StressSpec { std::string tAxis; uint tType; int tIndex; };
        Vector< StressSpec > tStressSpecs = {
            { "XX", (uint)fem::IQI_Type::NORMAL_STRESS, 0 },
            { "YY", (uint)fem::IQI_Type::NORMAL_STRESS, 1 },
            { "XY", (uint)fem::IQI_Type::SHEAR_STRESS, 2 },
        };
        for ( auto const &tSide : tSides )
        {
            for ( auto const &tSpec : tStressSpecs )
            {
                pl = prm::create_IQI_parameter_list();
                pl.set( "IQI_name", "IQIStress" + tSide + tSpec.tAxis );
                pl.set( "IQI_type", tSpec.tType );
                pl.set( "leader_dof_dependencies", "UX,UY" );
                pl.set( "leader_constitutive_models", "Material" + tSide + ",ElastLinIso" );
                pl.set( "leader_phase_name", "Phase" + tSide );
                pl.set( "vectorial_field_index", tSpec.tIndex );
                aParameterLists( FEM::IQI ).add_parameter_list( pl );
            }

            // contact pressure (reference config). DEFINED for post-measurement of P; NOT output in
            // VIS (FACETED on the interface side set, incompatible with the bulk-only unzipped mesh).
            pl = prm::create_IQI_parameter_list();
            pl.set( "IQI_name", "IQIContactPressure" + tSide );
            pl.set( "IQI_type", (uint)fem::IQI_Type::CONTACT_PRESSURE_REFERENCE );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_phase_name", "Phase" + tSide );
            pl.set( "leader_constitutive_models", "Material" + tSide + ",TractionCM" );
            aParameterLists( FEM::IQI ).add_parameter_list( pl );
        }

        // gap (ray length) across the non-conforming interface. Output on the SECOND (STANDARD-type)
        // VIS mesh over the nonconformal contact interface set (see VISParameterList) -- it CANNOT go
        // on the FULL_DISCONTINUOUS bulk mesh (nonconformal side set crashes populate_nodal_point_pairs).
        pl = prm::create_IQI_parameter_list();
        pl.set( "IQI_name", "IQIGap" );
        pl.set( "IQI_type", (uint)fem::IQI_Type::GAP );
        pl.set( "IQI_bulk_type", (uint)fem::Element_Type::NONCONFORMAL_SIDESET );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "neighbor_phases", "PhaseVoid" );
        pl.set( "follower_phase_name", "PhaseTop" );
        pl.set( "leader_dof_dependencies", "UX,UY" );
        pl.set( "follower_dof_dependencies", "UX,UY" );
        aParameterLists( FEM::IQI ).add_parameter_list( pl );

        /* --------------------------------------- Computation ---------------------------------------- */
        aParameterLists( FEM::COMPUTATION ).set( "is_analytical_forward", true );
        aParameterLists.set( "finite_difference_scheme_forward", (uint)( tFDScheme ) );
        aParameterLists.set( "finite_difference_perturbation_size_forward", tFDPerturbationSize );
        aParameterLists.set( "finite_difference_perturbation_strategy", (uint)tFDPerturbationStrategy );
        aParameterLists.set( "nonconformal_integration_order", static_cast< uint >( tNonconformalIntegrationOrder ) );
        aParameterLists.set( "nonconformal_max_negative_ray_length", tMaxNegativeRayLength );
        aParameterLists.set( "nonconformal_max_positive_ray_length", tMaxPositiveRayLength );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### SOL ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void SOLParameterList( Module_Parameter_Lists &aParameterLists )
    {
        /* ------------------------------------- Linear Algorithm ------------------------------------- */
        Parameter_List pl = prm::create_linear_algorithm_parameter_list( sol::SolverType::AMESOS_IMPL );
        pl.set( "Solver_Type", "Amesos_Umfpack" );
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( pl );

        /* -------------------------------------- Linear Solver --------------------------------------- */
        pl = prm::create_linear_solver_parameter_list();
        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list( pl );

        /* ------------------------------------ Nonlinear Algorithm ----------------------------------- */
        pl = prm::create_nonlinear_algorithm_parameter_list();
        // Neutral-Nitsche contact IWG implements compute_jacobian_and_residual -> combined assembly.
        pl.set( "NLA_combined_res_jac_assembly", true );
        pl.set( "NLA_Linear_solver", 0 );
        pl.set( "NLA_max_iter", tMaxIterations );
        pl.set( "NLA_rel_res_norm_drop", tRelResNormDrop );
        pl.set( "NLA_load_control_strategy", (uint)( sol::SolverLoadControlType::Linear ) );
        pl.set( "NLA_load_control_factor", tLoadControlFactor );
        pl.set( "NLA_load_control_steps", tLoadControlSteps );
        pl.set( "NLA_load_control_relres", tLoadControlRelRes );
        pl.set( "NLA_relaxation_parameter", tRelaxation );
        pl.set( "NLA_remap_strategy", (uint)tRaytracingStrategy );
        pl.set( "NLA_remap_load_step_frequency", tRemapLoadStepFrequency );
        pl.set( "NLA_remap_iteration_frequency", tRemapIterationFrequency );
        pl.set( "NLA_remap_residual_change_tolerance", tRemapResidualChange );
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list( pl );

        /* ------------------------------------- Nonlinear Solver ------------------------------------- */
        pl = prm::create_nonlinear_solver_parameter_list();
        pl.set( "NLA_DofTypes", "UX,UY" );
        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list( pl );

        /* --------------------------------------- Time Algorithm ------------------------------------- */
        pl = prm::create_time_solver_algorithm_parameter_list();
        pl.set( "TSA_Num_Time_Steps", 1 );
        pl.set( "TSA_Time_Frame", 10.0 );
        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list( pl );

        /* ---------------------------------------- Time Solver --------------------------------------- */
        pl = prm::create_time_solver_parameter_list();
        pl.set( "TSA_DofTypes", "UX,UY" );
        // Trigger BOTH VIS output meshes: index 0 = bulk stress (FULL_DISCONTINUOUS),
        // index 1 = contact gap (STANDARD). The time solver only writes indices listed here,
        // and SOL_Warehouse hard-asserts len(Indices)==len(Criteria), so both lists grow together.
        pl.set( "TSA_Output_Indices", "0,1" );
        pl.set( "TSA_Output_Criteria", F2STR( Output_Criterion ) "," F2STR( Output_Criterion ) );
        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list( pl );

        /* -------------------------------------- Preconditioner -------------------------------------- */
        pl = prm::create_preconditioner_parameter_list( sol::PreconditionerType::NONE );
        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( pl );
    }

    void MSIParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
    }

    /* ---------------------------------------------------------------------------------------------- */
    /*                                           ### VIS ###                                          */
    /* ---------------------------------------------------------------------------------------------- */
    void VISParameterList( Module_Parameter_Lists &aParameterLists )
    {
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        // Unzipped (per-phase discontinuous) VIS mesh: separate vertices per cut cell/phase so the
        // XFEM-discontinuous stress is NOT averaged across the interface -> exact per-side interface
        // stress on the indenter vs half-space. sigma_nn is recovered in post from these fields.
        aParameterLists.set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::FULL_DISCONTINUOUS );

        // Bulk domains ONLY (no nonconformal interface side set): FULL_DISCONTINUOUS + a nonconformal
        // side set in Set_Names crashes VIS_Factory::populate_nodal_point_pairs (core bug on the
        // discontinuous vertex path). The per-side interface stress still comes from the bulk fields
        // on the unzipped mesh (top/bottom blocks are separate meshes with their own cut cells).
        std::string tSetNames = tDomain;

        // displacement components (nodal)
        std::string tFieldNames = "UX,UY";
        std::string tFieldTypes = "NODAL,NODAL";
        std::string tIQINames   = "IQIDispX,IQIDispY";

        // FULL stress tensor per body, both NODAL (STRESSNodal<Side><Comp>) and ELEMENTAL_AVG
        // (STRESSBulk<Side><Comp>) -- so sigma_nn = nx^2*sxx + 2*nx*ny*sxy + ny^2*syy can be
        // recovered in post from the same IQIs with the analytic parabola normal.
        Vector< std::string > tSides      = { "Top", "Bottom" };
        Vector< std::string > tComponents = { "XX", "YY", "XY" };

        for ( auto const &tSide : tSides )
        {
            for ( auto const &tComp : tComponents )
            {
                // elemental-averaged stress
                tFieldNames += ",STRESSBulk" + tSide + tComp;
                tFieldTypes += ",ELEMENTAL_AVG";
                tIQINames += ",IQIStress" + tSide + tComp;
            }
        }
        for ( auto const &tSide : tSides )
        {
            for ( auto const &tComp : tComponents )
            {
                // nodal-projected stress (cleaner on cut cells; consumed by check_results by name)
                tFieldNames += ",STRESSNodal" + tSide + tComp;
                tFieldTypes += ",NODAL";
                tIQINames += ",IQIStress" + tSide + tComp;
            }
        }

        aParameterLists.set( "Set_Names", tSetNames );
        aParameterLists.set( "Field_Names", tFieldNames );
        aParameterLists.set( "Field_Type", tFieldTypes );
        aParameterLists.set( "IQI_Names", tIQINames );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 0.1 );
        // Explicit output index: the VIS Output Manager keys meshes by Output_Index
        // (mOutputData(mMeshIndex)); it defaults to 0, so a second mesh MUST get its own
        // index or it overwrites this one (cl_VIS_Output_Manager.cpp:190).
        aParameterLists.set( "Output_Index", 0 );

        /* ------------------------------------------------------------------------------------------ */
        /* SECOND output mesh: contact GAP on the nonconformal contact interface set.                 */
        /*   The bulk mesh above is FULL_DISCONTINUOUS and MUST stay bulk-only (a nonconformal side   */
        /*   set on a FULL_DISCONTINUOUS mesh crashes VIS_Factory::populate_nodal_point_pairs). So the */
        /*   interface-only GAP field lives on its own STANDARD-type mesh, written to <name>_gap.exo. */
        /*   Set_Names = bulk blocks (tDomain) + the ncss|iside_b0_2_b1_0|iside_b0_1_b1_0 double-sided  */
        /*   nonconformal pair the neutral-theta contact Nitsche IWG runs on. The bulk blocks MUST be   */
        /*   included so the double-sided contact cluster's LEADER cell (which lives in a bulk block     */
        /*   set) resolves to a VIS cell -- omitting them trips                                          */
        /*   VIS_Factory::create_visualization_dbl_side_clusters (cl_VIS_Factory.cpp:1121). GAP is only  */
        /*   defined on tContactInterface; the bulk blocks carry no GAP field, exactly as in the proven  */
        /*   Parabolic_Indenter STANDARD reference (Set_Names = tDomain + "," + tContactInterface, a     */
        /*   single GAP/FACETED_AVG/IQIGap triple). GAP (leader=PhaseBottom, follower=PhaseTop) renders  */
        /*   on the bottom-leader half of the pair.                                                      */
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName + "_gap" ) );
        aParameterLists.set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tDomain + "," + tContactInterface );
        aParameterLists.set( "Field_Names", "GAP" );
        aParameterLists.set( "Field_Type", "FACETED_AVG" );
        aParameterLists.set( "IQI_Names", "IQIGap" );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 0.1 );
        // Distinct output index so this gap mesh does not overwrite the bulk mesh above
        // (both default to Output_Index 0 -> mOutputData collision, only one survives).
        aParameterLists.set( "Output_Index", 1 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists &aParameterLists ) {}
}    // namespace moris
#ifdef __cplusplus
}
#endif
