/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Struc_Contact_Patch_Test.cpp
 *
 * Small-strain constant-stress contact patch test.
 *
 * Two FLAT linear-elastic blocks (E=1000, nu=0) are stacked with a FLAT immersed
 * non-conforming interface and compressed uniformly in the vertical direction.
 * With nu=0 the exact solution is a spatially UNIFORM uniaxial normal stress
 * sigma_yy = applied traction (zero tangential). The section asserts sigma_yy is
 * uniform across the non-conforming interface (patch-test criterion).
 *
 * Adapted from Parabolic_Indenter/Parabolic_Indenter_Linear.cpp (phase-based FEM API).
 *
 * The contact METHOD/BIAS and stabilization gamma are selected by the deck-visible global
 * gCaseIndex = method*4 + gamma  (0..15), a 4 (method) x 4 (gamma) matrix:
 *
 *   method = gCaseIndex / 4:
 *     0 -> nitsche, symmetric   (STRUC_LINEAR_CONTACT_NORMAL_SYMMETRIC_NITSCHE_UNBIASED)
 *     1 -> nitsche, unsymmetric (STRUC_LINEAR_CONTACT_NORMAL_UNSYMMETRIC_NITSCHE_UNBIASED)
 *     2 -> nitsche, neutral     (STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED)
 *     3 -> penalty              (STRUC_LINEAR_CONTACT_PENALTY)
 *
 *   gamma  = gCaseIndex % 4:  0 -> 10, 1 -> 100, 2 -> 1000, 3 -> 10000
 *     gamma (tContactStabilization) drives BOTH contact SP function_parameters:
 *     the Nitsche interface SP ("STAB;1.0") and the penalty SP ("STAB"). Sweeping gamma
 *     is thesis Argument B: post-gap-closure the only term coupling the two bodies' stress
 *     is gamma*jump, so larger gamma should collapse the sigma_yy split toward uniform.
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

// global variable for case index (selects contact method/bias, see header comment)
extern uint gCaseIndex;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
namespace moris
{
    /* -------------------------------------- Geometry Setting -------------------------------------- */
    real tHeight = 1.0;
    real tWidth  = 1.0;

    uint tNumElems  = 10;
    uint tNumElemsX = tNumElems;
    uint tNumElemsY = tNumElems;

    // small initial gap band (PhaseVoid) between the two flat blocks; small penetration
    // to keep the response in the small-strain regime.
    real tInitialGap           = 0.001;
    real tInitialGapMiddle     = 0.503 * tHeight;
    real tPenetration          = 0.005;
    real tVerticalDisplacement = -tInitialGap - tPenetration;

    /* --------------------------------------- Top Geometry ----------------------------------------- */
    // PATCH TEST: flat top block, flat horizontal interface.
    bool tTopIsPlane = true;

    real tInterfaceTopAngle = 0 * M_PI / 180.0;
    real tTopYShift         = tInitialGapMiddle + tInitialGap / 2;
    real tTopXShift         = 0.5;
    real tTopParabolaFactor = 1;

    /* --------------------------------------- Bottom Geometry -------------------------------------- */
    // PATCH TEST: flat bottom block, flat horizontal interface.
    bool tBottomIsPlane = true;

    real tInterfaceBottomAngle = 0 * M_PI / 180.0;
    real tBottomYShift         = tInitialGapMiddle - tInitialGap / 2;
    real tBottomXShift         = 0.5;
    real tBottomParabolaFactor = -0.7;

    /* ------------------------------------- Material Parameters ------------------------------------ */
    // nu = 0 so uniaxial compression yields a spatially constant sigma_yy (patch test).
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
            "Struc_Contact_Patch_Test_Case_" + std::to_string( gCaseIndex );

    /* ------------------------------------------- Solver ------------------------------------------- */
    int tMaxIterations = 200;

    real tRelResNormDrop     = 1e-8;
    real tRelaxation         = 1.0;
    real tFDPerturbationSize = 1e-8;

    fem::Perturbation_Type tFDPerturbationStrategy = fem::Perturbation_Type::ABSOLUTE;
    fem::FDScheme_Type     tFDScheme               = fem::FDScheme_Type::POINT_3_CENTRAL;

    int  tLoadControlSteps    = 5;
    real tLoadControlFactor   = 0.2;
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
    // gCaseIndex = method*4 + gamma  (0..15). Decode a 4 (method) x 4 (gamma) matrix.
    uint tContactMethodIndex = gCaseIndex / 4;    // 0..3 : method
    uint tContactGammaIndex  = gCaseIndex % 4;    // 0..3 : stabilization gamma

    // Contact method / bias are selected from the method index (see header comment).
    //   method 3   -> "penalty",  otherwise "nitsche"
    //   method 0/1/2 -> symmetric / unsymmetric / neutral bias (nitsche only)
    std::string tContactMethod =
            ( tContactMethodIndex == 3 ) ? std::string( "penalty" ) : std::string( "nitsche" );
    std::string tContactBias =
            ( tContactMethodIndex == 0 )   ? std::string( "symmetric" )
            : ( tContactMethodIndex == 1 ) ? std::string( "unsymmetric" )
            : ( tContactMethodIndex == 2 ) ? std::string( "neutral" )
                                           : std::string( "unsymmetric" );

    // gamma sweep {10, 100, 1000, 10000}, selected from the gamma index. This single
    // per-case value drives BOTH contact SP function_parameters (Nitsche "STAB;1.0" and
    // penalty "STAB") -- no hardcoded 10.0 in the SP definitions any more.
    real tContactStabilization =
            ( tContactGammaIndex == 0 )   ? 10.0
            : ( tContactGammaIndex == 1 ) ? 100.0
            : ( tContactGammaIndex == 2 ) ? 1000.0
                                          : 10000.0;    // Nitsche interface penalty scaling ("STAB")
    real tContactPenaltyFactor = tContactStabilization;    // penalty-contact SP factor (mParameters(0)), same per-case gamma

    // convenience flag
    bool tIsPenalty = ( tContactMethod == "penalty" );

    /* -------------------------------------- Control Variables ------------------------------------- */
    bool tOnlyGenerateMesh   = false;
    bool tUseGhost           = true;
    bool tTopHasDirichlet    = true;
    bool tBottomHasDirichlet = true;

    /* ---------------------------------------------------------------------------------------------- */
    /*                                         Field Functions                                        */
    /* ---------------------------------------------------------------------------------------------- */

    /* --------------------------------------- Phase Indexing --------------------------------------- */
    uint Phase_Index_Split( const Bitset< 2 > &aGeometrySigns )
    {
        if ( aGeometrySigns.test( 1 ) )
        {
            if ( aGeometrySigns.test( 0 ) )
            {
                return 1;    // upper domain
            }
            return 0;    // interface (thin void band)
        }
        return 2;    // bottom
    }

    /* ----------------------------------- Geometry Field Function ---------------------------------- */

    real tMinimumLevelSetValue = 1.0e-8;

    // Level set function for a parabola that splits the domain into a top and a bottom part.
    // Retained so tTopIsPlane/tBottomIsPlane can be toggled off for debugging; the patch
    // test itself uses the flat LINE geometries (both *IsPlane == true).
    real Parabola( const Matrix< DDRMat > &aCoordinates,
            const Vector< real >          &aGeometryParameters )
    {

        real tXShift = aGeometryParameters( 0 );
        real tYShift = aGeometryParameters( 1 );
        real tFactor = aGeometryParameters( 2 );

        real x = aCoordinates( 0 ) - tXShift;
        real y = aCoordinates( 1 ) - tYShift;

        // Compute Signed-Distance field
        real tVal = y - ( tFactor * x * x );

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

        /* ----------------------------------------- Top Plane ---------------------------------------- */
        if ( tTopIsPlane )
        {
            // calculate the components of the rotated normal
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

        /* --------------------------------------- Bottom Plane --------------------------------------- */
        if ( tBottomIsPlane )
        {
            // calculate the components of the rotated normal
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
        // uniform vertical compression on the top face (grows with the load factor)
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

        // Dirichlet Nitsche SPs (top/bottom loading faces)
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

        if ( !tIsPenalty )
        {
            /* ----------------------------- Contact Nitsche interface SPs ---------------------------- */
            // Wired into the contact IWG's "NitscheInterface" slot. DIRICHLET_NITSCHE-typed SP,
            // params "STAB;1.0" (interface penalty scaling), leader Material -> Youngs.
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
        }
        else
        {
            /* ----------------------------- Contact penalty SPs -------------------------------------- */
            // Wired into the penalty IWG's "PenaltyContact" slot. PENALTY_CONTACT-typed SP,
            // params = penalty factor (mParameters(0)); requires BOTH leader and follower
            // "Material" -> Youngs (SP_Penalty_Contact::eval_SP reads leader+follower E).
            pl = prm::create_stabilization_parameter_parameter_list();
            pl.set( "stabilization_name", "SPContactPenaltyTop" );
            pl.set( "stabilization_type", (uint)fem::Stabilization_Type::PENALTY_CONTACT );
            pl.set( "function_parameters", std::to_string( tContactPenaltyFactor ) );
            pl.set( "leader_phase_name", "PhaseTop" );
            pl.set( "follower_phase_name", "PhaseBottom" );
            pl.set( "leader_properties", "PropYoungsTop,Material" );
            pl.set( "follower_properties", "PropYoungsBottom,Material" );
            aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );

            pl = prm::create_stabilization_parameter_parameter_list();
            pl.set( "stabilization_name", "SPContactPenaltyBottom" );
            pl.set( "stabilization_type", (uint)fem::Stabilization_Type::PENALTY_CONTACT );
            pl.set( "function_parameters", std::to_string( tContactPenaltyFactor ) );
            pl.set( "leader_phase_name", "PhaseBottom" );
            pl.set( "follower_phase_name", "PhaseTop" );
            pl.set( "leader_properties", "PropYoungsBottom,Material" );
            pl.set( "follower_properties", "PropYoungsTop,Material" );
            aParameterLists( FEM::STABILIZATION ).add_parameter_list( pl );
        }

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
        // Determine the requested contact IWG type from the (method, bias) selection.
        IWG_Type tContactIWGType = IWG_Type::END_IWG_TYPE;
        if ( tContactMethod == "penalty" )
        {
            tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_PENALTY;
        }
        else    // nitsche
        {
            if ( tContactBias == "symmetric" )
                tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_SYMMETRIC_NITSCHE_UNBIASED;
            if ( tContactBias == "unsymmetric" )
                tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_UNSYMMETRIC_NITSCHE_UNBIASED;
            if ( tContactBias == "neutral" )
                tContactIWGType = IWG_Type::STRUC_LINEAR_CONTACT_NORMAL_NEUTRAL_NITSCHE_UNBIASED;
        }

        // SP name + IWG slot differ between methods:
        //   nitsche -> ("SPContactInterface*",  "NitscheInterface"), FD jacobian
        //   penalty -> ("SPContactPenalty*",    "PenaltyContact"),   analytic jacobian
        std::string tContactSPTop     = tIsPenalty ? "SPContactPenaltyTop,PenaltyContact" : "SPContactInterfaceTop,NitscheInterface";
        std::string tContactSPBottom  = tIsPenalty ? "SPContactPenaltyBottom,PenaltyContact" : "SPContactInterfaceBottom,NitscheInterface";
        bool        tAnalyticJacobian = tIsPenalty;    // penalty IWG implements an analytic jacobian; Nitsche uses FD

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
        pl.set( "stabilization_parameters", tContactSPTop );
        pl.set( "analytical_jacobian", tAnalyticJacobian );
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
        pl.set( "stabilization_parameters", tContactSPBottom );
        pl.set( "analytical_jacobian", tAnalyticJacobian );
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

        // Normal stress (vectorial index 1 == sigma_yy) + contact pressure, per side.
        // sigma_yy is the patch-test quantity: it must be spatially UNIFORM.
        Vector< std::string >                   tSides      = { "Top", "Bottom" };
        Vector< std::pair< int, std::string > > tComponents = { { 0, "X" }, { 1, "Y" } };
        for ( auto const &tSide : tSides )
        {
            for ( auto const &[ tDir, tAxis ] : tComponents )
            {
                pl = prm::create_IQI_parameter_list();
                pl.set( "IQI_name", "IQIStress" + tSide + tAxis );
                pl.set( "IQI_type", (uint)fem::IQI_Type::NORMAL_STRESS );
                pl.set( "leader_dof_dependencies", "UX,UY" );
                pl.set( "leader_constitutive_models", "Material" + tSide + ",ElastLinIso" );
                pl.set( "leader_phase_name", "Phase" + tSide );
                pl.set( "vectorial_field_index", tDir );
                aParameterLists( FEM::IQI ).add_parameter_list( pl );
            }

            pl = prm::create_IQI_parameter_list();
            pl.set( "IQI_name", "IQIContactPressure" + tSide );
            pl.set( "IQI_type", (uint)fem::IQI_Type::CONTACT_PRESSURE_REFERENCE );
            pl.set( "leader_dof_dependencies", "UX,UY" );
            pl.set( "leader_phase_name", "Phase" + tSide );
            pl.set( "leader_constitutive_models", "Material" + tSide + ",TractionCM" );
            aParameterLists( FEM::IQI ).add_parameter_list( pl );
        }

        // gap (ray length) across the non-conforming interface
        pl = prm::create_IQI_parameter_list();
        pl.set( "IQI_name", "IQIGap" );
        pl.set( "IQI_type", (uint)fem::IQI_Type::GAP );
        pl.set( "leader_phase_name", "PhaseBottom" );
        pl.set( "follower_phase_name", "PhaseTop" );
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
        // CRITICAL: the penalty contact IWG does NOT implement compute_jacobian_and_residual
        // (it errors). Assemble residual and jacobian separately for penalty; keep the combined
        // path for the Nitsche variants (which do implement it).
        pl.set( "NLA_combined_res_jac_assembly", !tIsPenalty );
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
        pl.set( "TSA_Output_Indices", "0" );
        pl.set( "TSA_Output_Criteria", F2STR( Output_Criterion ) );
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
        // XFEM-discontinuous stress is NOT averaged across the interface — gives the EXACT per-side
        // interface sigma_yy (top-phase vs bottom-phase) rather than a blended nodal value.
        aParameterLists.set( "Mesh_Type", (uint)vis::VIS_Mesh_Type::FULL_DISCONTINUOUS );

        // Output over the bulk domains AND the non-conforming contact interface so a
        // sigma_yy uniformity (max-min) check is possible across the contact region.
        // Bulk domains ONLY (no nonconformal interface side set): FULL_DISCONTINUOUS + a
        // nonconformal side set in Set_Names crashes VIS_Factory::populate_nodal_point_pairs
        // (core bug — the discontinuous vertex path leaves the FEM->VIS map empty). The exact
        // per-side interface sigma_yy still comes from the bulk fields on the unzipped mesh
        // (top/bottom blocks are separate meshes with their own interface-adjacent cut cells).
        std::string tSetNames   = tDomain;

        // displacement components (GAP dropped — it is FACETED on the interface side set)
        std::string tFieldNames = "UX,UY";
        std::string tFieldTypes = "NODAL,NODAL";
        std::string tIQINames   = "IQIDispX,IQIDispY";

        // normal stress components (STRESSBulk*X, STRESSBulk*Y == sigma_xx, sigma_yy) + contact pressure
        Vector< std::string > tSides      = { "Top", "Bottom" };
        Vector< std::string > tComponents = { "X", "Y" };

        for ( auto const &tSide : tSides )
        {
            for ( auto const &tComponent : tComponents )
            {
                tFieldNames += ",STRESSBulk" + tSide + tComponent;
                tFieldTypes += ",ELEMENTAL_AVG";
                tIQINames += ",IQIStress" + tSide + tComponent;
            }
            // ContactPressure dropped — FACETED on the interface side set (see Set_Names note).
        }

        // Cleaner sigma_yy for the patch-test metric: NODAL-projected normal stress.
        // There is no dedicated nodal/projected sigma_yy IQI type (the only recovered-stress
        // IQI, ZIENKIEWICZ_ZHU_VON_MISES_STRESS, is a von-Mises scalar), so we reuse the
        // NORMAL_STRESS IQIs IQIStress{Top,Bottom}Y but output them with NODAL projection,
        // which is more reliable on cut cells than the ELEMENTAL_AVG above. Exodus_IO_Helper
        // reads nodal fields only; example_test_case.cpp::check_results() consumes these by
        // name (STRESSNodalTopY / STRESSNodalBottomY) for the interior sigma_yy spread.
        for ( auto const &tSide : tSides )
        {
            tFieldNames += ",STRESSNodal" + tSide + "Y";
            tFieldTypes += ",NODAL";
            tIQINames += ",IQIStress" + tSide + "Y";
        }

        aParameterLists.set( "Set_Names", tSetNames );
        aParameterLists.set( "Field_Names", tFieldNames );
        aParameterLists.set( "Field_Type", tFieldTypes );
        aParameterLists.set( "IQI_Names", tIQINames );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Time_Offset", 0.1 );
    }

    void MORISGENERALParameterList( Module_Parameter_Lists &aParameterLists ) {}
}    // namespace moris
#ifdef __cplusplus
}
#endif
