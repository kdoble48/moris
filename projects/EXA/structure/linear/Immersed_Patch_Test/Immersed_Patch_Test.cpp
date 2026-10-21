/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Immersed_Patch_Test.cpp
 *
 * Linear patch test for the immersed (XTK) discretization. A 2 x 2 x 2 HEX8 background mesh on [0,2]^3 is cut by two
 * parallel planes tDelta apart, giving three phases with the SAME isotropic linear-elastic material. The exact linear
 * displacement u = A x + b is imposed weakly (symmetric Nitsche) on all outer faces and the phases are coupled by
 * symmetric Nitsche interface conditions, so the exact solution is u in every phase and the stress is uniform.
 *
 * With tDelta = 1e-6 the thin layer between the planes produces IG cells that the FEM cluster drops from the residual
 * (relative volume below 1e-6). The visualized element stress of those cells must still be the exact value.
 *
 */

#include <string>
#include <iostream>
#include <cmath>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_MSI_Equation_Object.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "cl_DLA_Solver_Interface.hpp"
#include "cl_DLA_Linear_Solver_Aztec.hpp"
#include "parameters.hpp"
#include "fn_PRM_MIG_Parameters.hpp"
#include "cl_HMR_Element.hpp"
#include "fn_equal_to.hpp"

#include "AztecOO.h"

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // geometry and stabilization

    real tDelta    = 1.0e-6;    // distance between the two planes
    bool tUseGhost = true;      // ghost stabilization

    /* ------------------------------------------------------------------------ */
    // exact solution u = A x + b

    real tA[ 3 ][ 3 ] = { { 1.0e-3, 2.0e-4, 0.0 }, { 1.0e-4, -5.0e-4, 3.0e-4 }, { 0.0, 2.0e-4, 4.0e-4 } };
    real tB[ 3 ]      = { 1.0e-3, -2.0e-3, 5.0e-4 };

    /* ------------------------------------------------------------------------ */
    // output

    std::string tOutputFileName = "Immersed_Patch_Test.exo";

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    // Level set A: signed distance to the plane through (0.95, 1, 1) with normal (1, 0.3, 0.2)
    moris::real
    Geometry_A(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        moris::real tNx   = 1.0;
        moris::real tNy   = 0.3;
        moris::real tNz   = 0.2;
        moris::real tNorm = std::sqrt( tNx * tNx + tNy * tNy + tNz * tNz );

        moris::real tDx = aCoordinates( 0 ) - 0.95;
        moris::real tDy = aCoordinates( 1 ) - 1.0;
        moris::real tDz = aCoordinates( 2 ) - 1.0;

        return ( tNx * tDx + tNy * tDy + tNz * tDz ) / tNorm;
    }

    /* ------------------------------------------------------------------------ */

    // Level set B: the parallel plane tDelta further along the normal
    moris::real
    Geometry_B(
            const moris::Matrix< DDRMat >& aCoordinates,
            const Vector< real >&          aGeometryParameters )
    {
        return Geometry_A( aCoordinates, aGeometryParameters ) - tDelta;
    }

    /* ------------------------------------------------------------------------ */

    void
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "print_phase_table", true );

        // The phase entry is 2 * sign_0 + sign_1, so geometry 0 is the high bit. Registering B before A gives
        // phase 0 = A-,B-; phase 1 = A+,B- (the layer); phase 2 = A-,B+ (empty); phase 3 = A+,B+.
        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "name", "GeomB" );
        aParameterLists.set( "field_function_name", "Geometry_B" );
        aParameterLists.set( "discretization_mesh_index", -1 );
        aParameterLists.set( "use_multilinear_interpolation", false );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );

        aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
        aParameterLists.set( "name", "GeomA" );
        aParameterLists.set( "field_function_name", "Geometry_A" );
        aParameterLists.set( "discretization_mesh_index", -1 );
        aParameterLists.set( "use_multilinear_interpolation", false );
        aParameterLists.set( "isocontour_threshold", 0.0 );
        aParameterLists.set( "isocontour_tolerance", 1.0e-12 );
        aParameterLists.set( "intersection_tolerance", 1.0e-12 );
    }

    /* ------------------------------------------------------------------------ */

    void
    Func_Const( moris::Matrix< moris::DDRMat >&           aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >& aParameters,
            moris::fem::Field_Interpolator_Manager*   aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    /* ------------------------------------------------------------------------ */

    // exact linear displacement u = A x + b (3x1) at the evaluation point
    void
    Func_Exact_Displacement( moris::Matrix< moris::DDRMat >& aPropMatrix,
            Vector< moris::Matrix< moris::DDRMat > >&       aParameters,
            moris::fem::Field_Interpolator_Manager*         aFIManager )
    {
        const moris::Matrix< moris::DDRMat >& tXp = aFIManager->get_IP_geometry_interpolator()->valx();

        aPropMatrix.set_size( 3, 1, 0.0 );

        for ( moris::uint iI = 0; iI < 3; iI++ )
        {
            moris::real tValue = tB[ iI ];
            for ( moris::uint iJ = 0; iJ < 3; iJ++ )
            {
                tValue += tA[ iI ][ iJ ] * tXp( iJ );
            }
            aPropMatrix( iI ) = tValue;
        }
    }

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();

        //------------------------------------------------------------------------------
        // mesh sets: every bulk phase gets the same physics (phase 2 is empty, its sets exist and are harmless)

        const Vector< std::string > tPhases = { "0", "1", "2", "3" };

        std::string tBulkSets;
        std::string tSideSets;
        std::string tGhostSets;

        for ( moris::uint iP = 0; iP < tPhases.size(); iP++ )
        {
            const std::string& tP = tPhases( iP );

            if ( iP > 0 )
            {
                tBulkSets += ",";
                tGhostSets += ",";
            }
            tBulkSets += "HMR_dummy_n_p" + tP + ",HMR_dummy_c_p" + tP;
            tGhostSets += "ghost_p" + tP;

            for ( moris::uint iS = 1; iS <= 6; iS++ )
            {
                if ( !tSideSets.empty() )
                {
                    tSideSets += ",";
                }
                std::string tS = std::to_string( iS );
                tSideSets += "SideSet_" + tS + "_n_p" + tP + ",SideSet_" + tS + "_c_p" + tP;
            }
        }

        // every phase pair, low-to-high (leader = lower phase), one orientation each
        std::string tInterfaceSets =
                "dbl_iside_p0_0_p1_1,dbl_iside_p0_0_p1_2,dbl_iside_p0_0_p1_3,"
                "dbl_iside_p0_1_p1_2,dbl_iside_p0_1_p1_3,dbl_iside_p0_2_p1_3";

        //------------------------------------------------------------------------------
        // PROPERTIES

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropYoung" );
        aParameterLists.set( "function_parameters", "1.0" );
        aParameterLists.set( "value_function", "Func_Const" );

        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropPoisson" );
        aParameterLists.set( "function_parameters", "0.3" );
        aParameterLists.set( "value_function", "Func_Const" );

        // exact displacement: Dirichlet value on the outer sides and the L2-error reference
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropExactDisp" );
        aParameterLists.set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Exact_Displacement" );

        //------------------------------------------------------------------------------
        // CONSTITUTIVE MODELS (one CM for every phase)

        aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
        aParameterLists.set( "constitutive_name", "CMStrucLinIso" );
        aParameterLists.set( "model_type", fem::Model_Type::FULL );
        aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
        aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
        aParameterLists.set( "properties", "PropYoung,YoungsModulus;PropPoisson,PoissonRatio" );

        //------------------------------------------------------------------------------
        // STABILIZATION PARAMETERS

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoung,Material" );

        aParameterLists( FEM::STABILIZATION ).add_parameter_list();
        aParameterLists.set( "stabilization_name", "SPInterfaceNitsche" );
        aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
        aParameterLists.set( "function_parameters", "100.0" );
        aParameterLists.set( "leader_properties", "PropYoung,Material" );
        aParameterLists.set( "follower_properties", "PropYoung,Material" );

        if ( tUseGhost )
        {
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGhost" );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", "0.01" );
            aParameterLists.set( "leader_properties", "PropYoung,Material" );
        }

        //------------------------------------------------------------------------------
        // IWGs

        // bulk elasticity in every phase (no body force)
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGBulk" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        // symmetric Nitsche Dirichlet: the exact field on all six outer sides of every phase
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGDirichlet" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_properties", "PropExactDisp,Dirichlet" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPNitsche,DirichletNitsche" );
        aParameterLists.set( "mesh_set_names", tSideSets );

        // symmetric Nitsche interface, same CM on leader and follower
        aParameterLists( FEM::IWG ).add_parameter_list();
        aParameterLists.set( "IWG_name", "IWGInterface" );
        aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
        aParameterLists.set( "dof_residual", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "follower_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "stabilization_parameters", "SPInterfaceNitsche,NitscheInterface" );
        aParameterLists.set( "mesh_set_names", tInterfaceSets );

        if ( tUseGhost )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGhost" );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY,UZ" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "stabilization_parameters", "SPGhost,GhostSP" );
            aParameterLists.set( "mesh_set_names", tGhostSets );
        }

        //------------------------------------------------------------------------------
        // IQIs

        // squared L2 displacement error vs the exact field
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDispErrorL2" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_properties", "PropExactDisp,L2Check" );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIUX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIUY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 1 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIUZ" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::DOF );
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "vectorial_field_index", 2 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        // sigma_xx: NORMAL_STRESS, index 0 -> Voigt component 0
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIStressXX" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::NORMAL_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 0 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        // sigma_xy: SHEAR_STRESS, index 2 -> Voigt component 5 (order xx,yy,zz,yz,xz,xy)
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIStressXY" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::SHEAR_STRESS );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
        aParameterLists.set( "vectorial_field_index", 2 );
        aParameterLists.set( "mesh_set_names", tBulkSets );

        //------------------------------------------------------------------------------
        // COMPUTATION (no sensitivities)

        aParameterLists( FEM::COMPUTATION );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // 2 x 2 x 2 elements on [0,2]^3, one pattern, linear Lagrange and linear B-spline, no refinement
        aParameterLists.set( "number_of_elements_per_dimension", Vector< uint >{ 2, 2, 2 } );
        aParameterLists.set( "domain_dimensions", Vector< real >{ 2.0, 2.0, 2.0 } );
        aParameterLists.set( "domain_offset", Vector< real >{ 0.0, 0.0, 0.0 } );

        aParameterLists.set( "lagrange_output_meshes", "0" );

        aParameterLists.set( "lagrange_orders", "1" );
        aParameterLists.set( "lagrange_pattern", "0" );
        aParameterLists.set( "bspline_orders", "1" );
        aParameterLists.set( "bspline_pattern", "0" );

        aParameterLists.set( "lagrange_to_bspline", "0" );

        aParameterLists.set( "pattern_initial_refinement", Vector< uint >{ 0 } );
    }

    /* ------------------------------------------------------------------------ */

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // displacement DOFs interpolated on B-spline mesh 0
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
        aParameterLists.set( "UZ", 0 );

        aParameterLists.set( "order_adofs_by_host", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // all bulk sets of all four phases (phase 2 is empty)
        std::string tVisSets =
                "HMR_dummy_n_p0,HMR_dummy_c_p0,"
                "HMR_dummy_n_p1,HMR_dummy_c_p1,"
                "HMR_dummy_n_p2,HMR_dummy_c_p2,"
                "HMR_dummy_n_p3,HMR_dummy_c_p3";

        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::FULL_DISCONTINUOUS );
        aParameterLists.set( "Set_Names", tVisSets );
        aParameterLists.set( "Field_Names", "IQIUX,IQIUY,IQIUZ,IQIDispErrorL2,IQIStressXX,IQIStressXY" );
        aParameterLists.set( "Field_Type", "NODAL,NODAL,NODAL,GLOBAL,ELEMENTAL_AVG,ELEMENTAL_AVG" );
        aParameterLists.set( "IQI_Names", "IQIUX,IQIUY,IQIUZ,IQIDispErrorL2,IQIStressXX,IQIStressXY" );
        aParameterLists.set( "Save_Frequency", 1 );
    }

    /* ------------------------------------------------------------------------ */

    void
    MIGParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // default MIG list: no periodic boundary conditions
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // conformal cut (HEX8 regular subdivision + node hierarchy), enrichment on B-spline mesh 0
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );

        aParameterLists.set( "ghost_stab", tUseGhost );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );

        // keep the sets of empty phases registered, so FEM and VIS can name them
        aParameterLists.set( "deactivate_empty_sets", false );

        // low-to-high interface names only (dbl_iside_p0_<lo>_p1_<hi>), as used by FEM
        aParameterLists.set( "high_to_low_dbl_side_sets", false );
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( moris::tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // direct solve: Amesos with KLU set explicitly
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aParameterLists.set( "Solver_Type", "Amesos_Klu" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        // Newton: one linear solve plus one residual check
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_max_iter", 2 );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-12 );
        aParameterLists.set( "NLA_tot_res_norm", 1.0e-20 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY,UZ" );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );

        // single static step
        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 0 );
        aParameterLists.set( "TSA_Num_Time_Steps", 1 );
        aParameterLists.set( "TSA_Time_Frame", 1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_Solver_algorithms", "0" );
        aParameterLists.set( "TSA_DofTypes", "UX,UY,UZ" );
        aParameterLists.set( "TSA_Output_Indices", "0" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
