/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * Stress_Measures_Patch_Test.cpp
 *
 * A patch test for the stress OUTPUTS: every stress measure has a value you can work out by hand, so each one is
 * checked against a closed form (example_test_case.cpp). Each check guards a defect that once gave a wrong number:
 *
 *   output                              exact value                              defect it guards
 *   PRINCIPAL_STRESS 1..3               eigenvalues of sigma                     MB-046: det(sigma) used s_yz for s_zz
 *   VON_MISES_STRESS                    sqrt( 3/2 s:s ), s = deviator            (reference, never wrong)
 *   MAX_NORMAL_STRESS (ref 1, exp 1)    int sigma_xx dV                          MB-049: aborted for every component
 *   JUMP_TRACTION on the interfaces     0 (traction is continuous)               MB-048: returned the SUM of tractions
 *   cases 1 and 2 (below)               identical results                        MB-047: one CM named on both sides
 *
 * Setup. [0,2]^3, a 2 x 2 x 2 HEX8 background mesh (linear Lagrange and B-splines), cut by two parallel planes with
 * unit normal n = (1, 0.3, 0.2)/|.|: plane A through p = (0.95, 1, 1) and plane B a distance tDelta = 0.1 further.
 * Phase 0 (behind A) is material 1 (E1 = 1); phases 1 and 3 (the layer and the far side) are material 2 (E2 = 10);
 * nu = 0.3 in both. The material-2 | material-2 interface (plane B) names the SAME constitutive model on its leader
 * and follower side, as many decks do; plane A is a bimaterial interface.
 *
 * Exact solution (case 0). With A a constant displacement gradient (all components nonzero, so every shear stress is
 * nonzero), u1 = A x + b behind plane A and u2 = u1 + a ( n . ( x - p ) ) beyond it. u is continuous across A, and
 * traction continuity C2 eps( A + a n^T ) n = C1 eps( A ) n fixes a:
 *     ( mu2 I + ( lambda2 + mu2 ) n n^T ) a = ( C1 - C2 ) eps( A ) n
 *     a = ( r - ( lambda2 + mu2 ) / ( lambda2 + 2 mu2 ) ( n . r ) n ) / mu2,   r = ( C1 - C2 ) eps( A ) n
 * (Sherman-Morrison, |n| = 1). Both stresses are uniform, so every stress output is exact to round-off.
 *
 * Cases (gCaseIndex, set by the EXA runner); one defect per case, so an old build fails each where it is broken:
 *   0  exact data above: the closed-form checks of the stress measures and the traction jumps
 *   1  the same boundary displacement plus a smooth non-polynomial perturbation, so the discrete stress differs across
 *      the interfaces; one CM named on both sides of plane B
 *   2  as 1, but with a second, identical CM on the follower side of plane B
 *   3  exact data, MAX_NORMAL_STRESS only (before the MB-049 fix it aborted the run, which would hide every other check)
 * Cases 1 and 2 must give the same strain energy: a shared CM object must behave like two copies (MB-047; before the
 * fix the shared object's traction cache and field-interpolator manager made the residual inconsistent with the
 * Jacobian, and Newton converged to a different answer).
 *
 * The 2D principal-stress formula (MB-051) and the Zienkiewicz-Zhu principal/normal/shear branches are covered by the
 * unit tests in projects/FEM/INT/test/UT_FEM_IQI_Stress_3D.cpp.
 */

#include <string>
#include <cmath>
#include <tuple>
#include <vector>
#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_TSA_Time_Solver.hpp"
#include "parameters.hpp"
#include "fn_PRM_MIG_Parameters.hpp"

//---------------------------------------------------------------

// case index, set by the EXA runner (see header comment)
extern uint gCaseIndex;

//---------------------------------------------------------------

#ifdef __cplusplus
extern "C" {
#endif
//------------------------------------------------------------------------------
namespace moris
{
    /* ------------------------------------------------------------------------ */
    // problem data

    const real tDelta = 0.1;    // distance between the planes (thickness of the material-2 layer)
    const real tE1    = 1.0;    // Young's modulus of material 1 (phase 0)
    const real tE2    = 10.0;   // Young's modulus of material 2 (phases 1 and 3)
    const real tNu    = 0.3;

    // displacement gradient and offset of u1 = A x + b
    const real tA[ 3 ][ 3 ] = { { 1.0e-3, 2.0e-4, 0.0 }, { 1.0e-4, -5.0e-4, 3.0e-4 }, { 0.0, 2.0e-4, 4.0e-4 } };
    const real tB[ 3 ]      = { 1.0e-3, -2.0e-3, 5.0e-4 };

    // plane A: point and (unnormalized) normal
    const real tP[ 3 ]  = { 0.95, 1.0, 1.0 };
    const real tN0[ 3 ] = { 1.0, 0.3, 0.2 };

    std::string tOutputFileName = "Stress_Measures_Patch_Test_" + std::to_string( gCaseIndex ) + ".exo";

    /* ------------------------------------------------------------------------ */
    // exact solution helpers

    // unit normal of the planes
    void
    unit_normal( real aN[ 3 ] )
    {
        real tNorm = std::sqrt( tN0[ 0 ] * tN0[ 0 ] + tN0[ 1 ] * tN0[ 1 ] + tN0[ 2 ] * tN0[ 2 ] );
        for ( uint i = 0; i < 3; i++ )
        {
            aN[ i ] = tN0[ i ] / tNorm;
        }
    }

    // jump vector a of the displacement gradient across plane A (header comment)
    void
    jump_vector( real aJump[ 3 ] )
    {
        real tN[ 3 ];
        unit_normal( tN );

        real tLambda1 = tE1 * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) ), tMu1 = tE1 / ( 2.0 * ( 1.0 + tNu ) );
        real tLambda2 = tE2 * tNu / ( ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu ) ), tMu2 = tE2 / ( 2.0 * ( 1.0 + tNu ) );

        // r = ( C1 - C2 ) eps( A ) n, with C eps = lambda tr( eps ) I + 2 mu eps
        real tTrace = tA[ 0 ][ 0 ] + tA[ 1 ][ 1 ] + tA[ 2 ][ 2 ];
        real tR[ 3 ];
        for ( uint i = 0; i < 3; i++ )
        {
            real tEpsN = 0.0;
            for ( uint j = 0; j < 3; j++ )
            {
                tEpsN += 0.5 * ( tA[ i ][ j ] + tA[ j ][ i ] ) * tN[ j ];
            }
            tR[ i ] = ( tLambda1 - tLambda2 ) * tTrace * tN[ i ] + 2.0 * ( tMu1 - tMu2 ) * tEpsN;
        }

        real tNR = tR[ 0 ] * tN[ 0 ] + tR[ 1 ] * tN[ 1 ] + tR[ 2 ] * tN[ 2 ];
        for ( uint i = 0; i < 3; i++ )
        {
            aJump[ i ] = ( tR[ i ] - ( tLambda2 + tMu2 ) / ( tLambda2 + 2.0 * tMu2 ) * tNR * tN[ i ] ) / tMu2;
        }
    }

    // signed distance to plane A
    real
    distance_to_plane_A( const Matrix< DDRMat >& aX )
    {
        real tN[ 3 ];
        unit_normal( tN );
        return tN[ 0 ] * ( aX( 0 ) - tP[ 0 ] ) + tN[ 1 ] * ( aX( 1 ) - tP[ 1 ] ) + tN[ 2 ] * ( aX( 2 ) - tP[ 2 ] );
    }

    /* ------------------------------------------------------------------------ */
    // geometry: A = signed distance to plane A, B = A - tDelta

    real
    Geometry_A( const Matrix< DDRMat >& aCoordinates, const Vector< real >& aGeometryParameters )
    {
        return distance_to_plane_A( aCoordinates );
    }

    real
    Geometry_B( const Matrix< DDRMat >& aCoordinates, const Vector< real >& aGeometryParameters )
    {
        return distance_to_plane_A( aCoordinates ) - tDelta;
    }

    /* ------------------------------------------------------------------------ */
    // property functions

    void
    Func_Const( Matrix< DDRMat >&    aPropMatrix,
            Vector< Matrix< DDRMat > >& aParameters,
            fem::Field_Interpolator_Manager* )
    {
        aPropMatrix = aParameters( 0 );
    }

    // boundary displacement: the exact field (case 0), plus a smooth perturbation (cases 1 and 2)
    void
    Func_Boundary_Displacement( Matrix< DDRMat >& aPropMatrix,
            Vector< Matrix< DDRMat > >&           aParameters,
            fem::Field_Interpolator_Manager*      aFIManager )
    {
        const Matrix< DDRMat >& tX = aFIManager->get_IP_geometry_interpolator()->valx();

        real tJump[ 3 ];
        jump_vector( tJump );
        real tSide = distance_to_plane_A( tX );

        aPropMatrix.set_size( 3, 1, 0.0 );
        for ( uint i = 0; i < 3; i++ )
        {
            real tValue = tB[ i ];
            for ( uint j = 0; j < 3; j++ )
            {
                tValue += tA[ i ][ j ] * tX( j );
            }
            if ( tSide > 0.0 )
            {
                tValue += tJump[ i ] * tSide;
            }
            if ( gCaseIndex == 1 or gCaseIndex == 2 )
            {
                tValue += 1.0e-3 * std::sin( 1.3 * tX( 0 ) + 0.7 * tX( 1 ) * ( i + 1 ) ) * std::cos( 0.9 * tX( 2 ) );
            }
            aPropMatrix( i ) = tValue;
        }
    }

    /* ------------------------------------------------------------------------ */

    void
    OPTParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "is_optimization_problem", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    HMRParameterList( Module_Parameter_Lists& aParameterLists )
    {
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
    GENParameterList( Module_Parameter_Lists& aParameterLists )
    {
        // geometry B first: phase = 2 * sign( B ) + sign( A ), so phase 0 = behind A, 1 = the layer, 3 = beyond B
        for ( std::string tName : { "B", "A" } )
        {
            aParameterLists( GEN::GEOMETRIES ).add_parameter_list( prm::create_level_set_geometry_parameter_list( gen::Field_Type::USER_DEFINED ) );
            aParameterLists.set( "name", "Geom" + tName );
            aParameterLists.set( "field_function_name", "Geometry_" + tName );
            aParameterLists.set( "discretization_mesh_index", -1 );
            aParameterLists.set( "use_multilinear_interpolation", false );
            aParameterLists.set( "isocontour_threshold", 0.0 );
        }
    }

    /* ------------------------------------------------------------------------ */
    // set names: material 1 = phase 0, material 2 = phases 1, 2 (empty) and 3

    std::string
    phase_sets( const std::string& aPrefix, const std::vector< int >& aPhases )
    {
        std::string tSets;
        for ( int iP : aPhases )
        {
            tSets += ( tSets.empty() ? "" : "," ) + aPrefix + std::to_string( iP );
        }
        return tSets;
    }

    std::string
    bulk_sets( const std::vector< int >& aPhases )
    {
        return phase_sets( "HMR_dummy_n_p", aPhases ) + "," + phase_sets( "HMR_dummy_c_p", aPhases );
    }

    std::string
    side_sets( const std::vector< int >& aPhases )
    {
        std::string tSets;
        for ( int iS = 1; iS <= 6; iS++ )
        {
            for ( std::string tCut : { "n", "c" } )
            {
                tSets += ( tSets.empty() ? "" : "," ) + phase_sets( "SideSet_" + std::to_string( iS ) + "_" + tCut + "_p", aPhases );
            }
        }
        return tSets;
    }

    const std::string tBimaterialInterfaces = "dbl_iside_p0_0_p1_1,dbl_iside_p0_0_p1_2,dbl_iside_p0_0_p1_3";
    const std::string tMaterial2Interfaces  = "dbl_iside_p0_1_p1_2,dbl_iside_p0_1_p1_3,dbl_iside_p0_2_p1_3";

    /* ------------------------------------------------------------------------ */

    void
    FEMParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.hack_for_legacy_fem();

        // properties
        for ( auto [ tName, tValue ] : std::vector< std::pair< std::string, std::string > >{
                      { "PropYoung1", std::to_string( tE1 ) }, { "PropYoung2", std::to_string( tE2 ) }, { "PropPoisson", std::to_string( tNu ) },
                      { "PropOne", "1.0" }, { "PropZero", "0.0" } } )
        {
            aParameterLists( FEM::PROPERTIES ).add_parameter_list();
            aParameterLists.set( "property_name", tName );
            aParameterLists.set( "function_parameters", tValue );
            aParameterLists.set( "value_function", "Func_Const" );
        }
        aParameterLists( FEM::PROPERTIES ).add_parameter_list();
        aParameterLists.set( "property_name", "PropBoundaryDisp" );
        aParameterLists.set( "function_parameters", "0.0;0.0;0.0" );
        aParameterLists.set( "value_function", "Func_Boundary_Displacement" );

        // constitutive models; CMStrucLinIso2F is an identical copy of CMStrucLinIso2, used by case 2 only
        for ( auto [ tName, tYoung ] : std::vector< std::pair< std::string, std::string > >{
                      { "CMStrucLinIso1", "PropYoung1" }, { "CMStrucLinIso2", "PropYoung2" }, { "CMStrucLinIso2F", "PropYoung2" } } )
        {
            aParameterLists( FEM::CONSTITUTIVE_MODELS ).add_parameter_list();
            aParameterLists.set( "constitutive_name", tName );
            aParameterLists.set( "model_type", fem::Model_Type::FULL );
            aParameterLists.set( "constitutive_type", fem::Constitutive_Type::STRUC_LIN_ISO );
            aParameterLists.set( "dof_dependencies", std::pair< std::string, std::string >( "UX,UY,UZ", "Displacement" ) );
            aParameterLists.set( "properties", tYoung + ",YoungsModulus;PropPoisson,PoissonRatio" );
        }

        // follower CM of the material-2 | material-2 interface: the same name as the leader, except in case 2
        const std::string tFollowerCM2 = gCaseIndex == 2 ? "CMStrucLinIso2F" : "CMStrucLinIso2";

        // stabilization parameters
        for ( int iM : { 1, 2 } )
        {
            std::string tM = std::to_string( iM );
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPNitsche" + tM );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::DIRICHLET_NITSCHE );
            aParameterLists.set( "function_parameters", "100.0" );
            aParameterLists.set( "leader_properties", "PropYoung" + tM + ",Material" );

            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", "SPGhost" + tM );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::GHOST_DISPL );
            aParameterLists.set( "function_parameters", "0.01" );
            aParameterLists.set( "leader_properties", "PropYoung" + tM + ",Material" );
        }
        for ( auto [ tName, tFollowerYoung ] : std::vector< std::pair< std::string, std::string > >{ { "SPInterface12", "PropYoung2" }, { "SPInterface22", "PropYoung2" } } )
        {
            aParameterLists( FEM::STABILIZATION ).add_parameter_list();
            aParameterLists.set( "stabilization_name", tName );
            aParameterLists.set( "stabilization_type", fem::Stabilization_Type::NITSCHE_INTERFACE );
            aParameterLists.set( "function_parameters", "100.0" );
            aParameterLists.set( "leader_properties", std::string( tName == "SPInterface12" ? "PropYoung1" : "PropYoung2" ) + ",Material" );
            aParameterLists.set( "follower_properties", tFollowerYoung + ",Material" );
        }

        // IWGs: bulk, Dirichlet (symmetric Nitsche) and ghost per material
        for ( int iM : { 1, 2 } )
        {
            std::string tM = std::to_string( iM );
            const std::vector< int > tPhases = iM == 1 ? std::vector< int >{ 0 } : std::vector< int >{ 1, 2, 3 };

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGBulk" + tM );
            aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_BULK );
            aParameterLists.set( "dof_residual", "UX,UY,UZ" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso" + tM + ",ElastLinIso" );
            aParameterLists.set( "mesh_set_names", bulk_sets( tPhases ) );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGDirichlet" + tM );
            aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_DIRICHLET_SYMMETRIC_NITSCHE );
            aParameterLists.set( "dof_residual", "UX,UY,UZ" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "leader_properties", "PropBoundaryDisp,Dirichlet" );
            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso" + tM + ",ElastLinIso" );
            aParameterLists.set( "stabilization_parameters", "SPNitsche" + tM + ",DirichletNitsche" );
            aParameterLists.set( "mesh_set_names", side_sets( tPhases ) );

            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", "IWGGhost" + tM );
            aParameterLists.set( "IWG_type", fem::IWG_Type::GHOST_NORMAL_FIELD );
            aParameterLists.set( "dof_residual", "UX,UY,UZ" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "stabilization_parameters", "SPGhost" + tM + ",GhostSP" );
            aParameterLists.set( "mesh_set_names", phase_sets( "ghost_p", tPhases ) );
        }

        // IWGs: interfaces (symmetric Nitsche); plane B names one CM on both sides (cases 0 and 1)
        for ( auto [ tName, tLeaderCM, tFollowerCM, tSP, tSets ] : std::vector< std::tuple< std::string, std::string, std::string, std::string, std::string > >{
                      { "IWGInterface12", "CMStrucLinIso1", "CMStrucLinIso2", "SPInterface12", tBimaterialInterfaces },
                      { "IWGInterface22", "CMStrucLinIso2", tFollowerCM2, "SPInterface22", tMaterial2Interfaces } } )
        {
            aParameterLists( FEM::IWG ).add_parameter_list();
            aParameterLists.set( "IWG_name", tName );
            aParameterLists.set( "IWG_type", fem::IWG_Type::STRUC_LINEAR_INTERFACE_SYMMETRIC_NITSCHE );
            aParameterLists.set( "dof_residual", "UX,UY,UZ" );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "leader_constitutive_models", tLeaderCM + ",ElastLinIso" );
            aParameterLists.set( "follower_constitutive_models", tFollowerCM + ",ElastLinIso" );
            aParameterLists.set( "stabilization_parameters", tSP + ",NitscheInterface" );
            aParameterLists.set( "mesh_set_names", tSets );
        }

        // IQIs: displacement error, and per material: the stress measures (ELEMENTAL_AVG), MAX_NORMAL_STRESS with
        // reference 1, exponent 1 and shift 0 (= int sigma_xx dV), volume and strain energy (GLOBAL)
        aParameterLists( FEM::IQI ).add_parameter_list();
        aParameterLists.set( "IQI_name", "IQIDispErrorL2" );
        aParameterLists.set( "IQI_type", fem::IQI_Type::L2_ERROR_ANALYTIC );
        aParameterLists.set( "dof_quantity", "UX,UY,UZ" );
        aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
        aParameterLists.set( "leader_properties", "PropBoundaryDisp,L2Check" );
        aParameterLists.set( "mesh_set_names", bulk_sets( { 0, 1, 2, 3 } ) );

        for ( int iM : { 1, 2 } )
        {
            std::string tM      = std::to_string( iM );
            std::string tCM     = "CMStrucLinIso" + tM + ",ElastLinIso";
            std::string tSets   = iM == 1 ? bulk_sets( { 0 } ) : bulk_sets( { 1, 2, 3 } );

            for ( auto [ tName, tType, tIndex ] : std::vector< std::tuple< std::string, fem::IQI_Type, sint > >{
                          { "SXX", fem::IQI_Type::NORMAL_STRESS, 0 },
                          { "SVM", fem::IQI_Type::VON_MISES_STRESS, -1 },
                          { "SP1", fem::IQI_Type::PRINCIPAL_STRESS, 0 },
                          { "SP2", fem::IQI_Type::PRINCIPAL_STRESS, 1 },
                          { "SP3", fem::IQI_Type::PRINCIPAL_STRESS, 2 } } )
            {
                aParameterLists( FEM::IQI ).add_parameter_list();
                aParameterLists.set( "IQI_name", tName + "_M" + tM );
                aParameterLists.set( "IQI_type", tType );
                aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
                aParameterLists.set( "leader_constitutive_models", tCM );
                if ( tIndex >= 0 )
                {
                    aParameterLists.set( "vectorial_field_index", tIndex );
                }
                aParameterLists.set( "mesh_set_names", tSets );
            }

            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "MAXSXX_M" + tM );
            aParameterLists.set( "IQI_type", fem::IQI_Type::MAX_NORMAL_STRESS );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "leader_constitutive_models", tCM );
            aParameterLists.set( "leader_properties", "PropOne,ReferenceValue;PropOne,Exponent;PropZero,Shift" );
            aParameterLists.set( "vectorial_field_index", 0 );
            aParameterLists.set( "mesh_set_names", tSets );

            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "VOL_M" + tM );
            aParameterLists.set( "IQI_type", fem::IQI_Type::VOLUME );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "mesh_set_names", tSets );

            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", "SE_M" + tM );
            aParameterLists.set( "IQI_type", fem::IQI_Type::STRAIN_ENERGY );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "leader_constitutive_models", "CMStrucLinIso" + tM + ",Elast" );
            aParameterLists.set( "mesh_set_names", tSets );
        }

        // traction jumps: t_leader( n ) - t_follower( n ) with the leader normal n; zero in an exact solution
        for ( auto [ tName, tLeaderCM, tFollowerCM, tSets ] : std::vector< std::tuple< std::string, std::string, std::string, std::string > >{
                      { "JT_BI", "CMStrucLinIso1", "CMStrucLinIso2", tBimaterialInterfaces },
                      { "JT_M2", "CMStrucLinIso2", tFollowerCM2, tMaterial2Interfaces } } )
        {
            aParameterLists( FEM::IQI ).add_parameter_list();
            aParameterLists.set( "IQI_name", tName );
            aParameterLists.set( "IQI_type", fem::IQI_Type::JUMP_TRACTION );
            aParameterLists.set( "leader_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "follower_dof_dependencies", "UX,UY,UZ" );
            aParameterLists.set( "leader_constitutive_models", tLeaderCM + ",TractionCM" );
            aParameterLists.set( "follower_constitutive_models", tFollowerCM + ",TractionCM" );
            aParameterLists.set( "mesh_set_names", tSets );
        }

        aParameterLists( FEM::COMPUTATION );
    }

    /* ------------------------------------------------------------------------ */

    void
    SOLParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists( SOL::LINEAR_ALGORITHMS ).add_parameter_list( sol::SolverType::AMESOS_IMPL );
        aParameterLists.set( "Solver_Type", "Amesos_Klu" );

        aParameterLists( SOL::LINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "DLA_Linear_solver_algorithms", "0" );

        // a linear problem: one Newton step solves it when the residual and the Jacobian agree. Cases 1 and 2 allow
        // more iterations, so a residual inconsistent with its Jacobian (MB-047) converges to its own, different root.
        aParameterLists( SOL::NONLINEAR_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "NLA_Linear_solver", 0 );
        aParameterLists.set( "NLA_max_iter", gCaseIndex == 1 or gCaseIndex == 2 ? 10 : 2 );
        aParameterLists.set( "NLA_rel_res_norm_drop", 1.0e-12 );
        aParameterLists.set( "NLA_tot_res_norm", 1.0e-20 );
        aParameterLists.set( "NLA_combined_res_jac_assembly", false );

        aParameterLists( SOL::NONLINEAR_SOLVERS ).add_parameter_list();
        aParameterLists.set( "NLA_DofTypes", "UX,UY,UZ" );
        aParameterLists.set( "NLA_Nonlinear_solver_algorithms", "0" );

        aParameterLists( SOL::TIME_SOLVER_ALGORITHMS ).add_parameter_list();
        aParameterLists.set( "TSA_Nonlinear_Solver", 0 );
        aParameterLists.set( "TSA_Num_Time_Steps", 1 );
        aParameterLists.set( "TSA_Time_Frame", 1.0 );

        aParameterLists( SOL::TIME_SOLVERS ).add_parameter_list();
        aParameterLists.set( "TSA_Solver_algorithms", "0" );
        aParameterLists.set( "TSA_DofTypes", "UX,UY,UZ" );
        aParameterLists.set( "TSA_Output_Indices", "0,1" );
        aParameterLists.set( "TSA_Output_Criteria", "Output_Criterion,Output_Criterion" );

        aParameterLists( SOL::PRECONDITIONERS ).add_parameter_list( sol::PreconditionerType::NONE );
    }

    /* ------------------------------------------------------------------------ */

    void
    MSIParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "UX", 0 );
        aParameterLists.set( "UY", 0 );
        aParameterLists.set( "UZ", 0 );
        aParameterLists.set( "order_adofs_by_host", false );
    }

    /* ------------------------------------------------------------------------ */

    void
    VISParameterList( Module_Parameter_Lists& aParameterLists )
    {
        const std::string tVisSets = bulk_sets( { 0, 1, 2, 3 } );

        // output 0: bulk fields per material (each IQI lives on its own material's sets)
        std::string tFields;
        std::string tTypes;
        for ( std::string tM : { "1", "2" } )
        {
            for ( std::string tS : { "SXX", "SVM", "SP1", "SP2", "SP3" } )
            {
                tFields += tS + "_M" + tM + ",";
                tTypes += "ELEMENTAL_AVG,";
            }
            for ( std::string tG : { "MAXSXX", "VOL", "SE" } )
            {
                if ( tG == "MAXSXX" and gCaseIndex != 3 )
                {
                    continue;    // only case 3 evaluates MAX_NORMAL_STRESS (header comment)
                }
                tFields += tG + "_M" + tM + ",";
                tTypes += "GLOBAL,";
            }
        }
        tFields += "IQIDispErrorL2";
        tTypes += "GLOBAL";

        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", tOutputFileName ) );
        aParameterLists.set( "Temp_Name", std::pair< std::string, std::string >( "./", "temp_" + std::to_string( gCaseIndex ) + "_0.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::FULL_DISCONTINUOUS );
        aParameterLists.set( "Set_Names", tVisSets );
        aParameterLists.set( "Field_Names", tFields );
        aParameterLists.set( "Field_Type", tTypes );
        aParameterLists.set( "IQI_Names", tFields );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Output_Index", 0 );

        // output 1: the traction jumps on the double-sided interface sets (faceted fields need a STANDARD mesh)
        aParameterLists( 0 ).add_parameter_list( prm::create_vis_parameter_list() );
        aParameterLists.set( "File_Name", std::pair< std::string, std::string >( "./", "Stress_Measures_Patch_Test_" + std::to_string( gCaseIndex ) + "_faceted.exo" ) );
        aParameterLists.set( "Temp_Name", std::pair< std::string, std::string >( "./", "temp_" + std::to_string( gCaseIndex ) + "_1.exo" ) );
        aParameterLists.set( "Mesh_Type", vis::VIS_Mesh_Type::STANDARD );
        aParameterLists.set( "Set_Names", tVisSets + "," + tBimaterialInterfaces + "," + tMaterial2Interfaces );
        aParameterLists.set( "Field_Names", "JT_BI,JT_M2" );
        aParameterLists.set( "Field_Type", "FACETED_AVG,FACETED_AVG" );
        aParameterLists.set( "IQI_Names", "JT_BI,JT_M2" );
        aParameterLists.set( "Save_Frequency", 1 );
        aParameterLists.set( "Output_Index", 1 );
    }

    /* ------------------------------------------------------------------------ */

    void
    XTKParameterList( Module_Parameter_Lists& aParameterLists )
    {
        aParameterLists.set( "decompose", true );
        aParameterLists.set( "decomposition_type", "conformal" );
        aParameterLists.set( "enrich", true );
        aParameterLists.set( "basis_rank", "bspline" );
        aParameterLists.set( "enrich_mesh_indices", "0" );
        aParameterLists.set( "ghost_stab", true );
        aParameterLists.set( "multigrid", false );
        aParameterLists.set( "verbose", false );
    }

    void
    MIGParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    void
    MORISGENERALParameterList( Module_Parameter_Lists& aParameterLists )
    {
    }

    /* ------------------------------------------------------------------------ */

    bool
    Output_Criterion( tsa::Time_Solver* aTimeSolver )
    {
        return true;
    }

    /* ------------------------------------------------------------------------ */
}    // namespace moris

//------------------------------------------------------------------------------
#ifdef __cplusplus
}
#endif
