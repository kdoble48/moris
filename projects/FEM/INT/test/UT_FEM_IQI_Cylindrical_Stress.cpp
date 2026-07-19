/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Cylindrical_Stress.cpp
 *
 */

#include <cmath>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#define protected public
#define private public
// FEM/INT/src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_IQI.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private

// FEM/INT/src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IQI_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"
// MTK/src
#include "cl_MTK_Enums.hpp"

using namespace moris;
using namespace fem;

//------------------------------------------------------------------------------
// This test builds a 2D plane-stress linear-isotropic-elasticity constitutive
// model that drives a known Cartesian stress state, constructs the new
// IQI_Cylindrical_Stress through the factory, and verifies that the cylindrical
// (polar) stress components ( sigma_rr, sigma_thetatheta, sigma_rtheta ) returned
// by the IQI match the closed-form rotation of the Cartesian stress about a
// user-specified center point, for several polar angles.
//------------------------------------------------------------------------------

TEST_CASE( "IQI_Cylindrical_Stress", "[moris],[fem],[IQI_Cylindrical_Stress]" )
{
    // tolerance for floating point comparisons
    real tEpsilon = 1.0E-10;

    // displacement dof types ( vector field )
    Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes     = tDispDofTypes;

    // create the elasticity properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    // define the constitutive model ( 2D plane stress linear isotropic elasticity )
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderElastLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderElastLinIso->set_dof_type_list( { tDispDofTypes } );
    tCMLeaderElastLinIso->set_property( tPropEMod, "YoungsModulus" );
    tCMLeaderElastLinIso->set_property( tPropNu, "PoissonRatio" );
    tCMLeaderElastLinIso->set_local_properties();

    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );
    tCMLeaderElastLinIso->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set unique dof type list
    tCMLeaderElastLinIso->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tCMLeaderElastLinIso->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderElastLinIso->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // set size and populate the set leader dof type map
    tCMLeaderElastLinIso->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tCMLeaderElastLinIso->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

    // build global dof type list
    tCMLeaderElastLinIso->get_global_dof_type_list();

    // set the model dimension and type
    tCMLeaderElastLinIso->set_space_dim( 2 );
    tCMLeaderElastLinIso->set_model_type( fem::Model_Type::PLANE_STRESS );

    // space and time geometry interpolators
    //------------------------------------------------------------------------------
    // QUAD geometry on the unit square so the element-center evaluation point is ( 0.5, 0.5 )
    Matrix< DDRMat > tXHat = { { 0.0, 0.0 },
        { 1.0, 0.0 },
        { 1.0, 1.0 },
        { 0.0, 1.0 } };

    // create a space geometry interpolation rule
    mtk::Interpolation_Rule tGIRule(
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // create a space time geometry interpolator
    Geometry_Interpolator tGI( tGIRule );

    // time coefficients
    Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

    // set the coefficients xHat, tHat
    tGI.set_coeff( tXHat, tTHat );

    // field interpolators
    //------------------------------------------------------------------------------
    // create a space time interpolation rule
    mtk::Interpolation_Rule tFIRule(
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );

    // fill displacement dof coefficients ( produces a non-trivial stress state )
    Matrix< DDRMat > tLeaderDOFHat;
    fill_uhat_Elast( tLeaderDOFHat, 2, 1 );

    // create the field interpolator for the displacement vector field
    Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );
    tLeaderFIs( 0 ) = new Field_Interpolator( 2, tFIRule, &tGI, tDispDofTypes( 0 ) );
    tLeaderFIs( 0 )->set_coeff( tLeaderDOFHat );

    // create a field interpolator manager
    Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
    Vector< Vector< enum mtk::Field_Type > > tDummyField;
    Field_Interpolator_Manager               tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

    // populate the field interpolator manager
    tFIManager.mFI                     = tLeaderFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;

    // set the interpolator manager to the set and the CM
    tCMLeaderElastLinIso->mSet->mLeaderFIManager = &tFIManager;
    tCMLeaderElastLinIso->set_field_interpolator_manager( &tFIManager );

    // create the IQI through the factory
    //------------------------------------------------------------------------------
    fem::IQI_Factory tIQIFactory;

    std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::CYLINDRICAL_STRESS );
    tIQI->set_constitutive_model( tCMLeaderElastLinIso, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tIQI->set_name( "Cylindrical Stress" );
    tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tIQI->set_field_interpolator_manager( &tFIManager );

    // evaluate everything at the element center ( xi = eta = 0, tau = 0 )
    Matrix< DDRMat > tParamPoint = { { 0.0 }, { 0.0 }, { 0.0 } };

    // reset CM flags and set the evaluation point
    tCMLeaderElastLinIso->reset_eval_flags();
    tFIManager.set_space_time( tParamPoint );

    // physical coordinates of the evaluation point ( unit square center )
    const Matrix< DDRMat >& tX = tGI.valx();
    real                    tx = tX( 0 );
    real                    ty = tX( 1 );
    CHECK( std::abs( tx - 0.5 ) < tEpsilon );
    CHECK( std::abs( ty - 0.5 ) < tEpsilon );

    // reference Cartesian stress from the constitutive model ( plane stress: [ sxx, syy, sxy ] )
    const Matrix< DDRMat >& tCMStress = tCMLeaderElastLinIso->flux();
    REQUIRE( tCMStress.numel() == 3 );
    real tSxx = tCMStress( 0 );
    real tSyy = tCMStress( 1 );
    real tSxy = tCMStress( 2 );

    // make sure the stress state is non-trivial so the rotation is meaningful
    REQUIRE( std::abs( tSxy ) > tEpsilon );
    REQUIRE( std::abs( tSxx - tSyy ) > tEpsilon );

    // storage for the IQI output
    Matrix< DDRMat > tQI;

    //------------------------------------------------------------------------------
    // theta = 0 : center to the -x side of the evaluation point ( c = 1, s = 0 )
    //   sigma_rr = sxx, sigma_thetatheta = syy, sigma_rtheta = sxy
    //------------------------------------------------------------------------------
    tIQI->set_parameters( { { { tx - 0.5 }, { ty } } } );

    tIQI->set_output_type_index( 0 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tSxx ) < tEpsilon );

    tIQI->set_output_type_index( 1 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tSyy ) < tEpsilon );

    tIQI->set_output_type_index( 2 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tSxy ) < tEpsilon );

    //------------------------------------------------------------------------------
    // theta = 90 deg : center directly below the evaluation point ( c = 0, s = 1 )
    //   sigma_rr = syy, sigma_thetatheta = sxx, sigma_rtheta = -sxy
    //------------------------------------------------------------------------------
    tIQI->set_parameters( { { { tx }, { ty - 0.5 } } } );

    tIQI->set_output_type_index( 0 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tSyy ) < tEpsilon );

    tIQI->set_output_type_index( 1 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tSxx ) < tEpsilon );

    tIQI->set_output_type_index( 2 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) + tSxy ) < tEpsilon );

    //------------------------------------------------------------------------------
    // theta = 45 deg : center at the origin ( c = s = sqrt(2)/2 ), exercising cross terms
    //   sigma_rr       = 0.5*( sxx + syy ) + sxy
    //   sigma_theta    = 0.5*( sxx + syy ) - sxy
    //   sigma_rtheta   = 0.5*( syy - sxx )
    //------------------------------------------------------------------------------
    tIQI->set_parameters( { { { 0.0 }, { 0.0 } } } );

    real tExpectedRR    = 0.5 * ( tSxx + tSyy ) + tSxy;
    real tExpectedTheta = 0.5 * ( tSxx + tSyy ) - tSxy;
    real tExpectedRT    = 0.5 * ( tSyy - tSxx );

    tIQI->set_output_type_index( 0 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tExpectedRR ) < tEpsilon );

    tIQI->set_output_type_index( 1 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tExpectedTheta ) < tEpsilon );

    tIQI->set_output_type_index( 2 );
    tIQI->compute_QI( tQI );
    CHECK( std::abs( tQI( 0 ) - tExpectedRT ) < tEpsilon );

    //------------------------------------------------------------------------------
    // invariant check: sigma_rr + sigma_thetatheta == sxx + syy ( trace preserved )
    // at an arbitrary center / angle
    //------------------------------------------------------------------------------
    tIQI->set_parameters( { { { 0.123 }, { -0.456 } } } );

    tIQI->set_output_type_index( 0 );
    tIQI->compute_QI( tQI );
    real tRR = tQI( 0 );

    tIQI->set_output_type_index( 1 );
    tIQI->compute_QI( tQI );
    real tTT = tQI( 0 );

    CHECK( std::abs( ( tRR + tTT ) - ( tSxx + tSyy ) ) < tEpsilon );

    // note: the field interpolator manager owns and deletes the field interpolators on scope exit;
    // the fem set is intentionally left for process teardown, matching the other FEM/INT unit tests

} /*END_TEST_CASE*/
