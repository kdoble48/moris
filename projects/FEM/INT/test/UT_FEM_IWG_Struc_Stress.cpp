/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IWG_Struc_Stress.cpp
 *
 */

#include <string>
#include <catch.hpp>
#include <memory>
#include "assert.hpp"

#define protected public
#define private   public
//FEM//INT//src
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IWG.hpp"
#include "cl_FEM_Set.hpp"
#undef protected
#undef private
//MTK/src
#include "cl_MTK_Enums.hpp"
//LINALG/src
#include "op_equal_equal.hpp"
//FEM//INT//src
#include "cl_FEM_Enums.hpp"
#include "cl_FEM_Field_Interpolator.hpp"
#include "cl_FEM_Property.hpp"
#include "cl_FEM_CM_Factory.hpp"
#include "cl_FEM_IWG_Factory.hpp"
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_Elasticity_UT.cpp"

using namespace moris;
using namespace fem;

// L2 projection of the von Mises stress, R = int v ( s* - sigma_vm( u ) ):
// Jacobian wrt s* (mass matrix) and wrt u (through fn_FEM_Stress_Measures) checked by finite differences
TEST_CASE( "IWG_Struc_Von_Mises_Stress", "[moris],[fem],[IWG_Struc_Stress]" )
{
    // define an epsilon environment
    real tEpsilon = 1E-6;

    // define a perturbation relative size
    real tPerturbation = 1E-6;

    // create list of interpolation orders
    Vector< mtk::Interpolation_Order > tInterpolationOrders = {
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Order::QUADRATIC };

    // create list of integration orders
    Vector< mtk::Integration_Order > tIntegrationOrders = {
            mtk::Integration_Order::QUAD_2x2,
            mtk::Integration_Order::HEX_2x2x2 };

    // create list with number of coeffs
    Matrix< DDRMat > tNumCoeffs = { { 8, 18 }, { 16, 54 } };

    // dof type list: displacement and the projected stress measure
    Vector< MSI::Dof_Type >           tStressDofTypes = { MSI::Dof_Type::STRESS_DOF };
    Vector< Vector< MSI::Dof_Type > > tDispDofTypes   = { { MSI::Dof_Type::UX } };
    Vector< Vector< MSI::Dof_Type > > tDofTypes       = { tDispDofTypes( 0 ), tStressDofTypes };

    // init IWG
    //------------------------------------------------------------------------------
    // create the properties
    std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );
    tPropEMod->set_val_function( tConstValFunc_Elast );

    std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );
    tPropNu->set_val_function( tConstValFunc_Elast );

    // define constitutive models
    fem::CM_Factory tCMFactory;

    std::shared_ptr< fem::Constitutive_Model > tCMLeaderStrucLinIso =
            tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCMLeaderStrucLinIso->set_dof_type_list( tDispDofTypes );
    tCMLeaderStrucLinIso->set_property( tPropEMod, "YoungsModulus" );
    tCMLeaderStrucLinIso->set_property( tPropNu, "PoissonRatio" );
    tCMLeaderStrucLinIso->set_local_properties();

    // define the IWGs
    fem::IWG_Factory tIWGFactory;

    std::shared_ptr< fem::IWG > tIWG =
            tIWGFactory.create_IWG( fem::IWG_Type::STRUC_VON_MISES_STRESS );
    tIWG->set_residual_dof_type( { tStressDofTypes } );
    tIWG->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIWG->set_constitutive_model( tCMLeaderStrucLinIso, "ElastLinIso" );

    // init set info
    //------------------------------------------------------------------------------
    // set a fem set pointer
    MSI::Equation_Set* tSet = new fem::Set();
    static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );
    tIWG->set_set_pointer( static_cast< fem::Set* >( tSet ) );

    // set size for the set EqnObjDofTypeList
    tIWG->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );

    // set size and populate the set dof type map
    tIWG->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIWG->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::STRESS_DOF ) ) = 1;

    // set size and populate the set leader dof type map
    tIWG->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) )         = 0;
    tIWG->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::STRESS_DOF ) ) = 1;

    // 2D plane stress (3 flux components), 2D plane strain (4), 3D (6)
    Vector< std::pair< uint, fem::Model_Type > > tCases = {
        { 2, fem::Model_Type::PLANE_STRESS },
        { 2, fem::Model_Type::PLANE_STRAIN },
        { 3, fem::Model_Type::PLANE_STRESS }
    };

    for ( const auto& [ iSpaceDim, tModelType ] : tCases )
    {
        // set geometry inputs
        mtk::Geometry_Type tGeometryType = iSpaceDim == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;

        // set space dimension to CM
        tCMLeaderStrucLinIso->set_model_type( tModelType );
        tCMLeaderStrucLinIso->set_space_dim( iSpaceDim );

        // loop on the interpolation order
        for ( uint iInterpOrder = 1; iInterpOrder < 3; iInterpOrder++ )
        {
            // create an interpolation order
            mtk::Interpolation_Order tInterpolationOrder = tInterpolationOrders( iInterpOrder - 1 );

            // space and time geometry interpolators
            //------------------------------------------------------------------------------
            mtk::Interpolation_Rule tGIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            Geometry_Interpolator tGI = Geometry_Interpolator( tGIRule );

            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };
            Matrix< DDRMat > tXHat;
            fill_xhat_Elast( tXHat, iSpaceDim, iInterpOrder );
            tGI.set_coeff( tXHat, tTHat );

            // integration points
            //------------------------------------------------------------------------------
            mtk::Integration_Rule tIntegrationRule(
                    tGeometryType,
                    mtk::Integration_Type::GAUSS,
                    tIntegrationOrders( iSpaceDim - 2 ),
                    mtk::Geometry_Type::LINE,
                    mtk::Integration_Type::GAUSS,
                    mtk::Integration_Order::BAR_1 );

            mtk::Integrator tIntegrator( tIntegrationRule );

            Matrix< DDRMat > tIntegPoints;
            tIntegrator.get_points( tIntegPoints );

            // field interpolators
            //------------------------------------------------------------------------------
            uint tNumCoeff   = tNumCoeffs( iSpaceDim - 2, iInterpOrder - 1 );
            int  tNumDofDisp = tNumCoeff * iSpaceDim;
            int  tNumDofSig  = tNumCoeff;

            mtk::Interpolation_Rule tFIRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    tInterpolationOrder,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            // fill coefficients for leader FIs
            Matrix< DDRMat > tLeaderDofHatDisp;
            fill_uhat_Elast( tLeaderDofHatDisp, iSpaceDim, iInterpOrder );
            Matrix< DDRMat > tLeaderDofHatSig;
            fill_phat_Elast( tLeaderDofHatSig, iSpaceDim, iInterpOrder );

            Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );
            tLeaderFIs( 0 ) = new Field_Interpolator( iSpaceDim, tFIRule, &tGI, tDispDofTypes( 0 ) );
            tLeaderFIs( 0 )->set_coeff( tLeaderDofHatDisp );
            tLeaderFIs( 1 ) = new Field_Interpolator( 1, tFIRule, &tGI, tStressDofTypes );
            tLeaderFIs( 1 )->set_coeff( tLeaderDofHatSig );

            // set size and fill the set residual assembly map
            tIWG->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDofDisp - 1 } };
            tIWG->mSet->mResDofAssemblyMap( 1 ) = { { tNumDofDisp, tNumDofDisp + tNumDofSig - 1 } };

            // set size and fill the set jacobian assembly map
            Matrix< DDSMat > tJacAssembly = {
                { 0, tNumDofDisp - 1 },
                { tNumDofDisp, tNumDofDisp + tNumDofSig - 1 }
            };
            tIWG->mSet->mJacDofAssemblyMap.resize( tDofTypes.size() );
            tIWG->mSet->mJacDofAssemblyMap( 0 ) = tJacAssembly;
            tIWG->mSet->mJacDofAssemblyMap( 1 ) = tJacAssembly;

            // set size and init the set residual and jacobian
            tIWG->mSet->mResidual.resize( 1 );
            tIWG->mSet->mResidual( 0 ).set_size( tNumDofDisp + tNumDofSig, 1, 0.0 );
            tIWG->mSet->mJacobian.set_size( tNumDofDisp + tNumDofSig, tNumDofDisp + tNumDofSig, 0.0 );

            // build global dof type list
            tIWG->get_global_dof_type_list();

            // populate the requested leader dof type
            tIWG->mRequestedLeaderGlobalDofTypes = tDofTypes;

            // create a field interpolator manager
            Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
            Vector< Vector< enum mtk::Field_Type > > tDummyField;
            Field_Interpolator_Manager               tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );

            // populate the field interpolator manager
            tFIManager.mFI                     = tLeaderFIs;
            tFIManager.mIPGeometryInterpolator = &tGI;
            tFIManager.mIGGeometryInterpolator = &tGI;

            // set the interpolator manager to the set and the IWG
            tIWG->mSet->mLeaderFIManager = &tFIManager;
            tIWG->set_field_interpolator_manager( &tFIManager );

            // loop over integration points
            uint tNumGPs = tIntegPoints.n_cols();
            for ( uint iGP = 0; iGP < tNumGPs; iGP++ )
            {
                // reset IWG evaluation flags
                tIWG->reset_eval_flags();

                // set integration point
                tIWG->mSet->mLeaderFIManager->set_space_time( tIntegPoints.get_column( iGP ) );

                // residual
                tIWG->mSet->mResidual( 0 ).fill( 0.0 );
                tIWG->compute_residual( 1.0 );

                // jacobian by FD
                tIWG->mSet->mJacobian.fill( 0.0 );

                Matrix< DDRMat > tJacobian;
                Matrix< DDRMat > tJacobianFD;

                bool tCheckJacobian = tIWG->check_jacobian(
                        tPerturbation,
                        tEpsilon,
                        1.0,
                        tJacobian,
                        tJacobianFD,
                        true );

                if ( !tCheckJacobian )
                {
                    std::cout << "Case: dim " << iSpaceDim << " model " << static_cast< uint >( tModelType )
                              << " order " << iInterpOrder << " iGP " << iGP << '\n';
                }

                REQUIRE( tCheckJacobian );
            }

            // clean up
            tLeaderFIs.clear();
        }
    }
}
