/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Max_Stress.cpp
 *
 */

#include <cmath>
#include <string>
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
#include "FEM_Test_Proxy/cl_FEM_Inputs_for_NS_Incompressible_UT.cpp"
// MTK/src
#include "cl_MTK_Enums.hpp"
#include "cl_MTK_Integrator.hpp"

using namespace moris;
using namespace fem;

namespace
{
    //------------------------------------------------------------------------------
    // von Mises stress from a CM flux vector (plane stress: 3, plane strain: 4, 3D: 6 components),
    // written independently of IQI_Max_Stress
    real
    reference_von_mises( const Matrix< DDRMat >& aFlux )
    {
        real tS11 = 0.0, tS22 = 0.0, tS33 = 0.0, tS23 = 0.0, tS13 = 0.0, tS12 = 0.0;
        if ( aFlux.numel() == 3 )
        {
            tS11 = aFlux( 0 );
            tS22 = aFlux( 1 );
            tS12 = aFlux( 2 );
        }
        else if ( aFlux.numel() == 4 )
        {
            tS11 = aFlux( 0 );
            tS22 = aFlux( 1 );
            tS33 = aFlux( 2 );
            tS12 = aFlux( 3 );
        }
        else
        {
            tS11 = aFlux( 0 );
            tS22 = aFlux( 1 );
            tS33 = aFlux( 2 );
            tS23 = aFlux( 3 );
            tS13 = aFlux( 4 );
            tS12 = aFlux( 5 );
        }
        return std::sqrt( 0.5 * ( ( tS11 - tS22 ) * ( tS11 - tS22 ) + ( tS22 - tS33 ) * ( tS22 - tS33 ) + ( tS33 - tS11 ) * ( tS33 - tS11 ) )
                          + 3.0 * ( tS23 * tS23 + tS13 * tS13 + tS12 * tS12 ) );
    }

    //------------------------------------------------------------------------------
    // one aggregation case: P-norm integrand ( sigma / sigma_ref )^p or KS integrand exp( rho ( sigma / sigma_ref - 1 ) )
    struct Aggregation_Case
    {
        std::string mName;
        bool        mKS;
        real        mParameter;    // p (P-norm) or rho (KS)
    };
}    // namespace

//------------------------------------------------------------------------------
// IQI_Max_Stress (MAX_VON_MISES_STRESS), P-norm and Kreisselmeier-Steinhauser (KSParameter) integrands:
//   (1) the pointwise value matches the closed form evaluated from the CM flux;
//   (2) the analytic dQI/du matches the finite-difference dQI/du (check_dQIdu_FD), for plane stress,
//       plane strain and 3D, linear and quadratic interpolation, at every Gauss point.
//------------------------------------------------------------------------------
TEST_CASE( "IQI_Max_Stress P-norm and KS", "[moris],[fem],[IQI_Max_Stress],[KS]" )
{
    real tPerturbation = 1.0E-6;
    real tEpsilon      = 1.0E-5;
    real tValueTol     = 1.0E-12;

    Vector< Aggregation_Case > tCases = {
        { "pnorm_p8", false, 8.0 },
        { "ks_rho20", true, 20.0 },
        { "ks_rho60", true, 60.0 }
    };

    // ( space dim, model type ) combinations: plane stress, plane strain, 3D
    Vector< uint >             tSpaceDims  = { 2, 2, 3 };
    Vector< fem::Model_Type >  tModelTypes = { fem::Model_Type::PLANE_STRESS, fem::Model_Type::PLANE_STRAIN, fem::Model_Type::PLANE_STRESS };
    Vector< std::string >      tModelNames = { "plane stress", "plane strain", "3D" };
    Matrix< DDRMat >           tNumCoeffs  = { { 8, 18 }, { 16, 54 } };

    Vector< mtk::Interpolation_Order > tInterpolationOrders = { mtk::Interpolation_Order::LINEAR, mtk::Interpolation_Order::QUADRATIC };

    for ( uint iCase = 0; iCase < tCases.size(); iCase++ )
    {
        const Aggregation_Case& tCase = tCases( iCase );

        for ( uint iModel = 0; iModel < tSpaceDims.size(); iModel++ )
        {
            uint tSpaceDim = tSpaceDims( iModel );

            mtk::Geometry_Type                tGeometryType = tSpaceDim == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;
            Vector< Vector< MSI::Dof_Type > > tDofTypes     = tSpaceDim == 2
                                                                    ? Vector< Vector< MSI::Dof_Type > >( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } )
                                                                    : Vector< Vector< MSI::Dof_Type > >( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } } );

            for ( uint iOrder = 1; iOrder < 3; iOrder++ )
            {
                // properties
                std::shared_ptr< fem::Property > tPropEMod = std::make_shared< fem::Property >();
                tPropEMod->set_parameters( { { { 1.0 } } } );
                std::shared_ptr< fem::Property > tPropNu = std::make_shared< fem::Property >();
                tPropNu->set_parameters( { { { 0.3 } } } );

                // elasticity CM
                fem::CM_Factory                            tCMFactory;
                std::shared_ptr< fem::Constitutive_Model > tCM = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
                tCM->set_dof_type_list( tDofTypes );
                tCM->set_property( tPropEMod, "YoungsModulus" );
                tCM->set_property( tPropNu, "PoissonRatio" );
                tCM->set_local_properties();
                // model type before space dim: set_space_dim sizes the strain/stress storage for the stored plane type
                // (cf. cl_FEM_Model_Initializer_Phasebased.cpp:202, "must come before set_space_dim")
                tCM->set_model_type( tModelTypes( iModel ) );
                tCM->set_space_dim( tSpaceDim );

                // set
                MSI::Equation_Set* tSet = new fem::Set();
                static_cast< fem::Set* >( tSet )->set_set_type( fem::Element_Type::BULK );

                // IQI (reference value is set below, once the stress level of the field is known)
                fem::IQI_Factory           tIQIFactory;
                std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::MAX_VON_MISES_STRESS );
                tIQI->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
                tIQI->set_constitutive_model( tCM, "ElastLinIso" );
                tIQI->set_name( "StressAgg" );
                tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );

                tIQI->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );
                tIQI->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
                tIQI->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
                tIQI->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
                tIQI->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

                // geometry interpolator
                mtk::Interpolation_Order tOrder = tInterpolationOrders( iOrder - 1 );
                mtk::Interpolation_Rule  tGIRule( tGeometryType, mtk::Interpolation_Type::LAGRANGE, tOrder, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );
                Geometry_Interpolator    tGI( tGIRule );
                Matrix< DDRMat >         tXHat;
                fill_xhat( tXHat, tSpaceDim, iOrder );
                Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };
                tGI.set_coeff( tXHat, tTHat );

                // integration points
                mtk::Integration_Rule tIntegrationRule(
                        tGeometryType,
                        mtk::Integration_Type::GAUSS,
                        tSpaceDim == 2 ? mtk::Integration_Order::QUAD_2x2 : mtk::Integration_Order::HEX_2x2x2,
                        mtk::Geometry_Type::LINE,
                        mtk::Integration_Type::GAUSS,
                        mtk::Integration_Order::BAR_1 );
                mtk::Integrator  tIntegrator( tIntegrationRule );
                Matrix< DDRMat > tIntegPoints;
                tIntegrator.get_points( tIntegPoints );

                // displacement field interpolator
                mtk::Interpolation_Rule       tFIRule( tGeometryType, mtk::Interpolation_Type::LAGRANGE, tOrder, mtk::Interpolation_Type::LAGRANGE, mtk::Interpolation_Order::LINEAR );
                Vector< Field_Interpolator* > tLeaderFIs( tDofTypes.size() );
                tLeaderFIs( 0 ) = new Field_Interpolator( tSpaceDim, tFIRule, &tGI, tDofTypes( 0 ) );
                Matrix< DDRMat > tUHat;
                fill_uhat_Elast( tUHat, tSpaceDim, iOrder );
                tLeaderFIs( 0 )->set_coeff( tUHat );

                int tNumDof = (int)tNumCoeffs( tSpaceDim - 2, iOrder - 1 ) * (int)tSpaceDim;
                tIQI->mSet->mResDofAssemblyMap.resize( tDofTypes.size() );
                tIQI->mSet->mResDofAssemblyMap( 0 ) = { { 0, tNumDof - 1 } };
                tIQI->mSet->mResidual.resize( 1 );
                tIQI->mSet->mResidual( 0 ).set_size( tNumDof, 1, 0.0 );

                moris::map< std::string, moris_index > tRequestedIQINamesAssemblyMap;
                tRequestedIQINamesAssemblyMap[ "StressAgg" ] = 0;
                tIQI->mSet->mRequestedIQINamesAssemblyMap   = tRequestedIQINamesAssemblyMap;
                tIQI->mSet->mQI.resize( 1 );
                tIQI->mSet->mQI( 0 ).set_size( 1, 1, 0.0 );

                tIQI->get_global_dof_type_list();
                tIQI->mRequestedLeaderGlobalDofTypes = tDofTypes;

                Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
                Vector< Vector< enum mtk::Field_Type > > tDummyField;
                Field_Interpolator_Manager               tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
                tFIManager.mFI                     = tLeaderFIs;
                tFIManager.mIPGeometryInterpolator = &tGI;
                tFIManager.mIGGeometryInterpolator = &tGI;
                tIQI->mSet->mLeaderFIManager       = &tFIManager;
                tIQI->set_field_interpolator_manager( &tFIManager );

                // reference stress = von Mises stress at the first Gauss point, so the integrand is O(1) and the KS
                // exponent stays within a few units of 0 over the element (a realistic, non-overflowing regime)
                tCM->reset_eval_flags();
                tFIManager.set_space_time( tIntegPoints.get_column( 0 ) );
                real tRef = reference_von_mises( tCM->flux() );
                REQUIRE( tRef > 0.0 );

                std::shared_ptr< fem::Property > tPropRef = std::make_shared< fem::Property >();
                tPropRef->set_parameters( { { { tRef } } } );
                tIQI->set_property( tPropRef, "ReferenceValue" );
                std::shared_ptr< fem::Property > tPropAgg = std::make_shared< fem::Property >();
                tPropAgg->set_parameters( { { { tCase.mParameter } } } );
                if ( tCase.mKS )
                {
                    tIQI->set_property( tPropAgg, "KSParameter" );    // shift defaults to 1.0
                }
                else
                {
                    std::shared_ptr< fem::Property > tPropShift = std::make_shared< fem::Property >();
                    tPropShift->set_parameters( { { { 0.0 } } } );
                    tIQI->set_property( tPropAgg, "Exponent" );
                    tIQI->set_property( tPropShift, "Shift" );
                }

                for ( uint iGP = 0; iGP < tIntegPoints.n_cols(); iGP++ )
                {
                    tIQI->reset_eval_flags();
                    tIQI->mSet->mLeaderFIManager->set_space_time( tIntegPoints.get_column( iGP ) );

                    // (1) value vs closed form
                    real tVM       = reference_von_mises( tCM->flux() );
                    real tExpected = tCase.mKS ? std::exp( tCase.mParameter * ( tVM / tRef - 1.0 ) ) : std::pow( tVM / tRef, tCase.mParameter );
                    Matrix< DDRMat > tQI;
                    tIQI->compute_QI( tQI );
                    CHECK( std::abs( tQI( 0 ) - tExpected ) <= tValueTol * std::max( 1.0, std::abs( tExpected ) ) );

                    tIQI->mSet->mQI( 0 ).fill( 0.0 );
                    tIQI->compute_QI( 1.0 );
                    CHECK( std::abs( tIQI->mSet->mQI( 0 )( 0 ) - tExpected ) <= tValueTol * std::max( 1.0, std::abs( tExpected ) ) );

                    // (2) analytic dQI/du vs finite differences
                    tIQI->mSet->mResidual( 0 ).fill( 0.0 );
                    Matrix< DDRMat > tdQIdu;
                    Matrix< DDRMat > tdQIduFD;
                    bool             tCheck = tIQI->check_dQIdu_FD( 1.0, tPerturbation, tEpsilon, tdQIdu, tdQIduFD, true );
                    if ( !tCheck )
                    {
                        std::cout << "IQI_Max_Stress dQIdu mismatch: " << tCase.mName << ", " << tModelNames( iModel )
                                  << ", order " << iOrder << ", GP " << iGP << std::endl;
                    }
                    REQUIRE( tCheck );
                }

                tLeaderFIs.clear();
            }
        }
    }
} /*END_TEST_CASE*/
