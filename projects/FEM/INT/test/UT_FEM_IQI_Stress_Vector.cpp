/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Stress_Vector.cpp
 *
 * STRESS_VECTOR as a global / criterion QI (MB-067): the QI is sized by get_matrix_dim and
 * compute_QI( real ) adds the stress vector into it. Both must use the standardized layout
 * ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy ) in every dimension and model type.
 * One element with u = A x, plane stress, plane strain, axisymmetric and 3D.
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
// MTK/src
#include "cl_MTK_Enums.hpp"

using namespace moris;
using namespace fem;

namespace
{
    // one linear elastic element ( QUAD4 on [1,2]x[0,1] or HEX8 unit cube ) with u = A x, and a set holding one QI
    struct Stress_Vector_Fixture
    {
        Stress_Vector_Fixture(
                Model_Type              aModelType,
                const Matrix< DDRMat >& aA )
        {
            uint tSpaceDim = aA.n_rows();

            mtk::Geometry_Type tGeometryType = tSpaceDim == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;
            mtk::Interpolation_Rule tRule( tGeometryType,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            Matrix< DDRMat > tXHat;
            if ( tSpaceDim == 2 )
            {
                tXHat = { { 1.0, 0.0 }, { 2.0, 0.0 }, { 2.0, 1.0 }, { 1.0, 1.0 } };
                mDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };
            }
            else
            {
                tXHat = { { 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0 }, { 1.0, 1.0, 0.0 }, { 0.0, 1.0, 0.0 },
                    { 0.0, 0.0, 1.0 }, { 1.0, 0.0, 1.0 }, { 1.0, 1.0, 1.0 }, { 0.0, 1.0, 1.0 } };
                mDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } };
            }
            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

            mGI = new Geometry_Interpolator( tRule );
            mGI->set_coeff( tXHat, tTHat );

            // nodal coefficients of u = A x at both time levels
            uint             tNumNodes = tXHat.n_rows();
            Matrix< DDRMat > tUHat( 2 * tNumNodes, tSpaceDim, 0.0 );
            for ( uint iNode = 0; iNode < 2 * tNumNodes; iNode++ )
            {
                for ( uint iI = 0; iI < tSpaceDim; iI++ )
                {
                    for ( uint iJ = 0; iJ < tSpaceDim; iJ++ )
                    {
                        tUHat( iNode, iI ) += aA( iI, iJ ) * tXHat( iNode % tNumNodes, iJ );
                    }
                }
            }

            Vector< Field_Interpolator* > tFIs( 1 );
            tFIs( 0 ) = new Field_Interpolator( tSpaceDim, tRule, mGI, mDofTypes( 0 ) );
            tFIs( 0 )->set_coeff( tUHat );

            // set with the displacement dofs and one QI
            mSet = new fem::Set();
            mSet->set_set_type( fem::Element_Type::BULK );
            mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );
            mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
            mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;

            moris::map< std::string, moris_index > tRequestedIQINamesAssemblyMap;
            tRequestedIQINamesAssemblyMap[ "IQI" ] = 0;
            mSet->mRequestedIQINamesAssemblyMap    = tRequestedIQINamesAssemblyMap;
            mSet->mQI.resize( 1 );

            Vector< Vector< gen::PDV_Type > >   tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            mFIManager                          = new Field_Interpolator_Manager( mDofTypes, tDummyDv, tDummyField, mSet );
            mFIManager->mFI                     = tFIs;
            mFIManager->mIPGeometryInterpolator = mGI;
            mFIManager->mIGGeometryInterpolator = mGI;
            mSet->mLeaderFIManager              = mFIManager;

            // properties: E = 1, nu = 0.3; axisymmetric: rotation axis data {{2 pi r},{r},{n1},{n2}} at r = 1.5
            mPropEMod = std::make_shared< fem::Property >();
            mPropEMod->set_parameters( { { { 1.0 } } } );
            mPropNu = std::make_shared< fem::Property >();
            mPropNu->set_parameters( { { { 0.3 } } } );
            mPropRotAxis = std::make_shared< fem::Property >();
            mPropRotAxis->set_parameters( { { { 3.0 * M_PI }, { 1.5 }, { 1.0 }, { 0.0 } } } );

            fem::CM_Factory tCMFactory;
            mCM = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            mCM->set_dof_type_list( { mDofTypes } );
            mCM->set_property( mPropEMod, "YoungsModulus" );
            mCM->set_property( mPropNu, "PoissonRatio" );
            if ( aModelType == Model_Type::AXISYMMETRIC )
            {
                mCM->set_property( mPropRotAxis, "AxisymRotationAxis" );
            }
            mCM->set_model_type( aModelType );
            mCM->set_local_properties();
            mCM->set_set_pointer( mSet );
            mCM->get_global_dof_type_list();
            mCM->set_space_dim( tSpaceDim );
            mCM->set_field_interpolator_manager( mFIManager );

            // evaluation point ( xi, eta, [zeta,] tau )
            Matrix< DDRMat > tParamPoint( tSpaceDim + 1, 1, 0.0 );
            tParamPoint( 0 ) = 0.1;
            tParamPoint( 1 ) = -0.2;
            mFIManager->set_space_time( tParamPoint );
            mCM->reset_eval_flags();
        }

        ~Stress_Vector_Fixture()
        {
            // the set owns the FI manager, which owns the FIs
            delete mSet;
            delete mGI;
        }

        std::shared_ptr< IQI >
        make_stress_vector_iqi()
        {
            fem::IQI_Factory       tFactory;
            std::shared_ptr< IQI > tIQI = tFactory.create_IQI( fem::IQI_Type::STRESS_VECTOR );
            tIQI->set_name( "IQI" );
            tIQI->set_constitutive_model( mCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
            tIQI->set_dof_type_list( mDofTypes, mtk::Leader_Follower::LEADER );
            tIQI->set_set_pointer( mSet );
            tIQI->set_field_interpolator_manager( mFIManager );
            return tIQI;
        }

        // the CM flux sorted by hand into ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy )
        Matrix< DDRMat >
        expected_stress_vector()
        {
            const Matrix< DDRMat >& tFlux = mCM->flux();

            Matrix< DDRMat > tExpected( 6, 1, 0.0 );
            switch ( tFlux.numel() )
            {
                case 3:
                    tExpected( 0 ) = tFlux( 0 );
                    tExpected( 1 ) = tFlux( 1 );
                    tExpected( 5 ) = tFlux( 2 );
                    break;
                case 4:
                    tExpected( 0 ) = tFlux( 0 );
                    tExpected( 1 ) = tFlux( 1 );
                    tExpected( 2 ) = tFlux( 2 );
                    tExpected( 5 ) = tFlux( 3 );
                    break;
                default:
                    tExpected = tFlux;
            }
            return tExpected;
        }

        // size the set QI from get_matrix_dim (as FEM_Model::initialize_IQIs and Set::initialize_mQI do),
        // assemble one point with weight aWStar through compute_QI( real ) and check it against the CM flux
        void
        check_global_stress_vector( uint aNumFluxComponents )
        {
            std::shared_ptr< IQI > tIQI = this->make_stress_vector_iqi();

            REQUIRE( mCM->flux().numel() == aNumFluxComponents );

            // the QI must hold the 6-component vector; checked first, so a smaller QI fails here and not in the assembly
            std::pair< uint, uint > tDim = tIQI->get_matrix_dim();
            REQUIRE( tDim.first == 6 );
            REQUIRE( tDim.second == 1 );

            mSet->mQI( 0 ).set_size( tDim.first, tDim.second, 0.0 );

            real tWStar = 0.5;
            tIQI->compute_QI( tWStar );

            Matrix< DDRMat > tExpected = this->expected_stress_vector();
            const Matrix< DDRMat >& tQI = mSet->mQI( 0 );
            REQUIRE( tQI.n_rows() == 6 );
            REQUIRE( tQI.n_cols() == 1 );
            for ( uint i = 0; i < 6; i++ )
            {
                CHECK( tQI( i ) == Approx( tWStar * tExpected( i ) ).margin( 1.0e-15 ) );
            }

            // the output path ( compute_QI( Matrix ) ) returns the same vector, unweighted
            Matrix< DDRMat > tQIPoint;
            tIQI->compute_QI( tQIPoint );
            REQUIRE( tQIPoint.numel() == 6 );
            for ( uint i = 0; i < 6; i++ )
            {
                CHECK( tQIPoint( i ) == Approx( tExpected( i ) ).margin( 1.0e-15 ) );
            }
        }

        Vector< Vector< MSI::Dof_Type > >     mDofTypes;
        Geometry_Interpolator*                mGI        = nullptr;
        fem::Set*                             mSet       = nullptr;
        Field_Interpolator_Manager*           mFIManager = nullptr;
        std::shared_ptr< fem::Property >      mPropEMod;
        std::shared_ptr< fem::Property >      mPropNu;
        std::shared_ptr< fem::Property >      mPropRotAxis;
        std::shared_ptr< Constitutive_Model > mCM;
    };

    // a displacement gradient with nonzero normal and shear strains
    const Matrix< DDRMat > tA2D = { { 1.0e-3, 2.0e-4 }, { 1.0e-4, -5.0e-4 } };

    const Matrix< DDRMat > tA3D = {
        { 1.0e-3, 2.0e-4, 0.0 },
        { 1.0e-4, -5.0e-4, 3.0e-4 },
        { 0.0, 2.0e-4, 4.0e-4 }
    };
}    // namespace

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Stress_Vector_Global", "[moris],[fem],[IQI_Stress_Vector],[MB-067]" )
{
    // strains of u = A x: e_xx = 1e-3, e_yy = -5e-4, gamma_xy = 3e-4; E = 1, nu = 0.3
    real tExx = 1.0e-3, tEyy = -5.0e-4, tGxy = 3.0e-4, tNu = 0.3;

    SECTION( "2D plane stress" )
    {
        Stress_Vector_Fixture tFix( Model_Type::PLANE_STRESS, tA2D );

        // the CM flux is the plane stress solution
        real tPre = 1.0 / ( 1.0 - tNu * tNu );
        CHECK( tFix.mCM->flux()( 0 ) == Approx( tPre * ( tExx + tNu * tEyy ) ) );
        CHECK( tFix.mCM->flux()( 1 ) == Approx( tPre * ( tEyy + tNu * tExx ) ) );
        CHECK( tFix.mCM->flux()( 2 ) == Approx( tGxy / ( 2.0 * ( 1.0 + tNu ) ) ) );

        tFix.check_global_stress_vector( 3 );
    }

    SECTION( "2D plane strain" )
    {
        Stress_Vector_Fixture tFix( Model_Type::PLANE_STRAIN, tA2D );

        // the CM flux is the plane strain solution, with s_zz = nu ( s_xx + s_yy )
        real tPre = 1.0 / ( 1.0 + tNu ) / ( 1.0 - 2.0 * tNu );
        CHECK( tFix.mCM->flux()( 0 ) == Approx( tPre * ( ( 1.0 - tNu ) * tExx + tNu * tEyy ) ) );
        CHECK( tFix.mCM->flux()( 1 ) == Approx( tPre * ( ( 1.0 - tNu ) * tEyy + tNu * tExx ) ) );
        CHECK( tFix.mCM->flux()( 2 ) == Approx( tPre * tNu * ( tExx + tEyy ) ) );
        CHECK( tFix.mCM->flux()( 3 ) == Approx( tGxy / ( 2.0 * ( 1.0 + tNu ) ) ) );

        tFix.check_global_stress_vector( 4 );
    }

    SECTION( "2D axisymmetric" )
    {
        Stress_Vector_Fixture tFix( Model_Type::AXISYMMETRIC, tA2D );

        // nonzero hoop stress in the out-of-plane slot
        REQUIRE( std::abs( tFix.mCM->flux()( 2 ) ) > 1.0e-6 );

        tFix.check_global_stress_vector( 4 );
    }

    SECTION( "3D (control)" )
    {
        Stress_Vector_Fixture tFix( Model_Type::FULL, tA3D );

        tFix.check_global_stress_vector( 6 );
    }
}
