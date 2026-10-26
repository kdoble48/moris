/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Stress_3D.cpp
 *
 */

#include <algorithm>
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
#include "cl_FEM_IQI_Zienkiewicz_Zhu.hpp"
// MTK/src
#include "cl_MTK_Enums.hpp"

using namespace moris;
using namespace fem;

namespace
{
    // HEX8 unit cube, standard Lagrange node order
    const Matrix< DDRMat > tXHatCube = {
        { 0.0, 0.0, 0.0 },
        { 1.0, 0.0, 0.0 },
        { 1.0, 1.0, 0.0 },
        { 0.0, 1.0, 0.0 },
        { 0.0, 0.0, 1.0 },
        { 1.0, 0.0, 1.0 },
        { 1.0, 1.0, 1.0 },
        { 0.0, 1.0, 1.0 }
    };

    // nodal coefficients of u = A x at both time levels (linear in time, 2 x 8 rows)
    Matrix< DDRMat >
    linear_displacement( const Matrix< DDRMat >& aA )
    {
        Matrix< DDRMat > tUHat( 16, 3, 0.0 );
        for ( uint iT = 0; iT < 2; iT++ )
        {
            for ( uint iNode = 0; iNode < 8; iNode++ )
            {
                for ( uint iI = 0; iI < 3; iI++ )
                {
                    for ( uint iJ = 0; iJ < 3; iJ++ )
                    {
                        tUHat( 8 * iT + iNode, iI ) += aA( iI, iJ ) * tXHatCube( iNode, iJ );
                    }
                }
            }
        }
        return tUHat;
    }

    // one side of an interface: a linear elastic CM on its own field interpolator manager
    struct Elastic_Side
    {
        Elastic_Side( MSI::Equation_Set* aSet, const Matrix< DDRMat >& aA )
                : mRule(
                          mtk::Geometry_Type::HEX,
                          mtk::Interpolation_Type::LAGRANGE,
                          mtk::Interpolation_Order::LINEAR,
                          mtk::Interpolation_Type::LAGRANGE,
                          mtk::Interpolation_Order::LINEAR )
                , mGI( mRule )
                , mFIManager( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } }, mDummyDv, mDummyField, aSet )
        {
            Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } };

            mPropEMod = std::make_shared< fem::Property >();
            mPropEMod->set_parameters( { { { 1.0 } } } );
            mPropNu = std::make_shared< fem::Property >();
            mPropNu->set_parameters( { { { 0.3 } } } );

            fem::CM_Factory tCMFactory;
            mCM = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
            mCM->set_dof_type_list( { tDispDofTypes } );
            mCM->set_property( mPropEMod, "YoungsModulus" );
            mCM->set_property( mPropNu, "PoissonRatio" );
            mCM->set_local_properties();
            mCM->set_set_pointer( static_cast< fem::Set* >( aSet ) );
            mCM->get_global_dof_type_list();
            mCM->set_space_dim( 3 );
            mCM->set_model_type( fem::Model_Type::FULL );

            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };
            mGI.set_coeff( tXHatCube, tTHat );

            // the field interpolator manager owns and deletes the field interpolator
            Vector< Field_Interpolator* > tFIs( 1 );
            tFIs( 0 ) = new Field_Interpolator( 3, mRule, &mGI, tDispDofTypes( 0 ) );
            tFIs( 0 )->set_coeff( linear_displacement( aA ) );

            mFIManager.mFI                     = tFIs;
            mFIManager.mIPGeometryInterpolator = &mGI;
            mFIManager.mIGGeometryInterpolator = &mGI;
            mCM->set_field_interpolator_manager( &mFIManager );

            // evaluation point, set on the interpolators directly (as in UT_FEM_IQI_Strain_Energy)
            Matrix< DDRMat > tParamPoint = { { 0.1 }, { -0.2 }, { 0.3 }, { 0.0 } };
            mGI.set_space_time( tParamPoint );
            tFIs( 0 )->set_space_time( tParamPoint );
            mCM->reset_eval_flags();
        }

        Vector< Vector< enum gen::PDV_Type > >   mDummyDv;
        Vector< Vector< enum mtk::Field_Type > > mDummyField;
        std::shared_ptr< fem::Property >         mPropEMod;
        std::shared_ptr< fem::Property >         mPropNu;
        std::shared_ptr< Constitutive_Model >    mCM;
        mtk::Interpolation_Rule                  mRule;
        Geometry_Interpolator                    mGI;
        Field_Interpolator_Manager               mFIManager;
    };

    MSI::Equation_Set*
    make_set()
    {
        MSI::Equation_Set* tSet = new fem::Set();
        auto*              tFemSet = static_cast< fem::Set* >( tSet );
        tFemSet->set_set_type( fem::Element_Type::BULK );
        tFemSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );
        tFemSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tFemSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
        tFemSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tFemSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
        // second DOF group, for the recovered stress field of the Zienkiewicz-Zhu test
        tFemSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::STRESS_DOF ) ) = 1;
        tFemSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::STRESS_DOF ) ) = 1;
        tFemSet->mFollowerDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
        tFemSet->mFollowerDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
        return tSet;
    }

    // a displacement gradient with every stress component nonzero and sigma_zz != sigma_yz
    const Matrix< DDRMat > tAGeneral = {
        { 1.0e-3, 2.0e-4, 0.0 },
        { 1.0e-4, -5.0e-4, 3.0e-4 },
        { 0.0, 2.0e-4, 4.0e-4 }
    };
}    // namespace

//------------------------------------------------------------------------------
// Principal stresses of a general 3D stress state (MB-046): their sum, pairwise-product sum and product must equal
// the three invariants of the CM stress, computed here from the full tensor, and they must be ordered.
//------------------------------------------------------------------------------

TEST_CASE( "IQI_Stress_Principal_3D", "[moris],[fem],[IQI_Stress_3D],[MB-046]" )
{
    MSI::Equation_Set* tSet = make_set();
    Elastic_Side       tSide( tSet, tAGeneral );

    const Matrix< DDRMat >& tS = tSide.mCM->flux();
    REQUIRE( tS.numel() == 6 );

    // tensor from Voigt (xx, yy, zz, yz, xz, xy)
    real tSxx = tS( 0 ), tSyy = tS( 1 ), tSzz = tS( 2 ), tSyz = tS( 3 ), tSxz = tS( 4 ), tSxy = tS( 5 );
    REQUIRE( std::abs( tSxy ) > 1.0e-6 );
    REQUIRE( std::abs( tSzz - tSyz ) > 1.0e-6 );

    real tI1 = tSxx + tSyy + tSzz;
    real tI2 = tSxx * tSyy + tSyy * tSzz + tSzz * tSxx - tSyz * tSyz - tSxz * tSxz - tSxy * tSxy;
    real tI3 = tSxx * ( tSyy * tSzz - tSyz * tSyz ) - tSxy * ( tSxy * tSzz - tSyz * tSxz ) + tSxz * ( tSxy * tSyz - tSyy * tSxz );

    fem::IQI_Factory tIQIFactory;
    Vector< real >   tPrincipal( 3 );
    for ( uint iP = 0; iP < 3; iP++ )
    {
        std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::PRINCIPAL_STRESS );
        tIQI->set_constitutive_model( tSide.mCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIQI->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } }, mtk::Leader_Follower::LEADER );
        tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
        tIQI->set_field_interpolator_manager( &tSide.mFIManager );
        tIQI->set_output_type_index( iP );

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        tPrincipal( iP ) = tQI( 0 );
    }

    real tScale = std::abs( tI1 ) + 1.0e-3;
    CHECK( tPrincipal( 0 ) >= tPrincipal( 1 ) );
    CHECK( tPrincipal( 1 ) >= tPrincipal( 2 ) );
    CHECK( std::abs( tPrincipal( 0 ) + tPrincipal( 1 ) + tPrincipal( 2 ) - tI1 ) < 1.0e-12 * tScale );
    CHECK( std::abs( tPrincipal( 0 ) * tPrincipal( 1 ) + tPrincipal( 1 ) * tPrincipal( 2 ) + tPrincipal( 2 ) * tPrincipal( 0 ) - tI2 )
            < 1.0e-12 * tScale * tScale );
    CHECK( std::abs( tPrincipal( 0 ) * tPrincipal( 1 ) * tPrincipal( 2 ) - tI3 ) < 1.0e-12 * tScale * tScale * tScale );

    // eigenvalue check: det( sigma - s I ) = 0 for each principal value s
    for ( uint iP = 0; iP < 3; iP++ )
    {
        real tA = tSxx - tPrincipal( iP ), tB = tSyy - tPrincipal( iP ), tC = tSzz - tPrincipal( iP );
        real tDet = tA * ( tB * tC - tSyz * tSyz ) - tSxy * ( tSxy * tC - tSyz * tSxz ) + tSxz * ( tSxy * tSyz - tB * tSxz );
        CHECK( std::abs( tDet ) < 1.0e-12 * tScale * tScale * tScale );
    }

    // the fem set is left for process teardown, as in the other FEM/INT unit tests
}

//------------------------------------------------------------------------------
// Hydrostatic state: all three principal stresses equal the mean stress, no NaN from acos( 0 / 0 )
//------------------------------------------------------------------------------

TEST_CASE( "IQI_Stress_Principal_3D_Hydrostatic", "[moris],[fem],[IQI_Stress_3D],[MB-046]" )
{
    MSI::Equation_Set* tSet = make_set();
    Matrix< DDRMat >   tA   = { { 1.0e-3, 0.0, 0.0 }, { 0.0, 1.0e-3, 0.0 }, { 0.0, 0.0, 1.0e-3 } };
    Elastic_Side       tSide( tSet, tA );

    real tMean = tSide.mCM->flux()( 0 );

    fem::IQI_Factory tIQIFactory;
    for ( uint iP = 0; iP < 3; iP++ )
    {
        std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::PRINCIPAL_STRESS );
        tIQI->set_constitutive_model( tSide.mCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIQI->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } }, mtk::Leader_Follower::LEADER );
        tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
        tIQI->set_field_interpolator_manager( &tSide.mFIManager );
        tIQI->set_output_type_index( iP );

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        CHECK( std::isfinite( tQI( 0 ) ) );
        CHECK( std::abs( tQI( 0 ) - tMean ) < 1.0e-12 * std::abs( tMean ) );
    }

    // the fem set is left for process teardown, as in the other FEM/INT unit tests
}

//------------------------------------------------------------------------------
// MAX_NORMAL_STRESS accepts a valid component index (MB-049): with reference 1, exponent 1 and shift 0 it returns the
// normal stress itself
//------------------------------------------------------------------------------

TEST_CASE( "IQI_Max_Normal_Stress_3D", "[moris],[fem],[IQI_Stress_3D],[MB-049]" )
{
    MSI::Equation_Set* tSet = make_set();
    Elastic_Side       tSide( tSet, tAGeneral );

    auto tPropOne  = std::make_shared< fem::Property >();
    auto tPropZero = std::make_shared< fem::Property >();
    tPropOne->set_parameters( { { { 1.0 } } } );
    tPropZero->set_parameters( { { { 0.0 } } } );

    fem::IQI_Factory tIQIFactory;
    for ( uint iC = 0; iC < 3; iC++ )
    {
        std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::MAX_NORMAL_STRESS );
        tIQI->set_constitutive_model( tSide.mCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIQI->set_property( tPropOne, "ReferenceValue", mtk::Leader_Follower::LEADER );
        tIQI->set_property( tPropOne, "Exponent", mtk::Leader_Follower::LEADER );
        tIQI->set_property( tPropZero, "Shift", mtk::Leader_Follower::LEADER );
        tIQI->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } }, mtk::Leader_Follower::LEADER );
        tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
        tIQI->set_field_interpolator_manager( &tSide.mFIManager );
        tIQI->set_output_type_index( iC );

        Matrix< DDRMat > tQI;
        REQUIRE_NOTHROW( tIQI->compute_QI( tQI ) );
        CHECK( std::abs( tQI( 0 ) - tSide.mCM->flux()( iC ) ) < 1.0e-14 );
    }

    // the fem set is left for process teardown, as in the other FEM/INT unit tests
}

//------------------------------------------------------------------------------
// Normal, shear and principal stresses need vectorial_field_index also when the first leader dof group has a single
// member (MB-079): NORMAL_STRESS read the standardized entry (uint) -1, SHEAR_STRESS returned s_zz
//------------------------------------------------------------------------------

TEST_CASE( "IQI_Stress_Index_Required_Single_Member_Dof_Group", "[moris],[fem],[IQI_Stress_3D],[MB-079]" )
{
    MSI::Equation_Set* tSet = make_set();
    Elastic_Side       tSide( tSet, tAGeneral );

    fem::IQI_Factory tIQIFactory;
    for ( fem::IQI_Type tType : { fem::IQI_Type::NORMAL_STRESS, fem::IQI_Type::SHEAR_STRESS, fem::IQI_Type::PRINCIPAL_STRESS } )
    {
        std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( tType );
        tIQI->set_constitutive_model( tSide.mCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIQI->set_dof_type_list( { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } },
                mtk::Leader_Follower::LEADER );
        tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
        tIQI->set_field_interpolator_manager( &tSide.mFIManager );

        // vectorial_field_index unset
        Matrix< DDRMat > tQI;
        CHECK_THROWS( tIQI->compute_QI( tQI ) );

        // out of range
        tIQI->set_output_type_index( 3 );
        CHECK_THROWS( tIQI->compute_QI( tQI ) );
    }

    // control: with the index set, NORMAL_STRESS 1 is s_yy and SHEAR_STRESS 2 is s_xy
    std::shared_ptr< fem::IQI > tNormal = tIQIFactory.create_IQI( fem::IQI_Type::NORMAL_STRESS );
    tNormal->set_constitutive_model( tSide.mCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tNormal->set_dof_type_list( { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } },
            mtk::Leader_Follower::LEADER );
    tNormal->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tNormal->set_field_interpolator_manager( &tSide.mFIManager );
    tNormal->set_output_type_index( 1 );

    Matrix< DDRMat > tQI;
    tNormal->compute_QI( tQI );
    CHECK( tQI( 0 ) == tSide.mCM->flux()( 1 ) );

    std::shared_ptr< fem::IQI > tShear = tIQIFactory.create_IQI( fem::IQI_Type::SHEAR_STRESS );
    tShear->set_constitutive_model( tSide.mCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tShear->set_dof_type_list( { { MSI::Dof_Type::TEMP }, { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ } },
            mtk::Leader_Follower::LEADER );
    tShear->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tShear->set_field_interpolator_manager( &tSide.mFIManager );
    tShear->set_output_type_index( 2 );

    tShear->compute_QI( tQI );
    CHECK( tQI( 0 ) == tSide.mCM->flux()( 5 ) );

    // the fem set is left for process teardown, as in the other FEM/INT unit tests
}

//------------------------------------------------------------------------------
// JUMP_TRACTION (MB-048): zero when both sides carry the same stress; (sigma_L - sigma_F) n when they differ
//------------------------------------------------------------------------------

TEST_CASE( "IQI_Jump_Traction_3D", "[moris],[fem],[IQI_Stress_3D],[MB-048]" )
{
    MSI::Equation_Set* tSet = make_set();

    Matrix< DDRMat > tNormal = { { 1.0 }, { 0.3 }, { 0.2 } };
    tNormal                  = tNormal / norm( tNormal );

    fem::IQI_Factory tIQIFactory;

    auto evaluate_jump = [ & ]( Elastic_Side& aLeader, Elastic_Side& aFollower, sint aIndex ) {
        std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::JUMP_TRACTION );
        tIQI->set_constitutive_model( aLeader.mCM, "TractionCM", mtk::Leader_Follower::LEADER );
        tIQI->set_constitutive_model( aFollower.mCM, "TractionCM", mtk::Leader_Follower::FOLLOWER );
        tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
        tIQI->set_normal( tNormal );
        tIQI->set_output_type_index( aIndex );
        aLeader.mCM->reset_eval_flags();
        aFollower.mCM->reset_eval_flags();

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        return tQI( 0 );
    };

    // equal stresses: the jump vanishes
    {
        Elastic_Side tLeader( tSet, tAGeneral );
        Elastic_Side tFollower( tSet, tAGeneral );
        real         tTraction = norm( tLeader.mCM->traction( tNormal ) );
        REQUIRE( tTraction > 1.0e-4 );
        CHECK( evaluate_jump( tLeader, tFollower, -1 ) < 1.0e-14 * tTraction );
    }

    // different stresses: each component equals ( sigma_L - sigma_F ) n
    {
        Matrix< DDRMat > tAFollower = 0.5 * tAGeneral;
        Elastic_Side     tLeader( tSet, tAGeneral );
        Elastic_Side     tFollower( tSet, tAFollower );
        Matrix< DDRMat > tExpected = tLeader.mCM->traction( tNormal ) - tFollower.mCM->traction( tNormal );
        for ( uint iC = 0; iC < 3; iC++ )
        {
            CHECK( std::abs( evaluate_jump( tLeader, tFollower, iC ) - tExpected( iC ) ) < 1.0e-14 );
        }
        CHECK( std::abs( evaluate_jump( tLeader, tFollower, -1 ) - norm( tExpected ) ) < 1.0e-14 );
    }

    // the fem set is left for process teardown, as in the other FEM/INT unit tests
}

//------------------------------------------------------------------------------
// 2D plane stress principal stresses (MB-051): sigma_1,2 = ( sxx + syy ) / 2 +- sqrt( ( ( sxx - syy ) / 2 )^2 + sxy^2 )
// with a nonzero in-plane shear
//------------------------------------------------------------------------------

TEST_CASE( "IQI_Stress_Principal_2D", "[moris],[fem],[IQI_Stress_3D],[MB-051]" )
{
    MSI::Equation_Set* tSet = make_set();

    Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };

    auto tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );
    auto tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    fem::CM_Factory                       tCMFactory;
    std::shared_ptr< Constitutive_Model > tCM = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCM->set_dof_type_list( { tDispDofTypes } );
    tCM->set_property( tPropEMod, "YoungsModulus" );
    tCM->set_property( tPropNu, "PoissonRatio" );
    tCM->set_local_properties();
    tCM->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tCM->get_global_dof_type_list();
    tCM->set_space_dim( 2 );
    tCM->set_model_type( fem::Model_Type::PLANE_STRESS );

    mtk::Interpolation_Rule tRule(
            mtk::Geometry_Type::QUAD,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );
    Geometry_Interpolator tGI( tRule );
    Matrix< DDRMat >      tXHat = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
    Matrix< DDRMat >      tTHat = { { 0.0 }, { 1.0 } };
    tGI.set_coeff( tXHat, tTHat );

    // u = A x with a dominant shear: A = [ 1e-3 4e-4; 2e-4 -5e-4 ]
    Matrix< DDRMat > tA    = { { 1.0e-3, 4.0e-4 }, { 2.0e-4, -5.0e-4 } };
    Matrix< DDRMat > tUHat( 8, 2, 0.0 );
    for ( uint iT = 0; iT < 2; iT++ )
    {
        for ( uint iNode = 0; iNode < 4; iNode++ )
        {
            for ( uint iI = 0; iI < 2; iI++ )
            {
                tUHat( 4 * iT + iNode, iI ) = tA( iI, 0 ) * tXHat( iNode, 0 ) + tA( iI, 1 ) * tXHat( iNode, 1 );
            }
        }
    }

    // the field interpolator manager owns and deletes the field interpolator
    Vector< Field_Interpolator* > tFIs( 1 );
    tFIs( 0 ) = new Field_Interpolator( 2, tRule, &tGI, tDispDofTypes( 0 ) );
    tFIs( 0 )->set_coeff( tUHat );

    Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
    Vector< Vector< enum mtk::Field_Type > > tDummyField;
    Field_Interpolator_Manager               tFIManager( tDispDofTypes, tDummyDv, tDummyField, tSet );
    tFIManager.mFI                     = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;
    tCM->set_field_interpolator_manager( &tFIManager );

    Matrix< DDRMat > tParamPoint = { { 0.2 }, { -0.1 }, { 0.0 } };
    tGI.set_space_time( tParamPoint );
    tFIs( 0 )->set_space_time( tParamPoint );
    tCM->reset_eval_flags();

    // plane stress flux: ( sxx, syy, sxy )
    const Matrix< DDRMat >& tS = tCM->flux();
    REQUIRE( tS.numel() == 3 );
    REQUIRE( std::abs( tS( 2 ) ) > 1.0e-5 );
    real tMean   = 0.5 * ( tS( 0 ) + tS( 1 ) );
    real tRadius = std::sqrt( std::pow( 0.5 * ( tS( 0 ) - tS( 1 ) ), 2.0 ) + std::pow( tS( 2 ), 2.0 ) );

    fem::IQI_Factory tIQIFactory;
    for ( uint iP = 0; iP < 2; iP++ )
    {
        std::shared_ptr< fem::IQI > tIQI = tIQIFactory.create_IQI( fem::IQI_Type::PRINCIPAL_STRESS );
        tIQI->set_constitutive_model( tCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
        tIQI->set_dof_type_list( tDispDofTypes, mtk::Leader_Follower::LEADER );
        tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
        tIQI->set_field_interpolator_manager( &tFIManager );
        tIQI->set_output_type_index( iP );

        Matrix< DDRMat > tQI;
        tIQI->compute_QI( tQI );
        real tExpected = iP == 0 ? tMean + tRadius : tMean - tRadius;
        CHECK( std::abs( tQI( 0 ) - tExpected ) < 1.0e-14 );
    }

    // the fem set is left for process teardown, as in the other FEM/INT unit tests
}

//------------------------------------------------------------------------------
// Zienkiewicz-Zhu indicator of the first principal stress (MB-046 in its own copy of eval_principal_stress). Decks can
// only request the von Mises ZZ IQI, so this branch is reachable from code only. With the recovered field s* set to
// the exact sigma_1 (an independent Jacobi eigenvalue computation), the indicator ( s* - sigma_1( u_h ) )^2 is 0.
//------------------------------------------------------------------------------

namespace
{
    // largest eigenvalue of a symmetric 3 x 3 matrix by Jacobi rotations
    real
    largest_eigenvalue( Matrix< DDRMat > aS )
    {
        for ( uint iSweep = 0; iSweep < 50; iSweep++ )
        {
            for ( uint p = 0; p < 2; p++ )
            {
                for ( uint q = p + 1; q < 3; q++ )
                {
                    if ( std::abs( aS( p, q ) ) < 1.0e-300 )
                    {
                        continue;
                    }
                    real tTheta = 0.5 * std::atan2( 2.0 * aS( p, q ), aS( q, q ) - aS( p, p ) );
                    real c = std::cos( tTheta ), s = std::sin( tTheta );
                    for ( uint k = 0; k < 3; k++ )
                    {
                        real tP = c * aS( p, k ) - s * aS( q, k ), tQ = s * aS( p, k ) + c * aS( q, k );
                        aS( p, k ) = tP, aS( q, k ) = tQ;
                    }
                    for ( uint k = 0; k < 3; k++ )
                    {
                        real tP = c * aS( k, p ) - s * aS( k, q ), tQ = s * aS( k, p ) + c * aS( k, q );
                        aS( k, p ) = tP, aS( k, q ) = tQ;
                    }
                }
            }
        }
        return std::max( { aS( 0, 0 ), aS( 1, 1 ), aS( 2, 2 ) } );
    }
}    // namespace

TEST_CASE( "IQI_Zienkiewicz_Zhu_Principal_3D", "[moris],[fem],[IQI_Stress_3D],[MB-046]" )
{
    MSI::Equation_Set* tSet = make_set();

    Vector< Vector< MSI::Dof_Type > > tDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY, MSI::Dof_Type::UZ }, { MSI::Dof_Type::STRESS_DOF } };

    auto tPropEMod = std::make_shared< fem::Property >();
    tPropEMod->set_parameters( { { { 1.0 } } } );
    auto tPropNu = std::make_shared< fem::Property >();
    tPropNu->set_parameters( { { { 0.3 } } } );

    fem::CM_Factory                       tCMFactory;
    std::shared_ptr< Constitutive_Model > tCM = tCMFactory.create_CM( fem::Constitutive_Type::STRUC_LIN_ISO );
    tCM->set_dof_type_list( { tDofTypes( 0 ) } );
    tCM->set_property( tPropEMod, "YoungsModulus" );
    tCM->set_property( tPropNu, "PoissonRatio" );
    tCM->set_local_properties();
    tCM->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tCM->get_global_dof_type_list();
    tCM->set_space_dim( 3 );
    tCM->set_model_type( fem::Model_Type::FULL );

    mtk::Interpolation_Rule tRule(
            mtk::Geometry_Type::HEX,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR,
            mtk::Interpolation_Type::LAGRANGE,
            mtk::Interpolation_Order::LINEAR );
    Geometry_Interpolator tGI( tRule );
    Matrix< DDRMat >      tTHat = { { 0.0 }, { 1.0 } };
    tGI.set_coeff( tXHatCube, tTHat );

    // the field interpolator manager owns and deletes the field interpolators
    Vector< Field_Interpolator* > tFIs( 2 );
    tFIs( 0 ) = new Field_Interpolator( 3, tRule, &tGI, tDofTypes( 0 ) );
    tFIs( 0 )->set_coeff( linear_displacement( tAGeneral ) );
    tFIs( 1 ) = new Field_Interpolator( 1, tRule, &tGI, tDofTypes( 1 ) );

    Vector< Vector< enum gen::PDV_Type > >   tDummyDv;
    Vector< Vector< enum mtk::Field_Type > > tDummyField;
    Field_Interpolator_Manager               tFIManager( tDofTypes, tDummyDv, tDummyField, tSet );
    tFIManager.mFI                     = tFIs;
    tFIManager.mIPGeometryInterpolator = &tGI;
    tFIManager.mIGGeometryInterpolator = &tGI;
    tCM->set_field_interpolator_manager( &tFIManager );

    Matrix< DDRMat > tParamPoint = { { 0.1 }, { -0.2 }, { 0.3 }, { 0.0 } };
    tGI.set_space_time( tParamPoint );
    tFIs( 0 )->set_space_time( tParamPoint );
    tCM->reset_eval_flags();

    // exact sigma_1 of the CM stress, independently of MORIS' invariant formula
    const Matrix< DDRMat >& tS     = tCM->flux();
    Matrix< DDRMat >        tSigma = { { tS( 0 ), tS( 5 ), tS( 4 ) }, { tS( 5 ), tS( 1 ), tS( 3 ) }, { tS( 4 ), tS( 3 ), tS( 2 ) } };
    real                    tSigma1 = largest_eigenvalue( tSigma );

    // recovered field s* = sigma_1 everywhere
    tFIs( 1 )->set_coeff( Matrix< DDRMat >( 16, 1, tSigma1 ) );
    tFIs( 1 )->set_space_time( tParamPoint );

    std::shared_ptr< fem::IQI > tIQI = std::make_shared< fem::IQI_Zienkiewicz_Zhu >( fem::Stress_Type::PRINCIPAL_STRESS );
    tIQI->set_constitutive_model( tCM, "ElastLinIso", mtk::Leader_Follower::LEADER );
    tIQI->set_dof_type_list( tDofTypes, mtk::Leader_Follower::LEADER );
    tIQI->set_quantity_dof_type( { MSI::Dof_Type::STRESS_DOF } );
    tIQI->set_set_pointer( static_cast< fem::Set* >( tSet ) );
    tIQI->set_field_interpolator_manager( &tFIManager );
    tIQI->set_output_type_index( 0 );

    Matrix< DDRMat > tQI;
    tIQI->compute_QI( tQI );

    // ( s* - sigma_1( u_h ) )^2, relative to sigma_1^2
    CHECK( tQI( 0 ) < 1.0e-24 * tSigma1 * tSigma1 );

    // the fem set is left for process teardown, as in the other FEM/INT unit tests
}
