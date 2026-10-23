/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Formulas.cpp
 *
 * Exact values of IQIs against closed forms on one QUAD4 (unit square) with linear fields:
 * J_INTEGRAL (MB-053), POWER_DISSIPATION_BULK (MB-056), DRAG/LIFT_COEFF (MB-059),
 * CONTACT_PRESSURE with a linear CM (MB-061), HEAT_METHOD_PENALTY (MB-063).
 *
 */

#include <string>
#include <cmath>
#include <catch.hpp>
#include <memory>

#include "assert.hpp"

#define protected public
#define private public
// FEM/INT/src
#include "cl_FEM_Model.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
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

using namespace moris;
using namespace fem;

namespace
{
    void
    tConstValFunction_UTIQIFORMULAS(
            Matrix< DDRMat >&           aPropMatrix,
            Vector< Matrix< DDRMat > >& aParameters,
            Field_Interpolator_Manager* aFIManager )
    {
        aPropMatrix = aParameters( 0 );
    }

    std::shared_ptr< Property >
    make_const_property( const Matrix< DDRMat >& aValue )
    {
        std::shared_ptr< Property > tProp = std::make_shared< Property >();
        tProp->set_parameters( { aValue } );
        tProp->set_val_function( tConstValFunction_UTIQIFORMULAS );
        return tProp;
    }

    // displacement gradient G (u = G x), simple shear rate and pressure of the flow v = ( gamma y, 0 ), p = p0
    const real tG00 = 1.0e-3, tG01 = 2.0e-3, tG10 = 3.0e-3, tG11 = -1.0e-3;
    const real tGamma = 0.8, tP0 = 0.7;

    // one QUAD4 (unit square), linear in time; FIs: {UX,UY}, {VX,VY}, {P}, {TEMP}
    struct Formulas_Fixture
    {
        Vector< Vector< MSI::Dof_Type > > mDofTypes = {
            { MSI::Dof_Type::UX, MSI::Dof_Type::UY },
            { MSI::Dof_Type::VX, MSI::Dof_Type::VY },
            { MSI::Dof_Type::P },
            { MSI::Dof_Type::TEMP }
        };

        Geometry_Interpolator*      mGI        = nullptr;
        fem::Set*                   mSet       = nullptr;
        Field_Interpolator_Manager* mFIManager = nullptr;

        Formulas_Fixture( real aTempOffset, real aTempSlopeX )
        {
            mtk::Interpolation_Rule tRule( mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            mGI = new Geometry_Interpolator( tRule );

            Matrix< DDRMat > tXHat = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };
            mGI->set_coeff( tXHat, tTHat );

            // nodal values, repeated for both time levels
            Matrix< DDRMat > tUHat( 8, 2, 0.0 ), tVHat( 8, 2, 0.0 ), tPHat( 8, 1, tP0 ), tTempHat( 8, 1, 0.0 );
            for ( uint iNode = 0; iNode < 8; iNode++ )
            {
                real tX = tXHat( iNode % 4, 0 );
                real tY = tXHat( iNode % 4, 1 );

                tUHat( iNode, 0 )    = tG00 * tX + tG01 * tY;
                tUHat( iNode, 1 )    = tG10 * tX + tG11 * tY;
                tVHat( iNode, 0 )    = tGamma * tY;
                tTempHat( iNode, 0 ) = aTempOffset + aTempSlopeX * tX;
            }

            Vector< Field_Interpolator* > tFIs( 4 );
            tFIs( 0 ) = new Field_Interpolator( 2, tRule, mGI, mDofTypes( 0 ) );
            tFIs( 0 )->set_coeff( tUHat );
            tFIs( 1 ) = new Field_Interpolator( 2, tRule, mGI, mDofTypes( 1 ) );
            tFIs( 1 )->set_coeff( tVHat );
            tFIs( 2 ) = new Field_Interpolator( 1, tRule, mGI, mDofTypes( 2 ) );
            tFIs( 2 )->set_coeff( tPHat );
            tFIs( 3 ) = new Field_Interpolator( 1, tRule, mGI, mDofTypes( 3 ) );
            tFIs( 3 )->set_coeff( tTempHat );

            mSet = new fem::Set();
            mSet->set_set_type( fem::Element_Type::BULK );
            mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );
            mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            for ( uint iGroup = 0; iGroup < mDofTypes.size(); iGroup++ )
            {
                for ( MSI::Dof_Type tType : mDofTypes( iGroup ) )
                {
                    mSet->mUniqueDofTypeMap( static_cast< int >( tType ) ) = iGroup;
                    mSet->mLeaderDofTypeMap( static_cast< int >( tType ) ) = iGroup;
                }
            }

            Vector< Vector< gen::PDV_Type > >   tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            mFIManager                          = new Field_Interpolator_Manager( mDofTypes, tDummyDv, tDummyField, mSet );
            mFIManager->mFI                     = tFIs;
            mFIManager->mIPGeometryInterpolator = mGI;
            mFIManager->mIGGeometryInterpolator = mGI;
            mSet->mLeaderFIManager              = mFIManager;

            // evaluate at the element centre, x = ( 0.5, 0.5 )
            Matrix< DDRMat > tParamPoint = { { 0.0 }, { 0.0 }, { 0.0 } };
            mFIManager->set_space_time( tParamPoint );
        }

        ~Formulas_Fixture()
        {
            // the set owns the FI manager, which owns the FIs
            delete mSet;
            delete mGI;
        }

        // finish the wiring of an IQI and evaluate it
        real
        evaluate( const std::shared_ptr< IQI >& aIQI )
        {
            aIQI->set_dof_type_list( mDofTypes, mtk::Leader_Follower::LEADER );
            aIQI->set_set_pointer( mSet );
            aIQI->get_global_dof_type_list();
            aIQI->set_field_interpolator_manager( mFIManager );
            aIQI->reset_eval_flags();

            Matrix< DDRMat > tQI( 1, 1, 0.0 );
            aIQI->compute_QI( tQI );
            return tQI( 0 );
        }
    };

    std::shared_ptr< Constitutive_Model >
    make_plane_stress_cm()
    {
        CM_Factory                            tCMFactory;
        std::shared_ptr< Constitutive_Model > tCM = tCMFactory.create_CM( Constitutive_Type::STRUC_LIN_ISO );
        tCM->set_dof_type_list( { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } } );
        tCM->set_property( make_const_property( { { 1.0 } } ), "YoungsModulus" );
        tCM->set_property( make_const_property( { { 0.3 } } ), "PoissonRatio" );
        tCM->set_space_dim( 2 );
        tCM->set_model_type( Model_Type::PLANE_STRESS );
        tCM->set_local_properties();
        return tCM;
    }

    std::shared_ptr< Constitutive_Model >
    make_fluid_cm( real aViscosity, real aDensity )
    {
        CM_Factory                            tCMFactory;
        std::shared_ptr< Constitutive_Model > tCM = tCMFactory.create_CM( Constitutive_Type::FLUID_INCOMPRESSIBLE );
        tCM->set_dof_type_list( { { MSI::Dof_Type::VX, MSI::Dof_Type::VY }, { MSI::Dof_Type::P } } );
        tCM->set_property( make_const_property( { { aViscosity } } ), "Viscosity" );
        tCM->set_property( make_const_property( { { aDensity } } ), "Density" );
        tCM->set_space_dim( 2 );
        tCM->set_local_properties();
        return tCM;
    }

    // plane-stress stress of u = G x, E = 1, nu = 0.3: ( sxx, syy, sxy )
    void
    plane_stress( real& aSxx, real& aSyy, real& aSxy )
    {
        real tNu = 0.3;
        real tC  = 1.0 / ( 1.0 - tNu * tNu );
        aSxx     = tC * ( tG00 + tNu * tG11 );
        aSyy     = tC * ( tG11 + tNu * tG00 );
        aSxy     = tC * ( 1.0 - tNu ) / 2.0 * ( tG01 + tG10 );
    }
    // HEAT_METHOD_PENALTY on TEMP with B = g = gamma = 1, all weights 1, zero references
    std::shared_ptr< IQI >
    make_heat_penalty_iqi( sint aIndex )
    {
        IQI_Factory            tFactory;
        std::shared_ptr< IQI > tIQI = tFactory.create_IQI( IQI_Type::HEAT_METHOD_PENALTY );
        tIQI->set_name( "IQI" );
        tIQI->set_quantity_dof_type( { MSI::Dof_Type::TEMP } );
        tIQI->set_parameters( { { { 1.0 } }, { { 1.0 } }, { { 1.0 } }, { { 1.0 } }, { { 1.0 } }, { { 1.0 } }, { { 1.0 } } } );
        tIQI->set_property( make_const_property( { { 0.0 } } ), "L2_Reference" );
        tIQI->set_property( make_const_property( { { 0.0 }, { 0.0 } } ), "H1S_Reference" );
        tIQI->set_output_type_index( aIndex );
        return tIQI;
    }
}    // namespace

//------------------------------------------------------------------------------

TEST_CASE( "IQI_J_Integral_Exact", "[moris],[fem],[IQI_Formulas],[MB-053]" )
{
    Formulas_Fixture tFix( 0.0, 0.0 );

    real tSxx, tSyy, tSxy;
    plane_stress( tSxx, tSyy, tSxy );

    // strain energy density W = 1/2 sigma : epsilon (engineering shear)
    real tW = 0.5 * ( tSxx * tG00 + tSyy * tG11 + tSxy * ( tG01 + tG10 ) );

    for ( const Matrix< DDRMat >& tNormal : { Matrix< DDRMat >{ { 1.0 }, { 0.0 } }, Matrix< DDRMat >{ { 0.6 }, { 0.8 } } } )
    {
        IQI_Factory            tFactory;
        std::shared_ptr< IQI > tIQI = tFactory.create_IQI( IQI_Type::J_INTEGRAL );
        tIQI->set_name( "IQI" );
        tIQI->set_constitutive_model( make_plane_stress_cm(), "ElastLinIso" );
        tIQI->set_normal( tNormal );

        // J = W n_1 - t . du/dx_1, t = sigma n
        real tT0 = tSxx * tNormal( 0 ) + tSxy * tNormal( 1 );
        real tT1 = tSxy * tNormal( 0 ) + tSyy * tNormal( 1 );
        real tJ  = tW * tNormal( 0 ) - ( tT0 * tG00 + tT1 * tG10 );

        CHECK( tFix.evaluate( tIQI ) == Approx( tJ ).epsilon( 1e-10 ) );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Power_Dissipation_Bulk_Exact", "[moris],[fem],[IQI_Formulas],[MB-056]" )
{
    Formulas_Fixture tFix( 0.0, 0.0 );

    real tMu = 1.5, tRho = 1.0;

    IQI_Factory            tFactory;
    std::shared_ptr< IQI > tIQI = tFactory.create_IQI( IQI_Type::POWER_DISSIPATION_BULK );
    tIQI->set_name( "IQI" );
    tIQI->set_constitutive_model( make_fluid_cm( tMu, tRho ), "Fluid" );

    // simple shear v = ( gamma y, 0 ): 2 mu eps : eps = mu gamma^2; plus the rho |v|^2 term at y = 0.5
    real tExpected = tMu * tGamma * tGamma + tRho * std::pow( 0.5 * tGamma, 2 );

    CHECK( tFix.evaluate( tIQI ) == Approx( tExpected ).epsilon( 1e-10 ) );
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Drag_Lift_Sign", "[moris],[fem],[IQI_Formulas],[MB-059]" )
{
    Formulas_Fixture tFix( 0.0, 0.0 );

    real tMu = 1.5;

    // bottom wall: the fluid's outward normal is -e_y; sigma = [ -p0, mu gamma; mu gamma, -p0 ]
    // force on the body = -sigma n = ( mu gamma, -p0 )
    Matrix< DDRMat > tNormal = { { 0.0 }, { -1.0 } };

    for ( IQI_Type tType : { IQI_Type::DRAG_COEFF, IQI_Type::LIFT_COEFF } )
    {
        IQI_Factory            tFactory;
        std::shared_ptr< IQI > tIQI = tFactory.create_IQI( tType );
        tIQI->set_name( "IQI" );
        tIQI->set_constitutive_model( make_fluid_cm( tMu, 1.0 ), "Fluid" );
        tIQI->set_property( make_const_property( { { 1.0 } } ), "RefPressure" );
        tIQI->set_normal( tNormal );

        real tExpected = tType == IQI_Type::DRAG_COEFF ? tMu * tGamma : -tP0;

        CHECK( tFix.evaluate( tIQI ) == Approx( tExpected ).epsilon( 1e-10 ) );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Contact_Pressure_Linear", "[moris],[fem],[IQI_Formulas],[MB-061]" )
{
    Formulas_Fixture tFix( 0.0, 0.0 );

    real tSxx, tSyy, tSxy;
    plane_stress( tSxx, tSyy, tSxy );

    // small strain: both IQIs are N . sigma N with the facet normal
    Matrix< DDRMat > tNormal = { { 0.6 }, { 0.8 } };
    real             tExpected =
            tNormal( 0 ) * tNormal( 0 ) * tSxx + 2.0 * tNormal( 0 ) * tNormal( 1 ) * tSxy + tNormal( 1 ) * tNormal( 1 ) * tSyy;

    for ( IQI_Type tType : { IQI_Type::CONTACT_PRESSURE_REFERENCE, IQI_Type::CONTACT_PRESSURE_CURRENT } )
    {
        IQI_Factory            tFactory;
        std::shared_ptr< IQI > tIQI = tFactory.create_IQI( tType );
        tIQI->set_name( "IQI" );
        tIQI->set_constitutive_model( make_plane_stress_cm(), "TractionCM" );
        tIQI->set_normal( tNormal );

        CHECK( tFix.evaluate( tIQI ) == Approx( tExpected ).epsilon( 1e-10 ) );

        // no facet normal (requested on a bulk set): 0, not an abort
        tIQI->mNormal.set_size( 0, 0 );
        CHECK( tFix.evaluate( tIQI ) == 0.0 );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Heat_Method_Penalty_Finite", "[moris],[fem],[IQI_Formulas],[MB-063]" )
{
    SECTION( "far from the interface: finite, phi tilde = -B" )
    {
        // exponent -2 g PHID / B = 800: the exp form gave inf / inf = NaN
        Formulas_Fixture tFix( -400.0, 0.0 );

        real tQI = tFix.evaluate( make_heat_penalty_iqi( 0 ) );
        REQUIRE( std::isfinite( tQI ) );

        // phi tilde = -1, alpha = exp( -1 ), weight 1: L2 = 1, H1 = 0
        CHECK( tQI == Approx( 1.0 ) );
        CHECK( tFix.evaluate( make_heat_penalty_iqi( 3 ) ) == Approx( -1.0 ) );
    }

    SECTION( "near the interface: same values as the exp form" )
    {
        // PHID = -0.3 + 0.4 x, at x = 0.5: -0.1
        Formulas_Fixture tFix( -0.3, 0.4 );

        real tVal      = std::exp( -2.0 * ( -0.1 ) );
        real tPhiTilde = 2.0 / ( 1.0 + tVal ) - 1.0;
        real tGradNorm = 4.0 * tVal / std::pow( 1.0 + tVal, 2 ) * 0.4;

        CHECK( tFix.evaluate( make_heat_penalty_iqi( 3 ) ) == Approx( tPhiTilde ).epsilon( 1e-12 ) );
        CHECK( tFix.evaluate( make_heat_penalty_iqi( 4 ) ) == Approx( tGradNorm ).epsilon( 1e-12 ) );
    }
}
