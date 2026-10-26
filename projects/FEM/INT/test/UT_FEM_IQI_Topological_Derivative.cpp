/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_IQI_Topological_Derivative.cpp
 *
 * IQI_Topological_Derivative (compliance topological derivative, Amstutz and Andrae 2006 eq. 9 in the local-stress
 * form) on one QUAD4 (unit square) with an affine displacement, so the strain is uniform and set exactly:
 *   U1 hole limit (r = 0, l = 0), plane stress, chi from the CM, against the independent Python reference
 *      td2d_common.td_compliance (studies/lunch_and_learn_xiga/scripts/td2d_common.py)
 *   U2 chi = 2 override: -3(r-1)/(2r+1) sigma:eps + delta l (Teixeira et al. 2025, eq. 15), r = 1e-4, both delta
 *   U3 chi from the model type: PLANE_STRAIN nu = 1/4 gives chi = 2 (U2 without override); nu = 1/3 gives
 *      chi = 5/3, whose hole-limit value is 8/9 of eq. 15 in pure shear and 4/3 of it under equibiaxial strain
 *   U4 r = 1 (no contrast) gives dJ = delta l
 * The IQI returns g = s delta dJ (s = output sign), so the checks below compare against s delta dJ.
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
    //------------------------------------------------------------------------------

    void
    tConstValFunction_UTIQITD(
            Matrix< DDRMat >&           aPropMatrix,
            Vector< Matrix< DDRMat > >& aParameters,
            Field_Interpolator_Manager* /*aFIManager*/ )
    {
        aPropMatrix = aParameters( 0 );
    }

    //------------------------------------------------------------------------------

    std::shared_ptr< Property >
    make_const_property_td( real aValue )
    {
        std::shared_ptr< Property > tProp = std::make_shared< Property >();
        tProp->set_parameters( { { { aValue } } } );
        tProp->set_val_function( tConstValFunction_UTIQITD );
        return tProp;
    }

    //------------------------------------------------------------------------------

    /**
     * One QUAD4 unit-square element with a linear elastic CM and a uniform strain ( exx, eyy, exy tensor shear ),
     * the IQI under test built through the factory.
     */
    class TD_Fixture
    {
      public:
        std::shared_ptr< Constitutive_Model > mCM;
        std::shared_ptr< IQI >                mIQI;
        Geometry_Interpolator*                mGI        = nullptr;
        Field_Interpolator_Manager*           mFIManager = nullptr;
        MSI::Equation_Set*                    mSet       = nullptr;

        TD_Fixture(
                real       aE,
                real       aNu,
                Model_Type aModelType,
                real       aExx,
                real       aEyy,
                real       aExy )
        {
            Vector< Vector< MSI::Dof_Type > > tDispDofTypes = { { MSI::Dof_Type::UX, MSI::Dof_Type::UY } };

            // constitutive model
            CM_Factory tCMFactory;
            mCM = tCMFactory.create_CM( Constitutive_Type::STRUC_LIN_ISO );
            mCM->set_dof_type_list( tDispDofTypes );
            mCM->set_property( make_const_property_td( aE ), "YoungsModulus" );
            mCM->set_property( make_const_property_td( aNu ), "PoissonRatio" );
            mCM->set_local_properties();

            // fem set
            mSet = new fem::Set();
            static_cast< fem::Set* >( mSet )->set_set_type( Element_Type::BULK );
            mCM->set_set_pointer( static_cast< fem::Set* >( mSet ) );

            mCM->mSet->mUniqueDofTypeList.resize( 100, MSI::Dof_Type::END_ENUM );
            mCM->mSet->mUniqueDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            mCM->mSet->mUniqueDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
            mCM->mSet->mLeaderDofTypeMap.set_size( static_cast< int >( MSI::Dof_Type::END_ENUM ) + 1, 1, -1 );
            mCM->mSet->mLeaderDofTypeMap( static_cast< int >( MSI::Dof_Type::UX ) ) = 0;
            mCM->get_global_dof_type_list();

            // model type before the space dimension (as in the FEM model initializer): set_space_dim sizes the
            // Voigt vectors for the model type
            mCM->set_model_type( aModelType );
            mCM->set_space_dim( 2 );

            // unit square geometry
            Matrix< DDRMat > tXHat = { { 0.0, 0.0 }, { 1.0, 0.0 }, { 1.0, 1.0 }, { 0.0, 1.0 } };
            Matrix< DDRMat > tTHat = { { 0.0 }, { 1.0 } };

            mtk::Interpolation_Rule tRule(
                    mtk::Geometry_Type::QUAD,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR,
                    mtk::Interpolation_Type::LAGRANGE,
                    mtk::Interpolation_Order::LINEAR );

            mGI = new Geometry_Interpolator( tRule );
            mGI->set_coeff( tXHat, tTHat );

            // affine displacement u = grad u * x with grad u = [ exx exy ; exy eyy ]: uniform strain,
            // engineering shear 2 exy; same coefficients at both time nodes
            Matrix< DDRMat > tUHat( 8, 2, 0.0 );
            for ( uint iTime = 0; iTime < 2; iTime++ )
            {
                for ( uint iNode = 0; iNode < 4; iNode++ )
                {
                    real tX = tXHat( iNode, 0 );
                    real tY = tXHat( iNode, 1 );

                    tUHat( 4 * iTime + iNode, 0 ) = aExx * tX + aExy * tY;
                    tUHat( 4 * iTime + iNode, 1 ) = aExy * tX + aEyy * tY;
                }
            }

            Vector< Field_Interpolator* > tFIs( 1 );
            tFIs( 0 ) = new Field_Interpolator( 2, tRule, mGI, tDispDofTypes( 0 ) );
            tFIs( 0 )->set_coeff( tUHat );

            Vector< Vector< gen::PDV_Type > >   tDummyDv;
            Vector< Vector< mtk::Field_Type > > tDummyField;
            mFIManager = new Field_Interpolator_Manager( tDispDofTypes, tDummyDv, tDummyField, mSet );

            mFIManager->mFI                     = tFIs;
            mFIManager->mIPGeometryInterpolator = mGI;
            mFIManager->mIGGeometryInterpolator = mGI;

            mCM->mSet->mLeaderFIManager = mFIManager;
            mCM->set_field_interpolator_manager( mFIManager );

            // IQI under test, through the factory
            IQI_Factory tIQIFactory;
            mIQI = tIQIFactory.create_IQI( IQI_Type::TOPOLOGICAL_DERIVATIVE );
            mIQI->set_constitutive_model( mCM, "Elast", mtk::Leader_Follower::LEADER );
            mIQI->set_name( "TD" );
            mIQI->set_set_pointer( static_cast< fem::Set* >( mSet ) );
            mIQI->set_field_interpolator_manager( mFIManager );

            // evaluation point (the strain is uniform, so any point gives the same state)
            Matrix< DDRMat > tParamPoint = { { 0.1 }, { -0.3 }, { 0.0 } };
            mCM->reset_eval_flags();
            mFIManager->set_space_time( tParamPoint );
        }

        ~TD_Fixture()
        {
            delete mFIManager;
            delete mGI;
        }

        // IQI value for the given parameters ( r / delta / l / chi / s )
        real
        evaluate( const Vector< Matrix< DDRMat > >& aParameters )
        {
            mIQI->set_parameters( aParameters );
            mCM->reset_eval_flags();

            Matrix< DDRMat > tQI;
            mIQI->compute_QI( tQI );
            REQUIRE( tQI.numel() == 1 );
            return tQI( 0 );
        }

        // in-plane sigma:eps from the CM (independent of the IQI)
        real
        flux_strain()
        {
            const Matrix< DDRMat >& tFlux   = mCM->flux();
            const Matrix< DDRMat >& tStrain = mCM->strain();
            uint                    tShear  = tFlux.numel() - 1;
            return tFlux( 0 ) * tStrain( 0 ) + tFlux( 1 ) * tStrain( 1 ) + tFlux( tShear ) * tStrain( tShear );
        }
    };

    //------------------------------------------------------------------------------

    bool
    close_relative( real aValue, real aReference, real aTolerance )
    {
        return std::abs( aValue - aReference ) <= aTolerance * std::max( 1.0, std::abs( aReference ) );
    }

    //------------------------------------------------------------------------------
}    // namespace

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Topological_Derivative_Factory", "[moris],[fem],[IQI_Topological_Derivative]" )
{
    // enum string round trip (the string list is generated by ENUM_MACRO)
    CHECK( IQI_Type_String::values( static_cast< uint >( IQI_Type::TOPOLOGICAL_DERIVATIVE ) ) == "TOPOLOGICAL_DERIVATIVE" );

    // the factory builds the class and sets its type
    IQI_Factory            tIQIFactory;
    std::shared_ptr< IQI > tIQI = tIQIFactory.create_IQI( IQI_Type::TOPOLOGICAL_DERIVATIVE );
    REQUIRE( tIQI != nullptr );
    CHECK( tIQI->get_fem_IQI_type() == IQI_Type::TOPOLOGICAL_DERIVATIVE );
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Topological_Derivative_U1_Python_Reference", "[moris],[fem],[IQI_Topological_Derivative]" )
{
    // U1: plane stress, hole limit r = 0, delta = -1, l = 0, chi from the CM (chi = (3 - nu)/(1 + nu)).
    // Reference: D_T J = td2d_common.td_compliance( exx, eyy, exy, E, nu ) (void-limit compliance TD, plane stress,
    // exy = tensor shear), file studies/lunch_and_learn_xiga/scripts/td2d_common.py:29, generated 2026-10-25 with
    //   rng = np.random.default_rng( 20261025 ); e = rng.normal( size = ( 3, 4 ) )
    //   td, _ = td_compliance( *e, E = 2.7, nu = 0.21 )
    // rows: exx, eyy, exy, D_T J (printed with %.17g)
    const real tE  = 2.7;
    const real tNu = 0.21;

    const Matrix< DDRMat > tReference = {
        { -0.47730393156071066, 1.6561217711207181, 0.54555488872667857, 27.189935310539976 },
        { 0.18852727277652609, 0.4253288505496946, 0.13503949858445136, 2.1060593257784253 },
        { 0.63606654904189575, 0.76086102380765563, -1.4119224949690927, 37.910330379599536 },
        { 0.73795631638967385, 0.70066761926686227, 0.095411746335976697, 9.0931647661518955 }
    };

    for ( uint iState = 0; iState < tReference.n_rows(); iState++ )
    {
        TD_Fixture tFixture( tE, tNu, Model_Type::PLANE_STRESS, tReference( iState, 0 ), tReference( iState, 1 ), tReference( iState, 2 ) );

        const real tDTJ = tReference( iState, 3 );

        // dJ = D_T J for the hole; default output sign s = +1 gives g = delta dJ = -D_T J
        CHECK( close_relative( tFixture.evaluate( { { { 0.0 } }, { { -1.0 } } } ), -tDTJ, 1.0e-10 ) );

        // s = -1 gives g = +D_T J
        CHECK( close_relative( tFixture.evaluate( { { { 0.0 } }, { { -1.0 } }, { { 0.0 } }, { { 0.0 } }, { { -1.0 } } } ), tDTJ, 1.0e-10 ) );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Topological_Derivative_U2_Chi_Override", "[moris],[fem],[IQI_Topological_Derivative]" )
{
    // U2: chi = 2 override reproduces eq. 15, dJ = -3(r-1)/(2r+1) sigma:eps + delta l, for both delta, on a
    // plane-strain CM with nu = 0.3 (whose own chi = 1.8 is ignored)
    TD_Fixture tFixture( 1.7, 0.3, Model_Type::PLANE_STRAIN, 0.31, -0.12, 0.27 );

    const real tR  = 1.0e-4;
    const real tL  = 0.7;
    const real tSE = tFixture.flux_strain();

    REQUIRE( tSE > 0.0 );

    for ( real tDelta : { -1.0, 1.0 } )
    {
        for ( real tSign : { 1.0, -1.0 } )
        {
            const real tDJ = -3.0 * ( tR - 1.0 ) / ( 2.0 * tR + 1.0 ) * tSE + tDelta * tL;

            const real tG = tFixture.evaluate( { { { tR } }, { { tDelta } }, { { tL } }, { { 2.0 } }, { { tSign } } } );

            CHECK( close_relative( tG, tSign * tDelta * tDJ, 1.0e-12 ) );
        }
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Topological_Derivative_U3_Chi_From_Model_Type", "[moris],[fem],[IQI_Topological_Derivative]" )
{
    SECTION( "plane strain nu = 1/4: chi = 2, equals the chi = 2 override" )
    {
        TD_Fixture tFixture( 1.3, 0.25, Model_Type::PLANE_STRAIN, 0.4, 0.15, -0.22 );

        for ( real tDelta : { -1.0, 1.0 } )
        {
            const real tFromCM   = tFixture.evaluate( { { { 1.0e-4 } }, { { tDelta } }, { { 0.5 } } } );
            const real tOverride = tFixture.evaluate( { { { 1.0e-4 } }, { { tDelta } }, { { 0.5 } }, { { 2.0 } } } );

            CHECK( close_relative( tFromCM, tOverride, 1.0e-12 ) );
        }
    }

    SECTION( "plane strain nu = 1/3: chi = 5/3, ratio to eq. 15 is 8/9 in pure shear" )
    {
        // pure shear: exx = -eyy, tr( eps ) = tr( sigma ) = 0
        TD_Fixture tFixture( 1.3, 1.0 / 3.0, Model_Type::PLANE_STRAIN, 0.2, -0.2, 0.0 );

        const real tFromCM = tFixture.evaluate( { { { 0.0 } }, { { -1.0 } } } );
        const real tEq15   = tFixture.evaluate( { { { 0.0 } }, { { -1.0 } }, { { 0.0 } }, { { 2.0 } } } );

        // independent value of eq. 15 at r = 0: 3 sigma:eps, output g = -dJ
        CHECK( close_relative( tEq15, -3.0 * tFixture.flux_strain(), 1.0e-12 ) );
        CHECK( close_relative( tFromCM / tEq15, 8.0 / 9.0, 1.0e-12 ) );
    }

    SECTION( "plane strain nu = 1/3: ratio to eq. 15 is 4/3 under equibiaxial strain" )
    {
        // equibiaxial in-plane stress: exx = eyy, no shear
        TD_Fixture tFixture( 1.3, 1.0 / 3.0, Model_Type::PLANE_STRAIN, 0.2, 0.2, 0.0 );

        const real tFromCM = tFixture.evaluate( { { { 0.0 } }, { { -1.0 } } } );
        const real tEq15   = tFixture.evaluate( { { { 0.0 } }, { { -1.0 } }, { { 0.0 } }, { { 2.0 } } } );

        CHECK( close_relative( tFromCM / tEq15, 4.0 / 3.0, 1.0e-12 ) );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "IQI_Topological_Derivative_U4_No_Contrast", "[moris],[fem],[IQI_Topological_Derivative]" )
{
    // U4: r = 1 means no change of material, so dJ = delta l whatever the strain and chi
    TD_Fixture tFixture( 2.0, 0.3, Model_Type::PLANE_STRESS, 0.5, -0.3, 0.4 );

    const real tL = 1.25;

    for ( real tDelta : { -1.0, 1.0 } )
    {
        // chi from the CM
        CHECK( close_relative( tFixture.evaluate( { { { 1.0 } }, { { tDelta } }, { { tL } } } ), tDelta * tDelta * tL, 1.0e-14 ) );

        // chi override
        CHECK( close_relative( tFixture.evaluate( { { { 1.0 } }, { { tDelta } }, { { tL } }, { { 2.0 } } } ), tDelta * tDelta * tL, 1.0e-14 ) );
    }
}

//------------------------------------------------------------------------------
