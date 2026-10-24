/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_FEM_Stress_Measures.cpp
 *
 */

#include <cmath>
#include <catch.hpp>

#include "fn_FEM_Stress_Measures.hpp"

using namespace moris;
using namespace fem;

namespace
{
    //------------------------------------------------------------------------------
    // central finite difference of the von Mises stress wrt each flux component
    Matrix< DDRMat >
    fd_dvon_mises_dflux( const Matrix< DDRMat >& aFlux, real aPerturbation )
    {
        Matrix< DDRMat > tFD( aFlux.numel(), 1, 0.0 );
        for ( uint i = 0; i < aFlux.numel(); i++ )
        {
            Matrix< DDRMat > tPlus  = aFlux;
            Matrix< DDRMat > tMinus = aFlux;
            tPlus( i ) += aPerturbation;
            tMinus( i ) -= aPerturbation;
            tFD( i ) = ( eval_von_mises_stress( tPlus ) - eval_von_mises_stress( tMinus ) ) / ( 2.0 * aPerturbation );
        }
        return tFD;
    }
}    // namespace

//------------------------------------------------------------------------------

TEST_CASE( "Stress_Measures_Von_Mises", "[moris],[fem],[Stress_Measures]" )
{
    // one general stress state per flux layout: plane stress, plane strain / axisymmetric, 3D
    Vector< Matrix< DDRMat > > tFluxes = {
        { { 3.0 }, { -1.5 }, { 0.7 } },
        { { 3.0 }, { -1.5 }, { 0.4 }, { 0.7 } },
        { { 3.0 }, { -1.5 }, { 0.4 }, { -0.9 }, { 1.2 }, { 0.7 } }
    };

    for ( const Matrix< DDRMat >& tFlux : tFluxes )
    {
        // value against the closed form on the standardized vector
        Matrix< DDRMat > tS = standardize_stress_vector( tFlux );
        real tExpected      = std::sqrt( 0.5 * ( std::pow( tS( 0 ) - tS( 1 ), 2 ) + std::pow( tS( 1 ) - tS( 2 ), 2 ) + std::pow( tS( 2 ) - tS( 0 ), 2 ) )
                                    + 3.0 * ( tS( 3 ) * tS( 3 ) + tS( 4 ) * tS( 4 ) + tS( 5 ) * tS( 5 ) ) );
        CHECK( eval_von_mises_stress( tFlux ) == Approx( tExpected ).epsilon( 1e-14 ) );
        CHECK( eval_stress_measure( Stress_Type::VON_MISES_STRESS, tFlux, 0, 3 ) == Approx( tExpected ).epsilon( 1e-14 ) );

        // derivative in the flux layout against central finite differences
        Matrix< DDRMat > tdVMdFlux = eval_dvon_mises_stress_dflux( tFlux );
        Matrix< DDRMat > tFD       = fd_dvon_mises_dflux( tFlux, 1e-6 );

        REQUIRE( tdVMdFlux.numel() == tFlux.numel() );
        for ( uint i = 0; i < tFlux.numel(); i++ )
        {
            CHECK( tdVMdFlux( i ) == Approx( tFD( i ) ).epsilon( 1e-7 ).margin( 1e-9 ) );
        }

        // the stress-type dispatch returns the same derivative
        Matrix< DDRMat > tdSdFlux = eval_dstress_measure_dflux( Stress_Type::VON_MISES_STRESS, tFlux );
        for ( uint i = 0; i < tFlux.numel(); i++ )
        {
            CHECK( tdSdFlux( i ) == tdVMdFlux( i ) );
        }
    }

    // stress-free point: the von Mises stress is zero and not differentiable; the derivative is defined as zero
    Matrix< DDRMat > tZero( 6, 1, 0.0 );
    CHECK( eval_von_mises_stress( tZero ) == 0.0 );
    Matrix< DDRMat > tdZero = eval_dvon_mises_stress_dflux( tZero );
    for ( uint i = 0; i < 6; i++ )
    {
        CHECK( tdZero( i ) == 0.0 );
    }
}

//------------------------------------------------------------------------------

TEST_CASE( "Stress_Measures_Components_And_Principal", "[moris],[fem],[Stress_Measures]" )
{
    // 3D flux ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy ) with a known principal frame: diag( 5, 2, -1 ) rotated about z by 30 deg
    real tC = std::cos( M_PI / 6.0 );
    real tSn = std::sin( M_PI / 6.0 );
    real tSxx = 5.0 * tC * tC + 2.0 * tSn * tSn;
    real tSyy = 5.0 * tSn * tSn + 2.0 * tC * tC;
    real tSxy = ( 5.0 - 2.0 ) * tC * tSn;
    Matrix< DDRMat > tFlux3D = { { tSxx }, { tSyy }, { -1.0 }, { 0.0 }, { 0.0 }, { tSxy } };

    CHECK( eval_principal_stress( tFlux3D, 1, 3 ) == Approx( 5.0 ).epsilon( 1e-12 ) );
    CHECK( eval_principal_stress( tFlux3D, 2, 3 ) == Approx( 2.0 ).epsilon( 1e-12 ) );
    CHECK( eval_principal_stress( tFlux3D, 3, 3 ) == Approx( -1.0 ).epsilon( 1e-12 ) );

    // plane stress: the in-plane principal stresses of the same rotated state
    Matrix< DDRMat > tFlux2D = { { tSxx }, { tSyy }, { tSxy } };
    CHECK( eval_principal_stress( tFlux2D, 1, 2 ) == Approx( 5.0 ).epsilon( 1e-12 ) );
    CHECK( eval_principal_stress( tFlux2D, 2, 2 ) == Approx( 2.0 ).epsilon( 1e-12 ) );

    // normal and shear components follow the standardized ordering
    CHECK( eval_normal_stress( tFlux3D, 1 ) == tSxx );
    CHECK( eval_normal_stress( tFlux3D, 3 ) == -1.0 );
    CHECK( eval_shear_stress( tFlux3D, 3 ) == tSxy );
    CHECK( eval_shear_stress( tFlux2D, 3 ) == tSxy );
    CHECK( eval_normal_stress( tFlux2D, 3 ) == 0.0 );

    // flux layout round trip
    for ( const Matrix< DDRMat >& tFlux : { tFlux2D, tFlux3D } )
    {
        Matrix< DDRMat > tBack = stress_vector_to_flux_layout( standardize_stress_vector( tFlux ), tFlux.numel() );
        REQUIRE( tBack.numel() == tFlux.numel() );
        for ( uint i = 0; i < tFlux.numel(); i++ )
        {
            CHECK( tBack( i ) == tFlux( i ) );
        }
    }
}
