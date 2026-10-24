/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Stress_Measures.cpp
 *
 */

#include "fn_FEM_Stress_Measures.hpp"

#include <cmath>

#include "moris_typedefs.hpp"
#include "fn_assert.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    standardize_stress_vector( const Matrix< DDRMat >& aFlux )
    {
        // create stress vector
        Matrix< DDRMat > tStressVector( 6, 1, 0.0 );

        // pull apart stress vector into components
        switch ( aFlux.numel() )
        {
            // 2D plane stress
            case 3:
                tStressVector( 0 ) = aFlux( 0 );
                tStressVector( 1 ) = aFlux( 1 );
                tStressVector( 5 ) = aFlux( 2 );
                break;

            // 2D plane strain and axisymmetric
            case 4:
                tStressVector( 0 ) = aFlux( 0 );
                tStressVector( 1 ) = aFlux( 1 );
                tStressVector( 2 ) = aFlux( 2 );
                tStressVector( 5 ) = aFlux( 3 );
                break;

            // 3D
            case 6:
                for ( uint i = 0; i < 6; i++ )
                {
                    tStressVector( i ) = aFlux( i );
                }
                break;

            // Unknown size - error
            default:
                MORIS_ERROR( false,
                        "standardize_stress_vector - CM stress vector of unknown size; 3, 4 or 6 components expected." );
        }

        // return filled stress vector
        return tStressVector;
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    stress_vector_to_flux_layout(
            const Matrix< DDRMat >& aStandardized,
            uint                    aNumFluxComponents )
    {
        switch ( aNumFluxComponents )
        {
            // 2D plane stress
            case 3:
                return { { aStandardized( 0 ) }, { aStandardized( 1 ) }, { aStandardized( 5 ) } };

            // 2D plane strain and axisymmetric
            case 4:
                return { { aStandardized( 0 ) }, { aStandardized( 1 ) }, { aStandardized( 2 ) }, { aStandardized( 5 ) } };

            // 3D
            case 6:
                return aStandardized;

            default:
                MORIS_ERROR( false,
                        "stress_vector_to_flux_layout - CM stress vector of unknown size; 3, 4 or 6 components expected." );
                return aStandardized;
        }
    }

    //------------------------------------------------------------------------------

    real
    eval_von_mises_stress( const Matrix< DDRMat >& aFlux )
    {
        // get standardized stress vector
        Matrix< DDRMat > tS = standardize_stress_vector( aFlux );

        // compute contributions to von mises stress
        real tNormalStressContribution =
                std::pow( tS( 0 ) - tS( 1 ), 2.0 ) +    //
                std::pow( tS( 1 ) - tS( 2 ), 2.0 ) +    //
                std::pow( tS( 2 ) - tS( 0 ), 2.0 );
        real tShearStressContribution =
                std::pow( tS( 3 ), 2.0 ) +    //
                std::pow( tS( 4 ), 2.0 ) +    //
                std::pow( tS( 5 ), 2.0 );

        // compute von mises stress value
        return std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    eval_dvon_mises_stress_dflux( const Matrix< DDRMat >& aFlux )
    {
        // standardized stress vector ( s11, s22, s33, s23, s13, s12 ) and the von Mises stress
        Matrix< DDRMat > tS        = standardize_stress_vector( aFlux );
        real             tVonMises = eval_von_mises_stress( aFlux );

        // d sigma_vm / d s_i in the standardized ordering
        Matrix< DDRMat > tdVMdS( 6, 1, 0.0 );
        if ( tVonMises > MORIS_REAL_EPS )
        {
            real tScale = 0.5 / tVonMises;
            tdVMdS( 0 ) = tScale * ( 2.0 * tS( 0 ) - tS( 1 ) - tS( 2 ) );
            tdVMdS( 1 ) = tScale * ( 2.0 * tS( 1 ) - tS( 2 ) - tS( 0 ) );
            tdVMdS( 2 ) = tScale * ( 2.0 * tS( 2 ) - tS( 0 ) - tS( 1 ) );
            tdVMdS( 3 ) = tScale * 6.0 * tS( 3 );
            tdVMdS( 4 ) = tScale * 6.0 * tS( 4 );
            tdVMdS( 5 ) = tScale * 6.0 * tS( 5 );
        }

        // map back to the CM flux components
        return stress_vector_to_flux_layout( tdVMdS, aFlux.numel() );
    }

    //------------------------------------------------------------------------------

    real
    eval_principal_stress(
            const Matrix< DDRMat >& aFlux,
            uint                    aPrincipalStressIndex,
            uint                    aNumSpaceDims )
    {
        // get stress vector
        Matrix< DDRMat > tStressVector = standardize_stress_vector( aFlux );

        // Standardized stress vector s = ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy ) (see standardize_stress_vector).
        // 2D: the in-plane principal stresses are
        //     sigma_1,2 = ( s_xx + s_yy ) / 2 +- sqrt( ( ( s_xx - s_yy ) / 2 )^2 + s_xy^2 )
        // (in-plane shear s_xy is entry 5 of the standardized vector, entry 2 is s_zz)

        // initialize stress value
        real tStressValue = ( tStressVector( 0 ) + tStressVector( 1 ) ) / 2.0;

        // get max shear stress
        real tMaxShearStress = std::sqrt( std::pow( ( tStressVector( 0 ) - tStressVector( 1 ) ) / 2.0, 2.0 ) +    //
                                          std::pow( tStressVector( 5 ), 2.0 ) );

        // switch between the 2D case and the 3D case
        // for 2D
        if ( aNumSpaceDims == 2 )
        {
            switch ( aPrincipalStressIndex )
            {
                case 1:
                    tStressValue += tMaxShearStress;
                    break;

                case 2:
                    tStressValue -= tMaxShearStress;
                    break;

                default:
                    MORIS_ERROR( false,
                            "eval_principal_stress - Only 1st and 2nd principal stresses known for 2D." );
            }
        }

        // for 3D
        else
        {
            // 3D: the principal stresses are the roots of det( sigma - lambda I ) = -lambda^3 + I1 lambda^2 - I2 lambda + I3 = 0,
            // with the stress invariants
            //     I1 = s_xx + s_yy + s_zz                                                          (trace)
            //     I2 = s_xx s_yy + s_yy s_zz + s_zz s_xx - s_yz^2 - s_xz^2 - s_xy^2               (sum of principal minors)
            //     I3 = s_xx s_yy s_zz + 2 s_yz s_xz s_xy - s_xx s_yz^2 - s_yy s_xz^2 - s_zz s_xy^2  (det sigma)
            // solved in closed form (trigonometric solution of the cubic):
            //     Q = ( I1^2 - 3 I2 ) / 9,   R = ( 2 I1^3 - 9 I1 I2 + 27 I3 ) / 54,   theta = acos( R / Q^(3/2) )
            //     lambda_k = 2 sqrt( Q ) cos( ( theta + 2 pi k ) / 3 ) + I1 / 3,   k = 0, 1, 2

            // compute invariants
            real tI1 = tStressVector( 0 ) + tStressVector( 1 ) + tStressVector( 2 );
            real tI2 = tStressVector( 0 ) * tStressVector( 1 ) +    //
                       tStressVector( 1 ) * tStressVector( 2 ) +    //
                       tStressVector( 2 ) * tStressVector( 0 ) -    //
                       std::pow( tStressVector( 3 ), 2.0 ) -        //
                       std::pow( tStressVector( 4 ), 2.0 ) -        //
                       std::pow( tStressVector( 5 ), 2.0 );
            real tI3 = tStressVector( 0 ) * tStressVector( 1 ) * tStressVector( 2 ) -    //
                       tStressVector( 0 ) * std::pow( tStressVector( 3 ), 2.0 ) -        //
                       tStressVector( 1 ) * std::pow( tStressVector( 4 ), 2.0 ) -        //
                       tStressVector( 2 ) * std::pow( tStressVector( 5 ), 2.0 ) +        //
                       2.0 * tStressVector( 3 ) * tStressVector( 4 ) * tStressVector( 5 );

            // help values
            real tQ = ( 1.0 / 9.0 ) * ( std::pow( tI1, 2.0 ) - 3.0 * tI2 );
            real tR = ( 1.0 / 54.0 ) * ( 2.0 * std::pow( tI1, 3.0 ) - 9.0 * tI1 * tI2 + 27.0 * tI3 );
            // Q = 0 is a hydrostatic state (all principal values I1/3); round-off can push Q below 0 or |R|/Q^1.5 above 1
            tQ      = std::max( tQ, 0.0 );
            real tT = tQ > 0.0 ? std::acos( std::max( -1.0, std::min( 1.0, tR / std::sqrt( std::pow( tQ, 3.0 ) ) ) ) ) : 0.0;

            // compute principal values
            real tE1 = 2.0 * std::sqrt( tQ ) * std::cos( ( tT ) / 3.0 ) + tI1 / 3.0;
            real tE2 = 2.0 * std::sqrt( tQ ) * std::cos( ( tT + 2.0 * M_PI ) / 3.0 ) + tI1 / 3.0;
            real tE3 = 2.0 * std::sqrt( tQ ) * std::cos( ( tT + 4.0 * M_PI ) / 3.0 ) + tI1 / 3.0;

            // figure out minimum and maximum values
            switch ( aPrincipalStressIndex )
            {
                case 1:
                    tStressValue = std::max( std::max( tE1, tE2 ), std::max( tE2, tE3 ) );
                    break;

                case 2:
                    tStressValue = std::max( std::min( tE1, tE2 ), std::min( std::max( tE1, tE2 ), tE3 ) );
                    break;

                case 3:
                    tStressValue = std::min( std::min( tE1, tE2 ), std::min( tE2, tE3 ) );
                    break;

                default:
                    MORIS_ERROR( false,
                            "eval_principal_stress - Only 1st, 2nd, and 3rd principal stresses known for 3D." );
            }
        }

        // return stress value
        return tStressValue;
    }

    //------------------------------------------------------------------------------

    real
    eval_normal_stress( const Matrix< DDRMat >& aFlux, uint aStressIndex )
    {
        // pick stress value
        return standardize_stress_vector( aFlux )( aStressIndex - 1 );
    }

    //------------------------------------------------------------------------------

    real
    eval_shear_stress( const Matrix< DDRMat >& aFlux, uint aStressIndex )
    {
        // pick stress value
        return standardize_stress_vector( aFlux )( aStressIndex + 2 );
    }

    //------------------------------------------------------------------------------

    real
    eval_stress_measure(
            Stress_Type             aStressType,
            const Matrix< DDRMat >& aFlux,
            uint                    aIndex,
            uint                    aNumSpaceDims )
    {
        // switch for different stress types
        switch ( aStressType )
        {
            case Stress_Type::VON_MISES_STRESS:
                return eval_von_mises_stress( aFlux );

            case Stress_Type::PRINCIPAL_STRESS:
                return eval_principal_stress( aFlux, aIndex, aNumSpaceDims );

            case Stress_Type::NORMAL_STRESS:
                return eval_normal_stress( aFlux, aIndex );

            case Stress_Type::SHEAR_STRESS:
                return eval_shear_stress( aFlux, aIndex );

            default:
                MORIS_ERROR( false, "eval_stress_measure - Unknown or non-scalar Stress Type." );
                return 0.0;
        }
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    eval_dstress_measure_dflux(
            Stress_Type             aStressType,
            const Matrix< DDRMat >& aFlux )
    {
        MORIS_ERROR( aStressType == Stress_Type::VON_MISES_STRESS,
                "eval_dstress_measure_dflux - the derivative is implemented for the von Mises stress only." );

        return eval_dvon_mises_stress_dflux( aFlux );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
