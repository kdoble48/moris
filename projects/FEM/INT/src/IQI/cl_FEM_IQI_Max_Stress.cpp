/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Max_Stress.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Max_Stress.hpp"
#include "cl_FEM_CM_Struc_Linear.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IQI_Max_Stress::IQI_Max_Stress( enum Stress_Type aStressType )
    {
        // assign stress type to evaluate
        mStressType = aStressType;

        // set FEM IQI type
        mFEMIQIType = fem::IQI_Type::MAX_STRESS;

        // set the property pointer cell size
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "ReferenceValue" ] = static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE );
        mPropertyMap[ "Exponent" ]       = static_cast< uint >( IQI_Property_Type::EXPONENT );
        mPropertyMap[ "Shift" ]          = static_cast< uint >( IQI_Property_Type::SHIFT );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::compute_QI( Matrix< DDRMat >& aQI )
    {
        // check if dof index was set
        if ( ( mLeaderDofTypes( 0 ).size() > 1 ) &&    //
                ( mStressType != Stress_Type::VON_MISES_STRESS ) )
        {
            MORIS_ERROR( mIQITypeIndex != -1,
                    "IQI_Max_Stress::compute_QI - mIQITypeIndex not set, but is needed for principal, normal and shear stresses." );

            MORIS_ERROR( mIQITypeIndex >= 0 && mIQITypeIndex <= 2,
                    "IQI_Max_Stress::compute_QI - mIQITypeIndex out of bounds, must be 0, 1, or 2 for the three spatial dimensions." );
        }

        // initialize stress value
        real tStressValue = 0.0;

        // switch for different stress types
        switch ( mStressType )
        {
            case Stress_Type::VON_MISES_STRESS:
                tStressValue = this->eval_Von_Mises_stress();
                break;

            case Stress_Type::PRINCIPAL_STRESS:
                tStressValue = this->eval_principal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::NORMAL_STRESS:
                tStressValue = this->eval_normal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::SHEAR_STRESS:
                tStressValue = this->eval_shear_stress( mIQITypeIndex + 1 );
                break;

            default:
                MORIS_ERROR( false, "IQI_Max_Stress::compute_QI - Unknown Stress Type." );
        }

        // check if properties set
        MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                "IQI_Max_Stress - no reference value set" );

        MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                "IQI_Max_Stress - no exponent set" );

        // get property values
        real tRefValue = mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
        real tExponent = mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );

        real tShift = 1.0;

        if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr )
            tShift = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );

        // evaluate the QI
        aQI = { { std::pow( 1.0 / tRefValue * tStressValue - tShift, tExponent ) } };
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // check if dof index was set
        if ( ( mLeaderDofTypes( 0 ).size() > 1 ) &&    //
                ( mStressType != Stress_Type::VON_MISES_STRESS ) )
        {
            MORIS_ERROR( mIQITypeIndex != -1,
                    "IQI_Max_Stress::compute_QI - mIQITypeIndex not set, but is needed for principal, normal and shear stresses." );

            MORIS_ERROR( mIQITypeIndex >= 0 && mIQITypeIndex <= 2,
                    "IQI_Max_Stress::compute_QI - mIQITypeIndex out of bounds, must be 0, 1, or 2 for the three spatial dimensions." );
        }

        // initialize stress value
        real tStressValue = 0.0;

        // switch for different stress types
        switch ( mStressType )
        {
            case Stress_Type::VON_MISES_STRESS:
                tStressValue = this->eval_Von_Mises_stress();
                break;

            case Stress_Type::PRINCIPAL_STRESS:
                tStressValue = this->eval_principal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::NORMAL_STRESS:
                tStressValue = this->eval_normal_stress( mIQITypeIndex + 1 );
                break;

            case Stress_Type::SHEAR_STRESS:
                tStressValue = this->eval_shear_stress( mIQITypeIndex + 1 );
                break;

            default:
                MORIS_ERROR( false, "IQI_Max_Stress::compute_QI - Unknown Stress Type." );
        }

        // check if properties set
        MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                "IQI_Max_Stress - no reference value set" );

        MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                "IQI_Max_Stress - no exponent set" );

        // get property values
        real tRefValue = mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
        real tExponent = mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );
        real tShift    = 1.0;
        if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr )
            tShift = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( std::pow( 1.0 / tRefValue * tStressValue - tShift, tExponent ) );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Max_Stress::eval_Von_Mises_stress()
    {
        // get standardized stress vector
        Matrix< DDRMat > tStressVector = this->get_stress_vector();

        // compute contributions to von mises stress
        real tNormalStressContribution =
                std::pow( tStressVector( 0 ) - tStressVector( 1 ), 2.0 ) +    //
                std::pow( tStressVector( 1 ) - tStressVector( 2 ), 2.0 ) +    //
                std::pow( tStressVector( 2 ) - tStressVector( 0 ), 2.0 );
        real tShearStressContribution =
                std::pow( tStressVector( 3 ), 2.0 ) +    //
                std::pow( tStressVector( 4 ), 2.0 ) +    //
                std::pow( tStressVector( 5 ), 2.0 );

        // compute von mises stress value
        return std::sqrt( 0.5 * tNormalStressContribution + 3.0 * tShearStressContribution );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Max_Stress::eval_principal_stress( uint aPrincipalStressIndex )
    {
        // get stress vector
        Matrix< DDRMat > tStressVector = this->get_stress_vector();

        // Standardized stress vector s = ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy ) (see get_stress_vector).
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
        if ( mLeaderFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions() == 2 )
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
                            "IQI_Max_Stress::eval_principal_stress - Only 1st and 2nd principal stresses known for 2D." );
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
            real tQ = ( 1. / 9. ) * ( std::pow( tI1, 2.0 ) - 3.0 * tI2 );
            real tR = ( 1. / 54. ) * ( 2.0 * std::pow( tI1, 3.0 ) - 9.0 * tI1 * tI2 + 27.0 * tI3 );
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
                            "IQI_Max_Stress::eval_principal_stress - Only 1st, 2nd, and 3rd principal stresses known for 3D." );
            }
        }

        // return stress value
        return tStressValue;
    }

    //------------------------------------------------------------------------------

    real
    IQI_Max_Stress::eval_normal_stress( uint aStressIndex )
    {
        // get stress vector
        Matrix< DDRMat > tStressVector = this->get_stress_vector();

        // pick stress value
        return tStressVector( aStressIndex - 1 );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Max_Stress::eval_shear_stress( uint aStressIndex )
    {
        // get stress vector
        Matrix< DDRMat > tStressVector = this->get_stress_vector();

        // pick stress value
        return tStressVector( aStressIndex + 2 );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    IQI_Max_Stress::get_stress_vector()
    {
        // get stress vector from Constitutive model and pull apart stress vector into components
        return standardize_stress_vector(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) )->flux() );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat >
    IQI_Max_Stress::standardize_stress_vector( const Matrix< DDRMat >& aCMStress )
    {
        // create stress vector
        Matrix< DDRMat > tStressVector( 6, 1, 0.0 );

        // pull apart stress vector into components
        uint tNumStressComponents = aCMStress.length();
        switch ( tNumStressComponents )
        {
            // 2D plane stress
            case 3:
                tStressVector( 0 ) = aCMStress( 0 );
                tStressVector( 1 ) = aCMStress( 1 );
                tStressVector( 5 ) = aCMStress( 2 );
                break;

                // 2D plane strain
            case 4:
                tStressVector( 0 ) = aCMStress( 0 );
                tStressVector( 1 ) = aCMStress( 1 );
                tStressVector( 2 ) = aCMStress( 2 );
                tStressVector( 5 ) = aCMStress( 3 );
                break;

                // 3D
            case 6:
                tStressVector = aCMStress;
                break;

                // Unknown size - error
            default:
                MORIS_ERROR( false,
                        "IQI_Max_Von_Mises_Stress::get_stress_vector - CM stress vector of unknown size; 3, 4 or 6 components expected." );
        }

        // return filled stress vector
        return tStressVector;
    }

    //------------------------------------------------------------------------------

    bool
    IQI_Max_Stress::supports_analytical_dQIdp_geometry()
    {
        // only the von Mises stress is differentiated
        if ( mStressType != Stress_Type::VON_MISES_STRESS )
        {
            return false;
        }

        // the spatial derivatives are implemented for the linear structural CMs
        const CM_Struc_Linear* tCMStrucLinear = dynamic_cast< const CM_Struc_Linear* >(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) ).get() );

        return tCMStrucLinear != nullptr && tCMStrucLinear->has_spatial_strain_derivatives();
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::compute_dQIdx_geometry_integrand(
            real              aWStar,
            Matrix< DDRMat >& adQIdx,
            Matrix< DDRMat >& adQIdn )
    {
        // get the elasticity CM
        CM_Struc_Linear* tCMElasticity = dynamic_cast< CM_Struc_Linear* >(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) ).get() );

        MORIS_ASSERT( tCMElasticity != nullptr,
                "IQI_Max_Stress::compute_dQIdx_geometry_integrand - elasticity CM is not a CM_Struc_Linear" );

        // number of space dimensions
        const uint tSpaceDim = mLeaderFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions();

        adQIdx.set_size( tSpaceDim, 1, 0.0 );
        adQIdn.set_size( 0, 0 );

        // get property values (held fixed in space)
        real tRefValue = mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );
        real tExponent = mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );
        real tShift    = 1.0;
        if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr )
        {
            tShift = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );
        }

        // von Mises stress and standardized stress vector
        Matrix< DDRMat > tStressVector = this->get_stress_vector();
        real             tVonMises     = this->eval_Von_Mises_stress();

        // derivative of the von Mises stress wrt the standardized stress vector (zero at zero stress)
        if ( tVonMises < MORIS_REAL_EPS )
        {
            return;
        }

        Matrix< DDRMat > tdVMdS( 6, 1, 0.0 );
        tdVMdS( 0 ) = ( 2.0 * tStressVector( 0 ) - tStressVector( 1 ) - tStressVector( 2 ) ) / ( 2.0 * tVonMises );
        tdVMdS( 1 ) = ( 2.0 * tStressVector( 1 ) - tStressVector( 2 ) - tStressVector( 0 ) ) / ( 2.0 * tVonMises );
        tdVMdS( 2 ) = ( 2.0 * tStressVector( 2 ) - tStressVector( 0 ) - tStressVector( 1 ) ) / ( 2.0 * tVonMises );
        tdVMdS( 3 ) = 3.0 * tStressVector( 3 ) / tVonMises;
        tdVMdS( 4 ) = 3.0 * tStressVector( 4 ) / tVonMises;
        tdVMdS( 5 ) = 3.0 * tStressVector( 5 ) / tVonMises;

        // derivative of the integrand ( sigma_vm / ref - shift )^p wrt sigma_vm
        const real tdQIdVM = tExponent * std::pow( tVonMises / tRefValue - tShift, tExponent - 1.0 ) / tRefValue;

        // constitutive matrix
        const Matrix< DDRMat >& tConst = tCMElasticity->constitutive();

        Matrix< DDRMat > tdStraindx;
        for ( uint iDir = 0; iDir < tSpaceDim; iDir++ )
        {
            // derivative of the flux wrt x_k, standardized
            tCMElasticity->eval_dStraindx( iDir, tdStraindx );
            Matrix< DDRMat > tdStressdx = standardize_stress_vector( tConst * tdStraindx );

            adQIdx( iDir ) = aWStar * tdQIdVM * dot( tdVMdS, tdStressdx );
        }
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
