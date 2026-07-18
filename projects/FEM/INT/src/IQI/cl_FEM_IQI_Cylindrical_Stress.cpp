/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Cylindrical_Stress.cpp
 *
 */

#include <cmath>

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Cylindrical_Stress.hpp"
#include "cl_FEM_Model.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Cylindrical_Stress::IQI_Cylindrical_Stress(
            enum CM_Function_Type aFluxType )
    {
        // assign IQI type
        mFEMIQIType = fem::IQI_Type::CYLINDRICAL_STRESS;

        // assign flux type for CM request
        mFluxType = aFluxType;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
    }

    //------------------------------------------------------------------------------

    void IQI_Cylindrical_Stress::compute_QI( Matrix< DDRMat >& aQI )
    {
        // evaluate the requested cylindrical stress component
        this->eval_cylindrical_stress( aQI );
    }

    //------------------------------------------------------------------------------

    void IQI_Cylindrical_Stress::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // evaluate the requested cylindrical stress component
        Matrix< DDRMat > tStressValue;
        this->eval_cylindrical_stress( tStressValue );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( tStressValue );
    }

    //------------------------------------------------------------------------------

    void IQI_Cylindrical_Stress::eval_cylindrical_stress( Matrix< DDRMat >& aStressValue )
    {
        // check that the component index was set
        MORIS_ERROR( mIQITypeIndex >= 0 && mIQITypeIndex <= 2,
                "IQI_Cylindrical_Stress::eval_cylindrical_stress - mIQITypeIndex out of bounds, "
                "must be 0 (sigma_rr), 1 (sigma_thetatheta), or 2 (sigma_rtheta)." );

        // check that the center point was provided through the function parameters
        MORIS_ERROR( mParameters.size() >= 1,
                "IQI_Cylindrical_Stress::eval_cylindrical_stress - mParameters not set; "
                "the center point ( cx, cy ) must be provided via function_parameters." );

        MORIS_ERROR( mParameters( 0 ).numel() == 2,
                "IQI_Cylindrical_Stress::eval_cylindrical_stress - the center point must contain "
                "exactly two entries ( cx, cy )." );

        // get standardized Cartesian stress vector
        Matrix< DDRMat > tStressVector;
        this->get_stress_vector( tStressVector );

        // in-plane Cartesian stress components
        real tSxx = tStressVector( 0 );
        real tSyy = tStressVector( 1 );
        real tSxy = tStressVector( 5 );

        // get the physical coordinates of the evaluation point
        const Matrix< DDRMat >& tX = mLeaderFIManager->get_IP_geometry_interpolator()->valx();

        // get the center point
        real tCx = mParameters( 0 )( 0 );
        real tCy = mParameters( 0 )( 1 );

        // polar angle of the evaluation point relative to the center
        real tTheta = std::atan2( tX( 1 ) - tCy, tX( 0 ) - tCx );
        real tC     = std::cos( tTheta );
        real tS     = std::sin( tTheta );

        // rotate the Cartesian stress into the polar frame and select the component
        real tValue = 0.0;
        switch ( mIQITypeIndex )
        {
            // radial stress sigma_rr
            case 0:
                tValue = tSxx * tC * tC + tSyy * tS * tS + 2.0 * tSxy * tS * tC;
                break;

            // hoop stress sigma_thetatheta
            case 1:
                tValue = tSxx * tS * tS + tSyy * tC * tC - 2.0 * tSxy * tS * tC;
                break;

            // shear stress sigma_rtheta
            case 2:
                tValue = ( tSyy - tSxx ) * tS * tC + tSxy * ( tC * tC - tS * tS );
                break;

            default:
                MORIS_ERROR( false,
                        "IQI_Cylindrical_Stress::eval_cylindrical_stress - Unknown component index." );
        }

        // return the scalar stress value
        aStressValue.set_size( 1, 1, tValue );
    }

    //------------------------------------------------------------------------------

    void IQI_Cylindrical_Stress::get_stress_vector( Matrix< DDRMat >& aStressVector )
    {
        // create stress vector
        aStressVector.set_size( 6, 1, 0.0 );

        // get stress vector from Constitutive model
        const Matrix< DDRMat >& tCMStress =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) )->flux( mFluxType );

        // pull apart stress vector into components
        uint tNumStressComponents = tCMStress.length();
        switch ( tNumStressComponents )
        {
            // 2D plane stress
            case 3:
            {
                aStressVector( 0 ) = tCMStress( 0 );
                aStressVector( 1 ) = tCMStress( 1 );
                aStressVector( 5 ) = tCMStress( 2 );
                break;
            }
            // 2D plane strain and axisymmetric
            case 4:
            {
                aStressVector( 0 ) = tCMStress( 0 );
                aStressVector( 1 ) = tCMStress( 1 );
                aStressVector( 2 ) = tCMStress( 2 );
                aStressVector( 5 ) = tCMStress( 3 );
                break;
            }
            // 3D
            case 6:
            {
                aStressVector = tCMStress;
                break;
            }
            // Unknown size - error
            default:
                MORIS_ERROR( false,
                        "IQI_Cylindrical_Stress::get_stress_vector - CM stress vector of unknown size; 3, 4 or 6 components expected." );
        }
    }

    //------------------------------------------------------------------------------

    std::pair< uint, uint > IQI_Cylindrical_Stress::get_matrix_dim()
    {
        // cylindrical stress component is a scalar
        return std::make_pair( 1, 1 );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
