/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Stress.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Stress.hpp"
#include "cl_FEM_Model.hpp"
#include "fn_FEM_Stress_Measures.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Stress::IQI_Stress(
            enum Stress_Type      aStressType,
            enum CM_Function_Type aFluxType )
    {
        // assign stress type to evaluate
        mStressType = aStressType;

        // assign flux type for CM request
        mFluxType = aFluxType;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
    }

    //------------------------------------------------------------------------------

    void IQI_Stress::compute_QI( Matrix< DDRMat >& aQI )
    {
        // the full stress vector, or one of the scalar stress measures
        if ( mStressType == Stress_Type::STRESS_VECTOR )
        {
            this->get_stress_vector( aQI );
        }
        else
        {
            aQI = { { this->eval_stress_value() } };
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Stress::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // evaluate the stress measure
        Matrix< DDRMat > tStressVector;
        this->compute_QI( tStressVector );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( tStressVector );
    }

    //------------------------------------------------------------------------------

    real IQI_Stress::eval_stress_value()
    {
        // principal, normal and shear stresses need vectorial_field_index ( 0, 1 or 2 ), whatever the dof groups (MB-079)
        uint tComponentIndex = 0;
        if ( mStressType != Stress_Type::VON_MISES_STRESS )
        {
            tComponentIndex = this->get_component_index( 3, -1 );
        }

        // the in-plane / 3D distinction is needed by the principal stresses only
        uint tNumSpaceDims = 0;
        if ( mStressType == Stress_Type::PRINCIPAL_STRESS )
        {
            tNumSpaceDims = mLeaderFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions();
        }

        return eval_stress_measure(
                mStressType,
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) )->flux( mFluxType ),
                tComponentIndex + 1,
                tNumSpaceDims );
    }

    //------------------------------------------------------------------------------

    Matrix< DDRMat > IQI_Stress::eval_dStressdFlux()
    {
        return eval_dstress_measure_dflux(
                mStressType,
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) )->flux( mFluxType ) );
    }

    //------------------------------------------------------------------------------

    void IQI_Stress::get_stress_vector( Matrix< DDRMat >& aStressVector )
    {
        aStressVector = standardize_stress_vector(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) )->flux( mFluxType ) );
    }

    //------------------------------------------------------------------------------

    std::pair< uint, uint > IQI_Stress::get_matrix_dim()
    {
        // if stress type not specified as a vector return and scaler
        if ( mStressType != Stress_Type::STRESS_VECTOR )
        {
            return std::make_pair( 1, 1 );
        }

        // the stress vector is assembled in the standardized layout ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy )
        // in every dimension and model type (see get_stress_vector), not in the layout of the CM flux (MB-067)
        return std::make_pair( 6, 1 );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
