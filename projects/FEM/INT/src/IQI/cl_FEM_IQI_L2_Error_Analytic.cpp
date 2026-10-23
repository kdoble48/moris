/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_L2_Error_Analytic.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_L2_Error_Analytic.hpp"

#include "fn_norm.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IQI_L2_Error_Analytic::IQI_L2_Error_Analytic()
    {
        // set FEM IQI type
        mFEMIQIType = fem::IQI_Type::L2_ERROR_ANALYTIC;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "L2Check" ] = static_cast< uint >( IQI_Property_Type::L2_CHECK );
    }

    //------------------------------------------------------------------------------

    void
    IQI_L2_Error_Analytic::compute_QI( Matrix< DDRMat >& aQI )
    {
        // check that mQuantityDofType is set
        MORIS_ERROR( mQuantityDofType.size() > 0,
                "IQI_L2_Error_Analytic::compute_QI - mQuantityDofType not set." );

        // get field interpolator
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        real tJumpNorm;

        if ( mLeaderFieldTypes.size() != 0 )
        {
            // reference solution given as a field
            Field_Interpolator* tFIField =
                    mLeaderFIManager->get_field_interpolators_for_type( mLeaderFieldTypes( 0 )( 0 ) );

            tJumpNorm = norm( tFI->val() - tFIField->val() );
        }
        else
        {
            // reference solution given as a property
            std::shared_ptr< Property >& tPropL2Check =
                    mLeaderProp( static_cast< uint >( IQI_Property_Type::L2_CHECK ) );

            MORIS_ERROR( tPropL2Check != nullptr,
                    "IQI_L2_Error_Analytic::compute_QI - IQI %s: neither an L2Check property nor a field is set.",
                    mName.c_str() );

            tJumpNorm = norm( tFI->val() - tPropL2Check->val() );
        }

        // evaluate the QI
        aQI = { { tJumpNorm * tJumpNorm } };
    }

    //------------------------------------------------------------------------------

    void
    IQI_L2_Error_Analytic::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // evaluate the QI
        Matrix< DDRMat > tQI;
        this->compute_QI( tQI );

        mSet->get_QI()( tQIIndex ) += aWStar * tQI( 0 );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
