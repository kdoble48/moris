/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Eigen_Value.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Eigen_Value.hpp"
#include "cl_FEM_Model.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Eigen_Value::IQI_Eigen_Value()
    {
        // set FEM IQI type
        mFEMIQIType = fem::IQI_Type::EIGEN_VECTOR;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Eigen_Value::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        IQI::set_parameters( aParameters );

        // size of parameter list
        uint tParamSize = mParameters.size();

        // check that eigen vector index is defined
        MORIS_ERROR( tParamSize == 1,
                "IQI_Eigen_Value::set_parameters - index of eigen vector needs to be specified.\n" );

        // extract spatial derivative order
        MORIS_ERROR( mParameters( 0 ).numel() == 1,
                "IQI_Eigen_Value::set_parameters - size of parameter list needs to be one.\n" );

        mEigenVectorIndex = mParameters( 0 )( 0 );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Eigen_Value::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        Matrix< DDRMat > tMat;

        this->evaluate_QI( tMat );

        mSet->get_QI()( tQIIndex ) += aWStar * tMat;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Eigen_Value::compute_QI( Matrix< DDRMat >& aQI )
    {
        // evaluate QI
        this->evaluate_QI( aQI );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Eigen_Value::evaluate_QI( Matrix< DDRMat >& aMat )
    {
        // get the eigen values
        std::shared_ptr< Vector< real > >& tEigenVals =
                mSet->get_fem_model()->get_eigen_values();

        // check that eigen values exist
        MORIS_ERROR( tEigenVals != nullptr && tEigenVals->size() != 0,
                "IQI_Eigen_Value::evaluate_QI - eigen values do not exist." );

        // the eigenvalue is selected by vectorial_field_index, independent of dof_quantity;
        // if that is unset, by the eigen vector index (function_parameters)
        uint tIndex = this->get_component_index( tEigenVals->size(), mEigenVectorIndex );

        aMat = { ( *tEigenVals )( tIndex ) };
    }

    //------------------------------------------------------------------------------

    void
    IQI_Eigen_Value::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false,
                "IQI_Eigen_Value::compute_dQIdu - not implemented." );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Eigen_Value::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false,
                "IQI_Eigen_Value::compute_dQIdu - not implemented." );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
