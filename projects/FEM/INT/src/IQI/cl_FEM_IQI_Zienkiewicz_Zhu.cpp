/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Zienkiewicz_Zhu.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Zienkiewicz_Zhu.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Zienkiewicz_Zhu::IQI_Zienkiewicz_Zhu( enum Stress_Type aStressType )
            : IQI_Stress( aStressType )
    {
        MORIS_ERROR( aStressType != Stress_Type::STRESS_VECTOR,
                "IQI_Zienkiewicz_Zhu - the error indicator needs a scalar stress measure, not the stress vector." );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Zienkiewicz_Zhu::eval_error()
    {
        // get field interpolator of the recovered stress measure s*
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // e = s* - s( sigma( u_h ) ): recovered (projected) stress measure minus the one of the discrete solution
        return tFI->val()( 0 ) - this->eval_stress_value();
    }

    //------------------------------------------------------------------------------

    void
    IQI_Zienkiewicz_Zhu::compute_QI( Matrix< DDRMat >& aQI )
    {
        real tL2error = this->eval_error();

        // evaluate the QI: e^2 (integrated over the set it gives eta^2, see the class documentation)
        aQI = { { tL2error * tL2error } };
    }

    //------------------------------------------------------------------------------

    void
    IQI_Zienkiewicz_Zhu::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        real tL2error = this->eval_error();

        // evaluate the QI: accumulate w * e^2, i.e. eta^2 = int ( s* - s )^2 dOmega
        mSet->get_QI()( tQIIndex ) += aWStar * ( tL2error * tL2error );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
