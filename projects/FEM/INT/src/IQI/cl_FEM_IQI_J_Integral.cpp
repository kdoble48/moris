/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_J_Integral.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_J_Integral.hpp"
#include "fn_dot.hpp"
#include "fn_trans.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_J_Integral::IQI_J_Integral()
    {
        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO );
    }

    //------------------------------------------------------------------------------

    void IQI_J_Integral::compute_QI( Matrix< DDRMat > &aQI )
    {
        /* evaluate the QI, crack growing straight ahead in x_1:
         * G = J = int_{\Gamma} ( W n_1 - t_i \frac{\partial u_i}{\partial x_1} ),
         * with W = 1/2 sigma : epsilon, t = sigma n and n the outward facet normal
         */

        // get the elasticity CM
        const std::shared_ptr< Constitutive_Model > &tCM =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

        // get the displacement field interpolator
        Field_Interpolator *tFIDispl =
                mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

        // strain energy density (the Voigt strain carries the engineering shear strain)
        real tW = 0.5 * dot( tCM->flux(), tCM->strain() );

        // derivatives of the displacement components in x_1
        Matrix< DDRMat > tdUdX1 = trans( tFIDispl->gradx( 1 ).get_row( 0 ) );

        aQI = { { tW * mNormal( 0 ) - dot( tCM->traction( mNormal ), tdUdX1 ) } };
    }

    //------------------------------------------------------------------------------

    void IQI_J_Integral::compute_QI( real aWStar )
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
     // end moris namespace
