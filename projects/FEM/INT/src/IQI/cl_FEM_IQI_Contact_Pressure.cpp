/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Contact_Pressure.cpp
 *
 */

#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Contact_Pressure.hpp"
#include "cl_Matrix.hpp"
#include "moris_typedefs.hpp"
#include "fn_assert.hpp"
#include "fn_norm.hpp"
#include "fn_dot.hpp"
#include "cl_FEM_CM_Struc_Nonlinear_Isotropic.hpp"

#include <iostream>
#include <memory>

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Contact_Pressure::IQI_Contact_Pressure( bool aInReferenceConfiguration )
            : mInReferenceConfiguration( aInReferenceConfiguration )
    {
        mFEMIQIType = fem::IQI_Type::JUMP_TRACTION;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );
        mFollowerCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "TractionCM" ] = (uint)IQI_Constitutive_Type::TRACTION_CM;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Contact_Pressure::compute_QI( Matrix< DDRMat >& aQI )
    {
        // get the constitutive model for computing the fluxes/tractions
        const std::shared_ptr< Constitutive_Model >& tCMLeader = mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::TRACTION_CM ) );

        // no facet normal (e.g. the IQI requested on a bulk set for output): the pressure is not defined
        Field_Interpolator* tLeaderDofs = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
        if ( mNormal.numel() != tLeaderDofs->get_number_of_fields() )
        {
            aQI = { { 0.0 } };
            return;
        }

        // a small-strain CM has one stress and one configuration: p = N . sigma N with the reference normal
        if ( dynamic_cast< CM_Struc_Nonlinear_Isotropic* >( tCMLeader.get() ) == nullptr )
        {
            aQI = { dot( tCMLeader->traction( mNormal ), mNormal ) };
            return;
        }

        if ( mInReferenceConfiguration )
        {
            // nominal pressure: N . P N with the undeformed (IG) normal N
            aQI = { dot( tCMLeader->traction( mNormal, CM_Function_Type::PK1 ), mNormal ) };
        }
        else
        {
            // true pressure: n . sigma n with the Cauchy stress and the deformed normal n
            Geometry_Interpolator* tLeaderGeometry = mLeaderFIManager->get_IG_geometry_interpolator();

            Matrix< DDRMat > tNormal = tLeaderGeometry->get_normal_current( tLeaderDofs );

            // Cauchy stress in Voigt notation: ( xx, yy, xy ) in 2D, ( xx, yy, zz, yz, xz, xy ) in 3D
            // (the nonlinear CMs implement the Cauchy flux, not all of them the Cauchy traction)
            const Matrix< DDRMat >& tSigma = tCMLeader->flux( CM_Function_Type::CAUCHY );

            real tPressure = 0.0;
            if ( tNormal.numel() == 2 )
            {
                tPressure = tSigma( 0 ) * tNormal( 0 ) * tNormal( 0 ) + tSigma( 1 ) * tNormal( 1 ) * tNormal( 1 )
                          + 2.0 * tSigma( 2 ) * tNormal( 0 ) * tNormal( 1 );
            }
            else
            {
                tPressure = tSigma( 0 ) * tNormal( 0 ) * tNormal( 0 ) + tSigma( 1 ) * tNormal( 1 ) * tNormal( 1 )
                          + tSigma( 2 ) * tNormal( 2 ) * tNormal( 2 )
                          + 2.0 * ( tSigma( 3 ) * tNormal( 1 ) * tNormal( 2 ) + tSigma( 4 ) * tNormal( 0 ) * tNormal( 2 ) + tSigma( 5 ) * tNormal( 0 ) * tNormal( 1 ) );
            }

            aQI = { { tPressure } };
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Contact_Pressure::compute_QI( real aWStar )
    {
        Matrix< DDRMat > tMat( 1, 1 );
        this->compute_QI( tMat );

        // add the contribution
        sint tQIIndex = mSet->get_QI_assembly_index( mName );
        mSet->get_QI()( tQIIndex ) += aWStar * tMat;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Contact_Pressure::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Contact_Pressure::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        MORIS_ERROR( false, "Not Implemented for pseudo error for double sided set " );
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
