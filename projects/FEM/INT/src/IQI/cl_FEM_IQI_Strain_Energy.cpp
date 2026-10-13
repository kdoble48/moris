/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Strain_Energy.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Strain_Energy.hpp"
#include "cl_FEM_CM_Struc_Linear.hpp"
#include "fn_FEM_Geometric_Sensitivity.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Strain_Energy::IQI_Strain_Energy()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::STRAIN_ENERGY;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Bedding" ]   = static_cast< uint >( IWG_Property_Type::BEDDING );
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Elast" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST );
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_QI( Matrix< DDRMat > &aQI )
    {
        // check the point is inside the bounded box
        if ( !this->is_within_box_bounds() )
        {
            return;
        }

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

        // evaluate the QI
        aQI = 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain();

        // if bedding
        if ( tPropBedding != nullptr )
        {
            // get field interpolator for displacements
            Field_Interpolator *tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // check that field interpolator exists
            MORIS_ASSERT( tDisplacementFI != nullptr,
                    "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist." );

            // compute body load contribution
            aQI += 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val();
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_QI( real aWStar )
    {
        // check the point is inside the bounded box
        if ( !this->is_within_box_bounds() )
        {
            return;
        }

        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

        // get thickness property
        const std::shared_ptr< Property > &tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * ( 0.5 * trans( tCMElasticity->flux() ) * tCMElasticity->strain() );

        // if bedding
        if ( tPropBedding != nullptr )
        {
            // get field interpolator for displacements
            Field_Interpolator *tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // check that field interpolator exists
            MORIS_ASSERT( tDisplacementFI != nullptr,
                    "IQI_Strain_Energy::compute_QI - field interpolator for dof type UX does not exist." );

            // compute body load contribution
            mSet->get_QI()( tQIIndex ) += aWStar * ( 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->val() );
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_dQIdu( real aWStar )
    {
        // check the point is inside the bounded box
        if ( !this->is_within_box_bounds() )
        {
            return;
        }

        // get the column index to assemble in residual
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

        // get thickness property
        const std::shared_ptr< Property > &tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute dQIdu for indirect dof dependencies
        for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        {
            // get the treated dof type
            Vector< MSI::Dof_Type > &tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // if elasticity CM depends on dof type
            if ( tCMElasticity->check_dof_dependency( tDofType ) )
            {
                // compute dQIdu
                mSet->get_residual()( tQIIndex )(
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * 0.5 * ( trans( tCMElasticity->dFluxdDOF( tDofType ) ) * tCMElasticity->strain() + trans( tCMElasticity->dStraindDOF( tDofType ) ) * tCMElasticity->flux() );
            }

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // get field interpolator for displacements
                Field_Interpolator *tDisplacementFI =
                        mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

                // compute bedding contribution - displacements
                if ( tDofType( 0 ) == MSI::Dof_Type::UX )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val() );
                }

                if ( tPropBedding->check_dof_dependency( tDofType ) )
                {
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex } ) += aWStar * ( 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( tDofType ) );
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void IQI_Strain_Energy::compute_dQIdu(
            Vector< MSI::Dof_Type > &aDofType,
            Matrix< DDRMat >        &adQIdu )
    {
        // check the point is inside the bounded box
        if ( !this->is_within_box_bounds() )
        {
            return;
        }

        // get the elasticity CM
        std::shared_ptr< Constitutive_Model > &tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

        // initialize derivative
        adQIdu.fill( 0.0 );

        // if elasticity CM depends on dof type
        if ( tCMElasticity->check_dof_dependency( aDofType ) )
        {
            // compute dQIdu
            adQIdu = 0.5 * ( trans( tCMElasticity->dFluxdDOF( aDofType ) ) * tCMElasticity->strain() + trans( tCMElasticity->dStraindDOF( aDofType ) ) * tCMElasticity->flux() );
        }

        // if bedding
        if ( tPropBedding != nullptr )
        {
            // get field interpolator for displacements
            Field_Interpolator *tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // compute bedding contribution - displacements
            if ( aDofType( 0 ) == MSI::Dof_Type::UX )
            {
                adQIdu += trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val();
            }

            if ( tPropBedding->check_dof_dependency( aDofType ) )
            {
                adQIdu += 0.5 * trans( tDisplacementFI->val() ) * tDisplacementFI->val() * tPropBedding->dPropdDOF( aDofType );
            }
        }
    }

    //------------------------------------------------------------------------------

    bool IQI_Strain_Energy::is_within_box_bounds()
    {
        // check if the box bounds are empty then skip
        if ( mParameters.empty() )
        {
            return true;
        }

        // if the box bounds are not empty then check if it is inside the box
        else
        {
            // get the coordinate
            const Matrix< DDRMat > &tGaussPoint = mLeaderFIManager->get_IG_geometry_interpolator()->valx();

            // check if the calculation point coordinates are more then lower corner of the box
            bool tLowerBound = std::equal( mParameters( 0 ).begin(), mParameters( 0 ).end(), tGaussPoint.begin(), tGaussPoint.end(), []( real aA, real aB ) -> bool { return aA < aB; } );

            // check if the calculation point coordinates are less then upper corner of the box
            bool tUpperBound = std::equal( tGaussPoint.begin(), tGaussPoint.end(), mParameters( 1 ).begin(), mParameters( 1 ).end(), []( real aA, real aB ) -> bool { return aA < aB; } );

            // combine the two bounds that satisfy both
            return tUpperBound and tLowerBound;
        }
    }

    //------------------------------------------------------------------------------
    //------------------------------------------------------------------------------

    bool
    IQI_Strain_Energy::supports_analytical_dQIdp_geometry()
    {
        // the spatial derivatives are implemented for the linear structural CMs
        const CM_Struc_Linear *tCMStrucLinear = dynamic_cast< const CM_Struc_Linear * >(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) ).get() );

        return tCMStrucLinear != nullptr && tCMStrucLinear->has_spatial_strain_derivatives();
    }

    //------------------------------------------------------------------------------

    void
    IQI_Strain_Energy::compute_dQIdx_geometry_integrand(
            real              aWStar,
            Matrix< DDRMat > &adQIdx,
            Matrix< DDRMat > &adQIdn )
    {
        // number of space dimensions
        const uint tSpaceDim = mLeaderFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions();

        adQIdx.set_size( tSpaceDim, 1, 0.0 );
        adQIdn.set_size( 0, 0 );

        // outside the box bounds the QI contribution is zero (box bounds are held fixed in space)
        if ( !this->is_within_box_bounds() )
        {
            return;
        }

        // get the elasticity CM
        CM_Struc_Linear *tCMElasticity = dynamic_cast< CM_Struc_Linear * >(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) ).get() );

        MORIS_ASSERT( tCMElasticity != nullptr,
                "IQI_Strain_Energy::compute_dQIdx_geometry_integrand - elasticity CM is not a CM_Struc_Linear" );

        // get bedding property
        const std::shared_ptr< Property > &tPropBedding =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

        // get thickness property
        const std::shared_ptr< Property > &tPropThickness =
                mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

        // multiplying aWStar by user defined thickness
        aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

        // unperturbed flux and strain, constitutive matrix
        const Matrix< DDRMat > &tFlux   = tCMElasticity->flux();
        const Matrix< DDRMat > &tStrain = tCMElasticity->strain();
        const Matrix< DDRMat > &tConst  = tCMElasticity->constitutive();

        // displacement field interpolator (bedding)
        Field_Interpolator *tDisplacementFI = nullptr;
        Matrix< DDRMat >    tUHat;
        if ( tPropBedding != nullptr )
        {
            tDisplacementFI = mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );
            vectorize_coefficients( tDisplacementFI->get_coeff(), tUHat );
        }

        Matrix< DDRMat > tdStraindx;
        Matrix< DDRMat > tNk;
        for ( uint iDir = 0; iDir < tSpaceDim; iDir++ )
        {
            // d/dx_k ( 1/2 sigma^T eps ) = 1/2 ( C deps^T eps + sigma^T deps )
            tCMElasticity->eval_dStraindx( iDir, tdStraindx );
            Matrix< DDRMat > tdFluxdx = tConst * tdStraindx;

            adQIdx( iDir ) = aWStar * 0.5 * ( dot( tdFluxdx, tStrain ) + dot( tFlux, tdStraindx ) );

            // bedding: d/dx_k ( 1/2 k u^T u ) = k u^T du/dx_k
            if ( tPropBedding != nullptr )
            {
                const Matrix< DDRMat > &tdNdx  = tDisplacementFI->dnNdxn( 1 );
                Matrix< DDRMat >        tdNdxk = tdNdx( { iDir, iDir }, { 0, tdNdx.n_cols() - 1 } );
                build_vector_shape_function_matrix( tdNdxk, tDisplacementFI->get_number_of_fields(), tNk );

                Matrix< DDRMat > tdUdx = tNk * tUHat;

                adQIdx( iDir ) += aWStar * tPropBedding->val()( 0 ) * dot( tDisplacementFI->val(), tdUdx );
            }
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
