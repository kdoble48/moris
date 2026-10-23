/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Bulk.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_CM_Struc_Linear.hpp"
#include "fn_FEM_Geometric_Sensitivity.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Bulk::IWG_Isotropic_Struc_Linear_Bulk()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Load" ]      = static_cast< uint >( IWG_Property_Type::LOAD );
        mPropertyMap[ "Bedding" ]   = static_cast< uint >( IWG_Property_Type::BEDDING );
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );
    }

    //------------------------------------------------------------------------------

    void IWG_Isotropic_Struc_Linear_Bulk::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for dof type
            Field_Interpolator* tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( MSI::Dof_Type::UX );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::LOAD ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            MORIS_ASSERT( !(tCMElasticity->get_plane_type() == Model_Type::AXISYMMETRIC
                    and tPropThickness == nullptr),
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_residual - must define axis of rotation "
                    "using IWG \"Thickness\" property as {{x1,y1},{x2,y2}} if using axisymmetric formulation");

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            tRes += aWStar * ( trans( tCMElasticity->testStrain() ) * tCMElasticity->flux() );

            // if body load
            if ( tPropLoad != nullptr )
            {
                // compute body load contribution
                tRes -= aWStar * ( trans( tDisplacementFI->N() ) * tPropLoad->val() );
            }

            // if bedding
            if ( tPropBedding != nullptr )
            {
                // compute body load contribution
                tRes += aWStar * ( trans( tDisplacementFI->N() ) * tDisplacementFI->val() * tPropBedding->val() );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_residual - Residual contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------
        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for given dof type
            Field_Interpolator * tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get body load property
            const std::shared_ptr< Property > & tPropLoad =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::LOAD ) );

            // get bedding property
            const std::shared_ptr< Property > & tPropBedding =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get elasticity CM
            const std::shared_ptr< Constitutive_Model > & tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get thickness property
            const std::shared_ptr< Property > & tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            MORIS_ASSERT( !(tCMElasticity->get_plane_type() == Model_Type::AXISYMMETRIC
                    and tPropThickness == nullptr),
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian - must define axis of rotation "
                    "using IWG \"Thickness\" property as {{x1,y1},{x2,y2}} if using axisymmetric formulation");

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= (tPropThickness!=nullptr) ? tPropThickness->val()(0) : 1;

            // get the number of leader dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over the leader dof dependencies
            for( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type > & tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if body load
                if ( tPropLoad != nullptr )
                {
                    // if property depends on the dof type
                    if ( tPropLoad->check_dof_dependency( tDofType ) )
                    {
                        // compute the contribution to Jacobian
                        tJac -= aWStar * ( trans( tDisplacementFI->N() ) * tPropLoad->dPropdDOF( tDofType ) );
                    }
                }

                // if bedding
                if ( tPropBedding != nullptr )
                {
                    if( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                    {
                        // if dof type is displacement, add bedding contribution
                        tJac += aWStar * ( trans( tDisplacementFI->N() ) *  tDisplacementFI->N() * tPropBedding->val()(0) );
                    }

                    // consider contributions from dependency of bedding parameter on DOFs
                    if ( tPropBedding->check_dof_dependency( tDofType ) )
                    {
                        tJac += aWStar * (
                                trans( tDisplacementFI->N() ) *  tDisplacementFI->val() * tPropBedding->dPropdDOF( tDofType ) );
                    }
                }

                // if constitutive model depends on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // compute the contribution to Jacobian
                    tJac += aWStar * ( trans( tCMElasticity->testStrain() ) * tCMElasticity->dFluxdDOF( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ) ,
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian - Jacobian contains NAN or INF, exiting!");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_jacobian_and_residual - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        void IWG_Isotropic_Struc_Linear_Bulk::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Bulk::compute_dRdp - This function does nothing.");
        }

        //------------------------------------------------------------------------------

        bool
        IWG_Isotropic_Struc_Linear_Bulk::supports_analytical_dRdp_geometry()
        {
            // get elasticity CM
            const std::shared_ptr< Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // the spatial derivatives are implemented for the linear structural CMs
            CM_Struc_Linear* tCMStrucLinear = dynamic_cast< CM_Struc_Linear* >( tCMElasticity.get() );

            if ( tCMStrucLinear == nullptr || !tCMStrucLinear->has_spatial_strain_derivatives() )
            {
                return false;
            }

            // the spatial variation of the CM and IWG properties must be known (MB-062)
            return tCMStrucLinear->has_spatial_constitutive_derivatives()                                             //
                && has_known_spatial_variation( mLeaderProp( static_cast< uint >( IWG_Property_Type::LOAD ) ) )       //
                && has_known_spatial_variation( mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) ) )    //
                && has_known_spatial_variation_nonzero( mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) ) );
        }

        //------------------------------------------------------------------------------

        std::shared_ptr< Property >
        IWG_Isotropic_Struc_Linear_Bulk::get_geometric_sensitivity_thickness_property()
        {
            return mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Bulk::compute_dRdx_geometry_integrand(
                real                        aWStar,
                Vector< Matrix< DDRMat > >& adRdx,
                Vector< Matrix< DDRMat > >& adRdn )
        {
            // get field interpolator for dof type
            Field_Interpolator* tDisplacementFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get body load property
            const std::shared_ptr< Property >& tPropLoad =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::LOAD ) );

            // get bedding property
            const std::shared_ptr< Property >& tPropBedding =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::BEDDING ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // get elasticity CM
            CM_Struc_Linear* tCMElasticity = dynamic_cast< CM_Struc_Linear* >(
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) ).get() );

            MORIS_ASSERT( tCMElasticity != nullptr,
                    "IWG_Isotropic_Struc_Linear_Bulk::compute_dRdx_geometry_integrand - elasticity CM is not a CM_Struc_Linear" );

            // multiplying aWStar by user defined thickness
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // number of space dimensions and fields
            const uint tSpaceDim  = tDisplacementFI->get_space_dim();
            const uint tNumFields = tDisplacementFI->get_number_of_fields();

            // unperturbed quantities
            const Matrix< DDRMat >& tTestStrain = tCMElasticity->testStrain();
            const Matrix< DDRMat >& tFlux       = tCMElasticity->flux();
            const Matrix< DDRMat >& tConst      = tCMElasticity->constitutive();
            const Matrix< DDRMat >& tdNdx       = tDisplacementFI->dnNdxn( 1 );

            // vectorized displacement coefficients
            Matrix< DDRMat > tUHat;
            vectorize_coefficients( tDisplacementFI->get_coeff(), tUHat );

            adRdx.resize( tSpaceDim );
            adRdn.clear();

            Matrix< DDRMat > tdBdx;
            Matrix< DDRMat > tdConstdx;
            Matrix< DDRMat > tdPropdx;
            Matrix< DDRMat > tNk;
            for ( uint iDir = 0; iDir < tSpaceDim; iDir++ )
            {
                // derivative of the strain-displacement matrix and the flux wrt x_k: d sigma = C dB u + dC eps
                tCMElasticity->eval_dTestStraindx( iDir, tdBdx );
                tCMElasticity->eval_dConstdx( iDir, tdConstdx );
                Matrix< DDRMat > tdFluxdx = tConst * ( tdBdx * tUHat ) + tdConstdx * tCMElasticity->strain();

                // d/dx_k ( B^T sigma ) = dB^T sigma + B^T C dB u
                adRdx( iDir ) = aWStar * ( trans( tdBdx ) * tFlux + trans( tTestStrain ) * tdFluxdx );

                // derivative of the shape function matrix wrt x_k
                if ( tPropLoad != nullptr || tPropBedding != nullptr )
                {
                    Matrix< DDRMat > tdNdxk = tdNdx( { iDir, iDir }, { 0, tdNdx.n_cols() - 1 } );
                    build_vector_shape_function_matrix( tdNdxk, tNumFields, tNk );
                }

                // if body load: d/dx_k ( N^T f ) = dN^T f + N^T df/dx_k
                if ( tPropLoad != nullptr )
                {
                    eval_property_space_derivative( tPropLoad, iDir, tdPropdx );

                    adRdx( iDir ) -= aWStar * ( trans( tNk ) * tPropLoad->val() + tDisplacementFI->N_trans() * tdPropdx );
                }

                // if bedding: d/dx_k ( N^T N u k ) = ( dN^T N + N^T dN ) u k + N^T N u dk/dx_k
                if ( tPropBedding != nullptr )
                {
                    eval_property_space_derivative( tPropBedding, iDir, tdPropdx );

                    adRdx( iDir ) += aWStar * tPropBedding->val()( 0 ) *    //
                                     ( trans( tNk ) * tDisplacementFI->val() + tDisplacementFI->N_trans() * ( tNk * tUHat ) );
                    adRdx( iDir ) += aWStar * tdPropdx( 0 ) * ( tDisplacementFI->N_trans() * tDisplacementFI->val() );
                }
            }
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
