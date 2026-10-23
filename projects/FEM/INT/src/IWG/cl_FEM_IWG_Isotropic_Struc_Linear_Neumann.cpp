/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Neumann.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "fn_FEM_Geometric_Sensitivity.hpp"

#include "fn_trans.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Neumann::IWG_Isotropic_Struc_Linear_Neumann()
    {
        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Traction" ]  = static_cast< uint >( IWG_Property_Type::TRACTION );
        mPropertyMap[ "Pressure" ]  = static_cast< uint >( IWG_Property_Type::PRESSURE );
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Linear_Neumann::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for the residual dof type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get traction load property
            const std::shared_ptr< Property >& tPropTraction =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get traction load property
            const std::shared_ptr< Property >& tPropPressure =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get sub-vector
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            if ( tPropTraction != nullptr )
            {
                tRes -= aWStar * ( tFI->N_trans() * tPropTraction->val() );
            }

            if ( tPropPressure != nullptr )
            {
                tRes -= aWStar * ( tFI->N_trans() * mNormal * tPropPressure->val() );
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Neumann::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------
        void
        IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for the residual dof type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get traction load property
            const std::shared_ptr< Property >& tPropTraction =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get traction load property
            const std::shared_ptr< Property >& tPropPressure =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute the Jacobian for indirect IWG dof dependencies through properties
            for ( uint iDOF = 0; iDOF < mRequestedLeaderGlobalDofTypes.size(); iDOF++ )
            {
                // get dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for the dof type
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if traction load depends on the dof type
                if ( tPropTraction != nullptr )
                {
                    if ( tPropTraction->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar * ( tFI->N_trans() * tPropTraction->dPropdDOF( tDofType ) );
                    }
                }

                // if pressure depends on the dof type
                if ( tPropPressure != nullptr )
                {
                    if ( tPropPressure->check_dof_dependency( tDofType ) )
                    {
                        // add contribution to Jacobian
                        tJac -= aWStar * ( tFI->N_trans() * mNormal * tPropPressure->dPropdDOF( tDofType ) );
                    }
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Neumann::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Neumann::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Neumann::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

        bool
        IWG_Isotropic_Struc_Linear_Neumann::supports_analytical_dRdp_geometry()
        {
            // the spatial variation of the loads must be known (MB-062)
            return has_known_spatial_variation( mLeaderProp( static_cast< uint >( IWG_Property_Type::TRACTION ) ) )    //
                && has_known_spatial_variation( mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) ) )    //
                && has_known_spatial_variation_nonzero( mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) ) );
        }

        //------------------------------------------------------------------------------

        std::shared_ptr< Property >
        IWG_Isotropic_Struc_Linear_Neumann::get_geometric_sensitivity_thickness_property()
        {
            return mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Neumann::compute_dRdx_geometry_integrand(
                real                        aWStar,
                Vector< Matrix< DDRMat > >& adRdx,
                Vector< Matrix< DDRMat > >& adRdn )
        {
            // get field interpolator for the residual dof type
            Field_Interpolator* tFI =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get traction load property
            const std::shared_ptr< Property >& tPropTraction =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::TRACTION ) );

            // get pressure load property
            const std::shared_ptr< Property >& tPropPressure =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::PRESSURE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // number of space dimensions and fields
            const uint tSpaceDim  = tFI->get_space_dim();
            const uint tNumFields = tFI->get_number_of_fields();

            // shape function gradients
            const Matrix< DDRMat >& tdNdx = tFI->dnNdxn( 1 );

            adRdx.resize( tSpaceDim );

            Matrix< DDRMat > tNk;
            Matrix< DDRMat > tdPropdx;
            for ( uint iDir = 0; iDir < tSpaceDim; iDir++ )
            {
                Matrix< DDRMat > tdNdxk = tdNdx( { iDir, iDir }, { 0, tdNdx.n_cols() - 1 } );
                build_vector_shape_function_matrix( tdNdxk, tNumFields, tNk );

                adRdx( iDir ).set_size( tFI->get_number_of_space_time_coefficients(), 1, 0.0 );

                // traction: d/dx_k ( N^T t ) = dN^T t + N^T dt/dx_k
                if ( tPropTraction != nullptr )
                {
                    eval_property_space_derivative( tPropTraction, iDir, tdPropdx );

                    adRdx( iDir ) -= aWStar * ( trans( tNk ) * tPropTraction->val() + tFI->N_trans() * tdPropdx );
                }

                // pressure: d/dx_k ( N^T n p ) = dN^T n p + N^T n dp/dx_k
                if ( tPropPressure != nullptr )
                {
                    eval_property_space_derivative( tPropPressure, iDir, tdPropdx );

                    adRdx( iDir ) -= aWStar * ( trans( tNk ) * mNormal * tPropPressure->val() + tFI->N_trans() * mNormal * tdPropdx );
                }
            }

            // normal dependence through the pressure only
            adRdn.clear();
            if ( tPropPressure != nullptr )
            {
                adRdn.resize( tSpaceDim );

                for ( uint iNormal = 0; iNormal < tSpaceDim; iNormal++ )
                {
                    Matrix< DDRMat > tUnitVector( tSpaceDim, 1, 0.0 );
                    tUnitVector( iNormal ) = 1.0;

                    adRdn( iNormal ) = -1.0 * aWStar * ( tFI->N_trans() * tUnitVector * tPropPressure->val() );
                }
            }
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
