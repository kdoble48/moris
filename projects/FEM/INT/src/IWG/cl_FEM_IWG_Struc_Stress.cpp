/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Struc_Stress.cpp
 *
 */

#include "cl_FEM_IWG_Struc_Stress.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "fn_FEM_Stress_Measures.hpp"
// LINALG/src
#include "fn_trans.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Struc_Stress::IWG_Struc_Stress( enum Stress_Type aStressType )
    {
        // the projected measure carries no component index, so only the von Mises stress is supported
        MORIS_ERROR( aStressType == Stress_Type::VON_MISES_STRESS,
                "IWG_Struc_Stress - only the von Mises stress is supported." );

        // set stress type
        mStressType = aStressType;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Struc_Stress::compute_residual( real aWStar )
    {
        // check leader field interpolators
#ifdef MORIS_HAVE_DEBUG
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 1 );

            // get residual dof type field interpolator
            Field_Interpolator* tFISig = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // stress measure of the discrete solution (see fn_FEM_Stress_Measures)
            real tStressVal = eval_stress_measure(
                    mStressType,
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) )->flux(),
                    0,
                    0 );

            // get sub-matrix
            auto tRes = mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } );

            // compute the residual
            tRes += aWStar * (                                     //
                            tFISig->N_trans() * tFISig->val() -    //
                            tFISig->N_trans() * tStressVal );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Diffusion_Bulk::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 1 );

            // get field interpolator for a given dof type
            Field_Interpolator* tFISig = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get elasticity CM
            std::shared_ptr< Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // get the number of leader dof type dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

            // loop over leader dof type dependencies
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the treated dof type
                const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDOF );

                // get the index for dof type, indices for assembly
                sint tDofDepIndex         = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
                uint tLeaderDepStartIndex = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 0 );
                uint tLeaderDepStopIndex  = mSet->get_jac_dof_assembly_map()( tLeaderDofIndex )( tDofDepIndex, 1 );

                // get sub-matrix
                auto tJac = mSet->get_jacobian()(
                        { tLeaderResStartIndex, tLeaderResStopIndex },
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } );

                // if dof type is residual dof type
                if ( tDofType( 0 ) == mResidualDofType( 0 )( 0 ) )
                {
                    tJac += aWStar * ( tFISig->N_trans() * tFISig->N() );
                }

                // if constitutive model has dependency on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // derivative of the stress measure wrt the dofs: d s / d flux ( flux layout ) contracted with d flux / d dof
                    Matrix< DDRMat > tDStressVal =
                            trans( eval_dstress_measure_dflux( mStressType, tCMElasticity->flux() ) ) * tCMElasticity->dFluxdDOF( tDofType );

                    // compute the Jacobian
                    tJac += ( -1.0 ) * aWStar * ( tFISig->N_trans() * tDStressVal );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Struc_Stress::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Diffusion_Bulk::compute_jacobian_and_residual - Not implemented." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Struc_Stress::compute_dRdp( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            MORIS_ERROR( false, "IWG_Isotropic_Spatial_Diffusion_Bulk::compute_dRdp - Not implemented." );
        }

        //------------------------------------------------------------------------------

}    // namespace moris::fem
