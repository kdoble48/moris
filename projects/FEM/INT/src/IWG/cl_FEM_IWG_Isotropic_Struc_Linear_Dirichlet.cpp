/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.cpp
 *
 */

#include "cl_FEM_IWG_Isotropic_Struc_Linear_Dirichlet.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_Set.hpp"
#include "cl_FEM_CM_Struc_Linear.hpp"
#include "fn_FEM_Geometric_Sensitivity.hpp"

#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IWG_Isotropic_Struc_Linear_Dirichlet::IWG_Isotropic_Struc_Linear_Dirichlet( sint aBeta )
    {
        // sign for symmetric/unsymmetric Nitsche
        mBeta = aBeta;

        // set size for the property pointer cell
        mLeaderProp.resize( static_cast< uint >( IWG_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "Dirichlet" ] = static_cast< uint >( IWG_Property_Type::DIRICHLET );
        mPropertyMap[ "Select" ]    = static_cast< uint >( IWG_Property_Type::SELECT );
        mPropertyMap[ "Thickness" ] = static_cast< uint >( IWG_Property_Type::THICKNESS );

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IWG_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "ElastLinIso" ] = static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO );

        // set size for the stabilization parameter pointer cell
        mStabilizationParam.resize( static_cast< uint >( IWG_Stabilization_Type::MAX_ENUM ), nullptr );

        // populate the stabilization map
        mStabilizationMap[ "DirichletNitsche" ] = static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE );
    }

    //------------------------------------------------------------------------------

    void
    IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual( real aWStar )
    {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tFIDispl = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tFIDispl->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();

                // skip computing residual if projection matrix is zero
                if ( norm( tM ) < MORIS_REAL_EPS ) return;
            }

            // get the imposed displacement property
            const std::shared_ptr< Property >& tPropDirichlet =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get CM for elasticity
            const std::shared_ptr< Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get SP for Nitsche
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute jump
            Matrix< DDRMat > tJump = tFIDispl->val() - tPropDirichlet->val();

            // compute the residual
            mSet->get_residual()( 0 )(
                    { tLeaderResStartIndex, tLeaderResStopIndex } ) +=
                    aWStar * ( -tFIDispl->N_trans() * tM * tCMElasticity->traction( mNormal )                                //
                               + mBeta * tCMElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) ) * tM * tJump    //
                               + tSPNitsche->val()( 0 ) * tFIDispl->N_trans() * tM * tJump );

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_residual()( 0 ) ),
                    "IWG_Isotropic_Struc_Linear_Dirichlet::compute_residual - Residual contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian( real aWStar )
        {
#ifdef MORIS_HAVE_DEBUG
            // check leader field interpolators, properties and constitutive models
            this->check_field_interpolators();
#endif

            // get leader index for residual dof type (here displacement), indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( mResidualDofType( 0 )( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderResStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderResStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // get field interpolator for residual dof type
            Field_Interpolator* tFIDispl =
                    mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                // get number of fields which should equal spatial dimension
                const uint tSpaceDim = tFIDispl->get_number_of_fields();

                // set selection matrix as identity
                eye( tSpaceDim, tSpaceDim, tM );
            }
            else
            {
                tM = tPropSelect->val();

                // skip computing Jacobian if projection matrix is null
                if ( norm( tM ) < MORIS_REAL_EPS ) return;
            }

            // get the imposed displacement property
            const std::shared_ptr< Property >& tPropDirichlet =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get CM for elasticity
            const std::shared_ptr< Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // get SP for Nitsche
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness (2*pi*r for axisymmetric)
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // compute jump
            Matrix< DDRMat > tJump = tFIDispl->val() - tPropDirichlet->val();

            // get number of dof dependencies
            uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();
            for ( uint iDOF = 0; iDOF < tNumDofDependencies; iDOF++ )
            {
                // get the dof type
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
                    // compute the jacobian for direct dof dependencies
                    tJac += aWStar * (                                                                                                  //
                                    mBeta * tCMElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) ) * tM * tFIDispl->N()    //
                                    + tSPNitsche->val()( 0 ) * tFIDispl->N_trans() * tM * tFIDispl->N() );
                }

                // if dependency on the dof type
                if ( tPropDirichlet->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac -= aWStar * (                                                                                                                          //
                                    mBeta * tCMElasticity->testTraction_trans( mNormal, mResidualDofType( 0 ) ) * tM * tPropDirichlet->dPropdDOF( tDofType )    //
                                    + tSPNitsche->val()( 0 ) * tFIDispl->N_trans() * tM * tPropDirichlet->dPropdDOF( tDofType ) );
                }

                // if dependency on the dof type
                if ( tCMElasticity->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * (                                                                               //
                                    -tFIDispl->N_trans() * tM * tCMElasticity->dTractiondDOF( tDofType, mNormal )    //
                                    + mBeta * tCMElasticity->dTestTractiondDOF( tDofType, mNormal, tM * tJump, mResidualDofType( 0 ) ) );
                }

                // if dependency on the dof type
                if ( tSPNitsche->check_dof_dependency( tDofType ) )
                {
                    // add contribution to jacobian
                    tJac += aWStar * ( tFIDispl->N_trans() * tM * tJump * tSPNitsche->dSPdLeaderDOF( tDofType ) );
                }
            }

            // check for nan, infinity
            MORIS_ASSERT( isfinite( mSet->get_jacobian() ),
                    "IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian - Jacobian contains NAN or INF, exiting!" );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_jacobian_and_residual - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdp( real aWStar )
        {
            MORIS_ERROR( false, "IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdp - This function does nothing." );
        }

        //------------------------------------------------------------------------------

        bool
        IWG_Isotropic_Struc_Linear_Dirichlet::supports_analytical_dRdp_geometry()
        {
            // get CM for elasticity
            const std::shared_ptr< Constitutive_Model >& tCMElasticity =
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) );

            // the spatial derivatives are implemented for the linear structural CMs
            CM_Struc_Linear* tCMStrucLinear = dynamic_cast< CM_Struc_Linear* >( tCMElasticity.get() );

            if ( tCMStrucLinear == nullptr || !tCMStrucLinear->has_spatial_strain_derivatives() )
            {
                return false;
            }

            // get SP for Nitsche
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // the spatial variation of the CM, SP and IWG properties must be known (MB-062)
            return tCMStrucLinear->has_spatial_constitutive_derivatives()                                               //
                && tSPNitsche != nullptr && tSPNitsche->has_spatial_derivatives()                                       //
                && has_known_spatial_variation( mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) ) )    //
                && is_spatially_constant( mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) ) )             //
                && has_known_spatial_variation_nonzero( mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) ) );
        }

        //------------------------------------------------------------------------------

        std::shared_ptr< Property >
        IWG_Isotropic_Struc_Linear_Dirichlet::get_geometric_sensitivity_thickness_property()
        {
            return mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );
        }

        //------------------------------------------------------------------------------

        void
        IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdx_geometry_integrand(
                real                        aWStar,
                Vector< Matrix< DDRMat > >& adRdx,
                Vector< Matrix< DDRMat > >& adRdn )
        {
            // get field interpolator for residual dof type
            Field_Interpolator* tFIDispl = mLeaderFIManager->get_field_interpolators_for_type( mResidualDofType( 0 )( 0 ) );

            // number of space dimensions and fields
            const uint tSpaceDim  = tFIDispl->get_space_dim();
            const uint tNumFields = tFIDispl->get_number_of_fields();

            // get the selection matrix property
            const std::shared_ptr< Property >& tPropSelect =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::SELECT ) );

            // set a default selection matrix if needed
            Matrix< DDRMat > tM;
            if ( tPropSelect == nullptr )
            {
                eye( tNumFields, tNumFields, tM );
            }
            else
            {
                tM = tPropSelect->val();
            }

            // initialize outputs (the residual is zero for a zero selection matrix, cf. compute_residual)
            adRdx.resize( tSpaceDim );
            adRdn.resize( tSpaceDim );
            for ( uint iDir = 0; iDir < tSpaceDim; iDir++ )
            {
                adRdx( iDir ).set_size( tFIDispl->get_number_of_space_time_coefficients(), 1, 0.0 );
                adRdn( iDir ).set_size( tFIDispl->get_number_of_space_time_coefficients(), 1, 0.0 );
            }

            if ( tPropSelect != nullptr && norm( tM ) < MORIS_REAL_EPS )
            {
                return;
            }

            // get the imposed displacement property
            const std::shared_ptr< Property >& tPropDirichlet =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::DIRICHLET ) );

            // get CM for elasticity
            CM_Struc_Linear* tCMElasticity = dynamic_cast< CM_Struc_Linear* >(
                    mLeaderCM( static_cast< uint >( IWG_Constitutive_Type::ELAST_LIN_ISO ) ).get() );

            MORIS_ASSERT( tCMElasticity != nullptr,
                    "IWG_Isotropic_Struc_Linear_Dirichlet::compute_dRdx_geometry_integrand - elasticity CM is not a CM_Struc_Linear" );

            // get SP for Nitsche (its cluster measure dependence is added by the caller)
            const std::shared_ptr< Stabilization_Parameter >& tSPNitsche =
                    mStabilizationParam( static_cast< uint >( IWG_Stabilization_Type::DIRICHLET_NITSCHE ) );

            // get thickness property
            const std::shared_ptr< Property >& tPropThickness =
                    mLeaderProp( static_cast< uint >( IWG_Property_Type::THICKNESS ) );

            // multiplying aWStar by user defined thickness
            aWStar *= ( tPropThickness != nullptr ) ? tPropThickness->val()( 0 ) : 1;

            // unperturbed quantities
            const Matrix< DDRMat >& tTestStrain = tCMElasticity->testStrain();
            const Matrix< DDRMat >& tFlux       = tCMElasticity->flux();
            const Matrix< DDRMat >& tConst      = tCMElasticity->constitutive();
            const Matrix< DDRMat >& tdNdx       = tFIDispl->dnNdxn( 1 );
            const real              tGamma      = tSPNitsche->val()( 0 );

            // jump
            Matrix< DDRMat > tJump = tFIDispl->val() - tPropDirichlet->val();

            // flattened normal: traction = F sigma, test traction = F C B
            Matrix< DDRMat > tFlatNormal;
            tCMElasticity->get_flattened_normal( mNormal, tFlatNormal );

            // vectorized displacement coefficients
            Matrix< DDRMat > tUHat;
            vectorize_coefficients( tFIDispl->get_coeff(), tUHat );

            // derivatives wrt the evaluation point location
            Matrix< DDRMat > tdBdx;
            Matrix< DDRMat > tdConstdx;
            Matrix< DDRMat > tdDirichletdx;
            Matrix< DDRMat > tdGammadx;
            Matrix< DDRMat > tNk;
            for ( uint iDir = 0; iDir < tSpaceDim; iDir++ )
            {
                // d sigma = C dB u + dC eps
                tCMElasticity->eval_dTestStraindx( iDir, tdBdx );
                tCMElasticity->eval_dConstdx( iDir, tdConstdx );
                Matrix< DDRMat > tdFluxdx = tConst * ( tdBdx * tUHat ) + tdConstdx * tCMElasticity->strain();

                Matrix< DDRMat > tdNdxk = tdNdx( { iDir, iDir }, { 0, tdNdx.n_cols() - 1 } );
                build_vector_shape_function_matrix( tdNdxk, tNumFields, tNk );

                // derivative of the jump: dN u - d u_D / dx_k
                eval_property_space_derivative( tPropDirichlet, iDir, tdDirichletdx );
                Matrix< DDRMat > tdJumpdx = tNk * tUHat - tdDirichletdx;

                // derivative of the Nitsche parameter at fixed cluster measures
                tSPNitsche->eval_dSPdx( iDir, tdGammadx );

                adRdx( iDir ) = aWStar * (                                                                                //
                                        -trans( tNk ) * tM * tFlatNormal * tFlux                                          //
                                        - tFIDispl->N_trans() * tM * tFlatNormal * tdFluxdx                               //
                                        + mBeta * trans( tFlatNormal * tConst * tdBdx ) * tM * tJump                      //
                                        + mBeta * trans( tFlatNormal * tConst * tTestStrain ) * tM * tdJumpdx             //
                                        + mBeta * trans( tFlatNormal * tdConstdx * tTestStrain ) * tM * tJump             //
                                        + tGamma * ( trans( tNk ) * tM * tJump + tFIDispl->N_trans() * tM * tdJumpdx )    //
                                        + tdGammadx( 0 ) * tFIDispl->N_trans() * tM * tJump );
            }

            // derivatives wrt the normal components: F is linear in the normal
            for ( uint iNormal = 0; iNormal < tSpaceDim; iNormal++ )
            {
                Matrix< DDRMat > tUnitVector( tSpaceDim, 1, 0.0 );
                tUnitVector( iNormal ) = 1.0;

                Matrix< DDRMat > tFlatUnit;
                tCMElasticity->get_flattened_normal( tUnitVector, tFlatUnit );

                adRdn( iNormal ) = aWStar * (                                               //
                                           -tFIDispl->N_trans() * tM * tFlatUnit * tFlux    //
                                           + mBeta * trans( tFlatUnit * tConst * tTestStrain ) * tM * tJump );
            }
        }

        //------------------------------------------------------------------------------
}    // namespace moris::fem
