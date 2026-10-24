/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Max_Stress.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Max_Stress.hpp"
#include "cl_FEM_CM_Struc_Linear.hpp"
#include "fn_FEM_Geometric_Sensitivity.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{

    //------------------------------------------------------------------------------

    IQI_Max_Stress::IQI_Max_Stress( enum Stress_Type aStressType )
            : IQI_Stress( aStressType )
    {
        MORIS_ERROR( aStressType != Stress_Type::STRESS_VECTOR,
                "IQI_Max_Stress - the stress vector cannot be aggregated; use a scalar stress measure." );

        // set FEM IQI type
        mFEMIQIType = fem::IQI_Type::MAX_STRESS;

        // set the property pointer cell size
        mLeaderProp.resize( static_cast< uint >( IQI_Property_Type::MAX_ENUM ), nullptr );

        // populate the property map
        mPropertyMap[ "ReferenceValue" ] = static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE );
        mPropertyMap[ "Exponent" ]       = static_cast< uint >( IQI_Property_Type::EXPONENT );
        mPropertyMap[ "Shift" ]          = static_cast< uint >( IQI_Property_Type::SHIFT );
        mPropertyMap[ "KSParameter" ]    = static_cast< uint >( IQI_Property_Type::KS_PARAMETER );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::compute_QI( Matrix< DDRMat >& aQI )
    {
        // evaluate the QI
        aQI = { { this->eval_integrand( this->eval_stress_value() ) } };
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * this->eval_integrand( this->eval_stress_value() );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::compute_dQIdu( real aWStar )
    {
        MORIS_ERROR( mStressType == Stress_Type::VON_MISES_STRESS,
                "IQI_Max_Stress::compute_dQIdu - analytic dQIdu is implemented for the von Mises stress only." );

        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the elasticity CM
        const std::shared_ptr< Constitutive_Model >& tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) );

        // stress measure and its derivative wrt the CM flux components (see IQI_Stress)
        real             tStressValue = this->eval_stress_value();
        Matrix< DDRMat > tdStressdFlux = this->eval_dStressdFlux();

        // derivative of the integrand wrt the stress measure
        real tdQIdStress = this->eval_dintegrand_dstress( tStressValue );

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        {
            // get the treated dof type
            const Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // the aggregation properties are not differentiated wrt the dofs
            for ( const std::shared_ptr< Property >& tProp : mLeaderProp )
            {
                MORIS_ERROR( tProp == nullptr || !tProp->check_dof_dependency( tDofType ),
                        "IQI_Max_Stress::compute_dQIdu - dof-dependent aggregation properties are not supported." );
            }

            // if the elasticity CM depends on the dof type
            if ( tCMElasticity->check_dof_dependency( tDofType ) )
            {
                mSet->get_residual()( tQIIndex )(
                        { tLeaderDepStartIndex, tLeaderDepStopIndex } ) +=
                        aWStar * tdQIdStress * trans( tCMElasticity->dFluxdDOF( tDofType, mFluxType ) ) * tdStressdFlux;
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::get_aggregation_parameters(
            real& aRefValue,
            real& aExponent,
            real& aShift,
            real& aKSParameter )
    {
        // check if properties set
        MORIS_ERROR( mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) ) != nullptr,
                "IQI_Max_Stress - no reference value set" );

        const std::shared_ptr< Property >& tPropKS = mLeaderProp( static_cast< uint >( IQI_Property_Type::KS_PARAMETER ) );

        // the exponent is required by the P-norm form only
        MORIS_ERROR( tPropKS != nullptr || mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr,
                "IQI_Max_Stress - no exponent set (needed unless KSParameter is set)" );

        aRefValue = mLeaderProp( static_cast< uint >( IQI_Property_Type::REFERENCE_VALUE ) )->val()( 0 );

        aExponent = 0.0;
        if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) ) != nullptr )
        {
            aExponent = mLeaderProp( static_cast< uint >( IQI_Property_Type::EXPONENT ) )->val()( 0 );
        }

        aShift = 1.0;
        if ( mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) ) != nullptr )
        {
            aShift = mLeaderProp( static_cast< uint >( IQI_Property_Type::SHIFT ) )->val()( 0 );
        }

        // KS parameter rho > 0 selects the Kreisselmeier-Steinhauser integrand; 0 = P-norm integrand
        aKSParameter = 0.0;
        if ( tPropKS != nullptr )
        {
            aKSParameter = tPropKS->val()( 0 );
            MORIS_ERROR( aKSParameter > 0.0, "IQI_Max_Stress - KSParameter must be positive, got %e", aKSParameter );
        }
    }

    //------------------------------------------------------------------------------

    real
    IQI_Max_Stress::eval_integrand( real aStressValue )
    {
        real tRefValue, tExponent, tShift, tKS;
        this->get_aggregation_parameters( tRefValue, tExponent, tShift, tKS );

        if ( tKS > 0.0 )
        {
            // shifted KS integrand exp( rho ( sigma / sigma_ref - shift ) ); with the default shift = 1 it is <= 1
            // wherever sigma <= sigma_ref, so the integral cannot overflow at a feasible design
            return std::exp( tKS * ( aStressValue / tRefValue - tShift ) );
        }

        // P-norm integrand ( sigma / sigma_ref - shift )^p
        return std::pow( aStressValue / tRefValue - tShift, tExponent );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Max_Stress::eval_dintegrand_dstress( real aStressValue )
    {
        real tRefValue, tExponent, tShift, tKS;
        this->get_aggregation_parameters( tRefValue, tExponent, tShift, tKS );

        if ( tKS > 0.0 )
        {
            return tKS / tRefValue * std::exp( tKS * ( aStressValue / tRefValue - tShift ) );
        }

        return tExponent / tRefValue * std::pow( aStressValue / tRefValue - tShift, tExponent - 1.0 );
    }

    //------------------------------------------------------------------------------

    bool
    IQI_Max_Stress::supports_analytical_dQIdp_geometry()
    {
        // only the von Mises stress is differentiated, and the spatial flux derivative below is that of the default flux
        if ( mStressType != Stress_Type::VON_MISES_STRESS || mFluxType != CM_Function_Type::DEFAULT )
        {
            return false;
        }

        // the spatial derivatives are implemented for the linear structural CMs
        CM_Struc_Linear* tCMStrucLinear = dynamic_cast< CM_Struc_Linear* >(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) ).get() );

        if ( tCMStrucLinear == nullptr || !tCMStrucLinear->has_spatial_strain_derivatives() )
        {
            return false;
        }

        // the spatial variation of the CM properties must be known, and the aggregation properties must be constant in
        // space (their spatial derivatives are not part of compute_dQIdx_geometry_integrand) (MB-062)
        bool tAggregationConstant = true;
        for ( const std::shared_ptr< Property >& tProp : mLeaderProp )
        {
            tAggregationConstant = tAggregationConstant && is_spatially_constant( tProp );
        }

        return tCMStrucLinear->has_spatial_constitutive_derivatives() && tAggregationConstant;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Max_Stress::compute_dQIdx_geometry_integrand(
            real              aWStar,
            Matrix< DDRMat >& adQIdx,
            Matrix< DDRMat >& adQIdn )
    {
        // get the elasticity CM
        CM_Struc_Linear* tCMElasticity = dynamic_cast< CM_Struc_Linear* >(
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST_LIN_ISO ) ).get() );

        MORIS_ASSERT( tCMElasticity != nullptr,
                "IQI_Max_Stress::compute_dQIdx_geometry_integrand - elasticity CM is not a CM_Struc_Linear" );

        // number of space dimensions
        const uint tSpaceDim = mLeaderFIManager->get_IP_geometry_interpolator()->get_number_of_space_dimensions();

        adQIdx.set_size( tSpaceDim, 1, 0.0 );
        adQIdn.set_size( 0, 0 );

        // stress measure, its derivative wrt the CM flux components and the derivative of the integrand
        // (P-norm or KS, see eval_integrand) wrt the stress measure; as in compute_dQIdu
        real             tStressValue  = this->eval_stress_value();
        Matrix< DDRMat > tdStressdFlux = this->eval_dStressdFlux();
        real             tdQIdStress   = this->eval_dintegrand_dstress( tStressValue );

        // constitutive matrix and strain
        const Matrix< DDRMat >& tConst  = tCMElasticity->constitutive();
        const Matrix< DDRMat >& tStrain = tCMElasticity->strain();

        Matrix< DDRMat > tdStraindx;
        Matrix< DDRMat > tdConstdx;
        for ( uint iDir = 0; iDir < tSpaceDim; iDir++ )
        {
            // derivative of the flux wrt x_k: dsigma = C deps + dC eps
            tCMElasticity->eval_dStraindx( iDir, tdStraindx );
            tCMElasticity->eval_dConstdx( iDir, tdConstdx );
            Matrix< DDRMat > tdFluxdx = tConst * tdStraindx + tdConstdx * tStrain;

            adQIdx( iDir ) = aWStar * tdQIdStress * dot( tdStressdFlux, tdFluxdx );
        }
    }

    //------------------------------------------------------------------------------

}    // namespace moris::fem
