/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Topological_Derivative.cpp
 *
 */

#include <cmath>

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Topological_Derivative.hpp"
#include "cl_FEM_Constitutive_Model.hpp"
#include "cl_FEM_Property.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Topological_Derivative::IQI_Topological_Derivative()
    {
        // set fem IQI type
        mFEMIQIType = fem::IQI_Type::TOPOLOGICAL_DERIVATIVE;

        // set size for the constitutive model pointer cell
        mLeaderCM.resize( static_cast< uint >( IQI_Constitutive_Type::MAX_ENUM ), nullptr );

        // populate the constitutive map
        mConstitutiveMap[ "Elast" ] = static_cast< uint >( IQI_Constitutive_Type::ELAST );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Topological_Derivative::set_parameters( const Vector< Matrix< DDRMat > >& aParameters )
    {
        // store the parameters on the base class
        IQI::set_parameters( aParameters );

        // the parameter list is passed even if empty; the check that r and delta are set happens at evaluation
        uint tNumParameters = aParameters.size();
        mParametersSet      = tNumParameters >= 2;

        if ( !mParametersSet )
        {
            return;
        }

        MORIS_ERROR( tNumParameters <= 5,
                "IQI_Topological_Derivative::set_parameters - IQI %s: at most 5 parameters ( r / delta / l / chi / s ), %d given.",
                mName.c_str(),
                static_cast< int >( tNumParameters ) );

        // contrast ratio
        mContrast = aParameters( 0 )( 0 );
        MORIS_ERROR( mContrast >= 0.0,
                "IQI_Topological_Derivative::set_parameters - IQI %s: contrast r must be >= 0, got %e.",
                mName.c_str(),
                mContrast );

        // hole (-1) or inclusion (+1)
        mDelta = aParameters( 1 )( 0 );
        MORIS_ERROR( mDelta == -1.0 || mDelta == 1.0,
                "IQI_Topological_Derivative::set_parameters - IQI %s: delta must be -1 (material set) or +1 (void set), got %e.",
                mName.c_str(),
                mDelta );

        // area weight
        mAreaWeight = tNumParameters > 2 ? aParameters( 2 )( 0 ) : 0.0;

        // polarisation constant override (0: from the CM)
        mChiOverride = tNumParameters > 3 ? aParameters( 3 )( 0 ) : 0.0;
        MORIS_ERROR( mChiOverride >= 0.0,
                "IQI_Topological_Derivative::set_parameters - IQI %s: chi override must be 0 (from the CM) or positive, got %e.",
                mName.c_str(),
                mChiOverride );

        // output sign
        mOutputSign = tNumParameters > 4 ? aParameters( 4 )( 0 ) : 1.0;
        MORIS_ERROR( mOutputSign == -1.0 || mOutputSign == 1.0,
                "IQI_Topological_Derivative::set_parameters - IQI %s: output sign must be +1 or -1, got %e.",
                mName.c_str(),
                mOutputSign );
    }

    //------------------------------------------------------------------------------

    real
    IQI_Topological_Derivative::get_chi( Model_Type aModelType )
    {
        // override
        if ( mChiOverride > 0.0 )
        {
            return mChiOverride;
        }

        // Poisson ratio of the elasticity CM
        const std::shared_ptr< Property >& tPropPoisson =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) )->get_property( "PoissonRatio" );

        MORIS_ERROR( tPropPoisson != nullptr,
                "IQI_Topological_Derivative::get_chi - IQI %s: the elasticity CM has no PoissonRatio property; set the chi override.",
                mName.c_str() );

        const real tNu = tPropPoisson->val()( 0 );

        // chi = ( lambda + 3 mu ) / ( lambda + mu ) for the 2D model
        switch ( aModelType )
        {
            case Model_Type::PLANE_STRAIN:
                return 3.0 - 4.0 * tNu;
            case Model_Type::PLANE_STRESS:
                return ( 3.0 - tNu ) / ( 1.0 + tNu );
            default:
                MORIS_ERROR( false,
                        "IQI_Topological_Derivative::get_chi - IQI %s: chi from the CM needs PLANE_STRESS or PLANE_STRAIN.",
                        mName.c_str() );
                return 0.0;
        }
    }

    //------------------------------------------------------------------------------

    real
    IQI_Topological_Derivative::evaluate_topological_derivative()
    {
        MORIS_ERROR( mParametersSet,
                "IQI_Topological_Derivative - IQI %s: function_parameters must give at least r / delta.",
                mName.c_str() );

        // get the elasticity CM
        const std::shared_ptr< Constitutive_Model >& tCMElasticity =
                mLeaderCM( static_cast< uint >( IQI_Constitutive_Type::ELAST ) );

        MORIS_ERROR( tCMElasticity != nullptr,
                "IQI_Topological_Derivative - IQI %s: leader constitutive model \"Elast\" not set.",
                mName.c_str() );

        const Matrix< DDRMat >& tFlux   = tCMElasticity->flux();
        const Matrix< DDRMat >& tStrain = tCMElasticity->strain();

        // in-plane components of the Voigt vectors (engineering shear strain):
        // PLANE_STRESS [ xx, yy, xy ], PLANE_STRAIN [ xx, yy, zz, xy ]
        const Model_Type tModelType  = tCMElasticity->get_plane_type();
        uint             tShearIndex = 0;

        if ( tFlux.numel() == 3 && tModelType == Model_Type::PLANE_STRESS )
        {
            tShearIndex = 2;
        }
        else if ( tFlux.numel() == 4 && tModelType == Model_Type::PLANE_STRAIN )
        {
            tShearIndex = 3;
        }
        else
        {
            MORIS_ERROR( false,
                    "IQI_Topological_Derivative - IQI %s: 2D PLANE_STRESS or PLANE_STRAIN only (flux size %d); no 3D or axisymmetric formula.",
                    mName.c_str(),
                    static_cast< int >( tFlux.numel() ) );
        }

        MORIS_ASSERT( tStrain.numel() == tFlux.numel(),
                "IQI_Topological_Derivative - flux and strain sizes differ." );

        // sigma:eps, tr( sigma ), tr( eps ), in plane
        const real tFluxStrain  = tFlux( 0 ) * tStrain( 0 ) + tFlux( 1 ) * tStrain( 1 ) + tFlux( tShearIndex ) * tStrain( tShearIndex );
        const real tTraceFlux   = tFlux( 0 ) + tFlux( 1 );
        const real tTraceStrain = tStrain( 0 ) + tStrain( 1 );

        // polarisation constant
        const real tChi = this->get_chi( tModelType );
        const real tR   = mContrast;

        const real tDenominator1 = tChi * tR + 1.0;
        const real tDenominator2 = tChi + 2.0 * tR - 1.0;

        MORIS_ERROR( tDenominator1 > 0.0 && tDenominator2 > 0.0,
                "IQI_Topological_Derivative - IQI %s: chi r + 1 = %e and chi + 2r - 1 = %e must be positive (chi = %e, r = %e).",
                mName.c_str(),
                tDenominator1,
                tDenominator2,
                tChi,
                tR );

        // dJ = -(r-1)/(chi r+1) (chi+1)/2 [ 2 sigma:eps + (r-1)(chi-2)/(chi+2r-1) tr(sigma) tr(eps) ] + delta l
        const real tLeading = -( tR - 1.0 ) / tDenominator1 * ( tChi + 1.0 ) / 2.0;
        const real tTrace   = ( tR - 1.0 ) * ( tChi - 2.0 ) / tDenominator2;

        const real tDJ = tLeading * ( 2.0 * tFluxStrain + tTrace * tTraceFlux * tTraceStrain ) + mDelta * mAreaWeight;

        // generalised topological derivative g = s delta dJ
        return mOutputSign * mDelta * tDJ;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Topological_Derivative::compute_QI( Matrix< DDRMat >& aQI )
    {
        aQI.set_size( 1, 1, this->evaluate_topological_derivative() );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Topological_Derivative::compute_QI( real aWStar )
    {
        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // evaluate the QI
        mSet->get_QI()( tQIIndex ) += aWStar * this->evaluate_topological_derivative();
    }

    //------------------------------------------------------------------------------

    void
    IQI_Topological_Derivative::compute_dQIdu( real /*aWStar*/ )
    {
        MORIS_ERROR( false,
                "IQI_Topological_Derivative::compute_dQIdu - not implemented: the topological derivative is not differentiated." );
    }

    //------------------------------------------------------------------------------

    void
    IQI_Topological_Derivative::compute_dQIdu(
            Vector< MSI::Dof_Type >& /*aDofType*/,
            Matrix< DDRMat >& /*adQIdu*/ )
    {
        MORIS_ERROR( false,
                "IQI_Topological_Derivative::compute_dQIdu - not implemented: the topological derivative is not differentiated." );
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
