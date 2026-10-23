/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_CM_Struc_Linear_Isotropic.cpp
 *
 */

#include "cl_FEM_CM_Struc_Linear_Isotropic.hpp"
#include "fn_FEM_Geometric_Sensitivity.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"

#include "fn_trans.hpp"
#include "fn_norm.hpp"
#include "fn_eye.hpp"
#include "fn_dot.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    CM_Struc_Linear_Isotropic::CM_Struc_Linear_Isotropic()
    {
        // set the property pointer cell size
        mProperties.resize( mProperties.size() + static_cast< uint >( CM_Property_Type_Iso::MAX_ENUM ), nullptr );

        // populate the map
        mPropertyMap[ "YoungsModulus" ] = static_cast< uint >( CM_Property_Type_Iso::EMOD )    //
                                        + static_cast< uint >( CM_Property_Type_Lin::MAX_ENUM );

        mPropertyMap[ "PoissonRatio" ] = static_cast< uint >( CM_Property_Type_Iso::NU )    //
                                       + static_cast< uint >( CM_Property_Type_Lin::MAX_ENUM );
    }

    //------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::set_local_properties()
    {
        // set the parent properties
        CM_Struc_Linear::set_local_properties();

        // set the Young's modulus property
        mPropEMod = this->get_property( "YoungsModulus" );

        // set the Poisson ratio property
        mPropPoisson = this->get_property( "PoissonRatio" );

        // check that essential properties exist
        MORIS_ASSERT( mPropEMod,
                "CM_Struc_Linear_Isotropic::set_local_properties - Young's modulus property does not exist.\n" );

        MORIS_ASSERT( mPropPoisson,
                "CM_Struc_Linear_Isotropic::set_local_properties - Poisson ratio property does not exist.\n" );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::eval_const()
    {
        // get the Poisson's ratio value
        const real tNu = mPropPoisson->val()( 0 );

        // get the Young's modulus value
        const real tEmod = mPropEMod->val()( 0 );

        // evaluate the constitutive matrix
        ( this->*mConstFunc )( { tEmod, tNu } );
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    CM_Struc_Linear_Isotropic::eval_inv_bulk_modulus()
    {
        // get Poisson ratio value
        const real tNu = mPropPoisson->val()( 0 );

        // get elasticity modulus value
        const real tEMod = mPropEMod->val()( 0 );

        // init inverse of the bulk modulus
        real tInvBulkModulus;

        // evaluate inverse of the bulk modulus
        ( this->*m_eval_inv_bulk_modulus )( tNu, tEMod, tInvBulkModulus );

        // return
        return tInvBulkModulus;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::eval_inv_bulk_modulus_generic(
            const real& aNu,
            const real& aEMod,
            real&       aInvBulkModulus )
    {
        aInvBulkModulus = 3.0 * ( 1.0 - 2.0 * aNu ) / aEMod;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::eval_inv_bulk_modulus_plane_stress(
            const real& aNu,
            const real& aEMod,
            real&       aInvBulkModulus )
    {
        aInvBulkModulus = 2.0 * ( 1.0 - aNu ) / aEMod;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    CM_Struc_Linear_Isotropic::eval_dInvBulkModulusdDOF(
            const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // get the dof FI
        Field_Interpolator* tFI =
                mFIManager->get_field_interpolators_for_type( aDofTypes( 0 ) );

        // init inverse of the bulk modulus
        Matrix< DDRMat > tdInvBulkModulusdDOF( 1, tFI->get_number_of_space_time_coefficients(), 0.0 );

        // if Young's modulus property depends on dof type
        if ( mPropEMod->check_dof_dependency( aDofTypes ) )
        {
            tdInvBulkModulusdDOF -=
                    eval_inv_bulk_modulus() * mPropEMod->dPropdDOF( aDofTypes ) / mPropEMod->val()( 0 );
        }

        // if Poisson ratio property depends on dof type
        if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dInvBulkModulusdDOF - Poisson's ratio depends on dof, not handled." );
        }

        // return
        return tdInvBulkModulusdDOF;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::eval_dFluxdDOF( const Vector< MSI::Dof_Type >& aDofTypes )
    {
        // call the parent contribution
        CM_Struc_Linear::eval_dFluxdDOF( aDofTypes );

        // get the dof type as a uint
        const uint tDofType = static_cast< uint >( aDofTypes( 0 ) );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( tDofType );

        // if elastic modulus depends on dof type
        if ( mPropEMod->check_dof_dependency( aDofTypes ) )
        {
            // compute derivative with indirect dependency through properties
            mdFluxdDof( tDofIndex ) +=
                    this->constitutive() * this->strain() * mPropEMod->dPropdDOF( aDofTypes ) / mPropEMod->val()( 0 );
        }

        // if Poisson ratio depends on dof type
        if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dFluxdDOF - Poisson's ratio depends on dof, not handled." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::eval_dTestTractiondDOF(
            const Vector< MSI::Dof_Type >& aDofTypes,
            const Matrix< DDRMat >&        aNormal,
            const Matrix< DDRMat >&        aJump,
            const Vector< MSI::Dof_Type >& aTestDofTypes )
    {
        // compute generic derivative of test traction
        CM_Struc_Linear::eval_dTestTractiondDOF( aDofTypes, aNormal, aJump, aTestDofTypes );

        // get test dof type index
        const uint tTestDofIndex = mDofTypeMap( static_cast< uint >( aTestDofTypes( 0 ) ) );

        // get the dof type index
        const uint tDofIndex = mGlobalDofTypeMap( static_cast< uint >( aDofTypes( 0 ) ) );

        // if elastic modulus depends on dof type
        if ( mPropEMod->check_dof_dependency( aDofTypes ) )
        {
            // compute derivative
            mdTestTractiondDof( tTestDofIndex )( tDofIndex ) =
                    trans( this->testTraction( aNormal, aTestDofTypes ) ) * aJump * mPropEMod->dPropdDOF( aDofTypes ) / mPropEMod->val()( 0 );
        }

        // if Poisson's ratio depends on dof type
        if ( mPropPoisson->check_dof_dependency( aDofTypes ) )
        {
            MORIS_ERROR( false, "CM_Struc_Linear_Isotropic::eval_dFluxdDOF - Poisson's ratio depends on dof, not handled." );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    CM_Struc_Linear_Isotropic::get_e_prime()
    {
        real tEPrime;

        // get elasticity modulus value
        const real tE = mPropEMod->val()( 0 );

        // get Poisson ratio value
        const real tNu = mPropPoisson->val()( 0 );

        switch ( mPlaneType )
        {
            case Model_Type::PLANE_STRESS:
            {
                // Eprime = E
                tEPrime = tE;
                break;
            }
            case Model_Type::PLANE_STRAIN:
            {
                // Eprime = E / ( 1 - nu^2 )
                tEPrime = tE / ( 1 - ( tNu * tNu ) );
                break;
            }
            default:
            {
                tEPrime = tE;
                MORIS_ASSERT( false, "CM_Struc_Linear_Isotropic::get_e_prime() - unknown model type " );
                break;
            }
        }
        return tEPrime;
    }

    //--------------------------------------------------------------------------------------------------------------

    bool
    CM_Struc_Linear_Isotropic::has_spatial_constitutive_derivatives()
    {
        // the Young's modulus and the Poisson ratio may vary in space with a known derivative
        if ( !has_known_spatial_variation( mPropEMod ) || !has_known_spatial_variation( mPropPoisson ) )
        {
            return false;
        }

        // the remaining properties (thermal, eigen strain, rotation axis) must be constant
        for ( const std::shared_ptr< Property >& tProperty : mProperties )
        {
            if ( tProperty != mPropEMod && tProperty != mPropPoisson && !is_spatially_constant( tProperty ) )
            {
                return false;
            }
        }
        return true;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::eval_dConstdx(
            uint              aDirection,
            Matrix< DDRMat >& adConstdx )
    {
        // constitutive matrix
        const Matrix< DDRMat >& tConst = this->constitutive();

        adConstdx.set_size( tConst.n_rows(), tConst.n_cols(), 0.0 );

        // spatial derivatives of the Young's modulus and the Poisson ratio
        Matrix< DDRMat > tdEdx;
        eval_property_space_derivative( mPropEMod, aDirection, tdEdx );

        Matrix< DDRMat > tdNudx;
        eval_property_space_derivative( mPropPoisson, aDirection, tdNudx );

        // C is linear in E
        if ( tdEdx( 0 ) != 0.0 )
        {
            adConstdx += tConst * ( tdEdx( 0 ) / mPropEMod->val()( 0 ) );
        }

        // Poisson ratio contribution
        if ( tdNudx( 0 ) != 0.0 )
        {
            Matrix< DDRMat > tdConstdNu;
            this->eval_dConstdNu( tdConstdNu );

            adConstdx += tdConstdNu * tdNudx( 0 );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::eval_dConstdNu( Matrix< DDRMat >& adConstdNu )
    {
        MORIS_ERROR( mPlaneType != Model_Type::AXISYMMETRIC,
                "CM_Struc_Linear_Isotropic::eval_dConstdNu - axisymmetric model not supported." );

        // constitutive matrix, Young's modulus and Poisson ratio
        const Matrix< DDRMat >& tConst = this->constitutive();

        const real tEmod = mPropEMod->val()( 0 );
        const real tNu   = mPropPoisson->val()( 0 );

        // deviatoric tensors: C = E / ( c ( 1 + nu ) ) * const, so dC/dnu = - C / ( 1 + nu )
        if ( mTensorType == Model_Type::DEVIATORIC )
        {
            adConstdNu = tConst * ( -1.0 / ( 1.0 + tNu ) );
            return;
        }

        // full tensors: C = E ( a(nu) on the normal diagonal, b(nu) off the normal diagonal, s(nu) on the shear diagonal )
        real tdA;
        real tdB;
        if ( mSpaceDim == 2 && mPlaneType == Model_Type::PLANE_STRESS )
        {
            // a = 1 / ( 1 - nu^2 ), b = nu / ( 1 - nu^2 )
            const real tDen = ( 1.0 - tNu * tNu ) * ( 1.0 - tNu * tNu );

            tdA = 2.0 * tNu / tDen;
            tdB = ( 1.0 + tNu * tNu ) / tDen;
        }
        else
        {
            // plane strain and 3D: a = ( 1 - nu ) / D, b = nu / D, D = ( 1 + nu ) ( 1 - 2 nu )
            const real tD   = ( 1.0 + tNu ) * ( 1.0 - 2.0 * tNu );
            const real tDen = tD * tD;

            tdA = 2.0 * tNu * ( 2.0 - tNu ) / tDen;
            tdB = ( 1.0 + 2.0 * tNu * tNu ) / tDen;
        }

        // s = 1 / ( 2 ( 1 + nu ) ) for all full tensors
        const real tdS = -0.5 / ( ( 1.0 + tNu ) * ( 1.0 + tNu ) );

        adConstdNu.set_size( tConst.n_rows(), tConst.n_cols(), 0.0 );

        for ( uint iRow = 0; iRow < mNumNormalStress; iRow++ )
        {
            for ( uint iCol = 0; iCol < mNumNormalStress; iCol++ )
            {
                adConstdNu( iRow, iCol ) = tEmod * ( iRow == iCol ? tdA : tdB );
            }
        }

        for ( uint iRow = mNumNormalStress; iRow < tConst.n_rows(); iRow++ )
        {
            adConstdNu( iRow, iRow ) = tEmod * tdS;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::full_plane_stress( std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        const real tPre = tEmod / ( 1 - std::pow( tNu, 2 ) );

        mConst( 0, 0 ) = tPre;
        mConst( 1, 1 ) = tPre;
        mConst( 0, 1 ) = tPre * tNu;
        mConst( 1, 0 ) = tPre * tNu;
        mConst( 2, 2 ) = tPre * 0.5 * ( 1.0 - tNu );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::deviatoric_plane_stress( std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        const real tPre = tEmod / ( ( 1 + tNu ) * 2.0 );

        mConst( 0, 0 ) = tPre;
        mConst( 1, 1 ) = tPre;
        mConst( 0, 1 ) = -tPre;
        mConst( 1, 0 ) = -tPre;
        mConst( 2, 2 ) = tPre;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::full_plane_strain( std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        const real tPre = tEmod / ( 1.0 + tNu ) / ( 1.0 - 2.0 * tNu );

        mConst( 0, 0 ) = tPre * ( 1.0 - tNu );
        mConst( 0, 1 ) = tPre * tNu;
        mConst( 0, 2 ) = tPre * tNu;
        mConst( 1, 0 ) = tPre * tNu;
        mConst( 1, 1 ) = tPre * ( 1.0 - tNu );
        mConst( 1, 2 ) = tPre * tNu;
        mConst( 2, 0 ) = tPre * tNu;
        mConst( 2, 1 ) = tPre * tNu;
        mConst( 2, 2 ) = tPre * ( 1.0 - tNu );
        mConst( 3, 3 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::deviatoric_plane_strain( std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        const real tPre = tEmod / ( 3.0 * ( 1.0 + tNu ) );

        mConst( 0, 0 ) = tPre * 4.0;
        mConst( 0, 1 ) = tPre;
        mConst( 0, 2 ) = tPre;
        mConst( 1, 0 ) = tPre;
        mConst( 1, 1 ) = tPre * 4.0;
        mConst( 1, 2 ) = tPre;
        mConst( 2, 0 ) = tPre;
        mConst( 2, 1 ) = tPre;
        mConst( 2, 2 ) = tPre * 4.0;
        mConst( 3, 3 ) = tPre * 3.0 / 2.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::full_axisymmetric( std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        const real tPre = tEmod / ( 1.0 + tNu ) / ( 1.0 - 2.0 * tNu );

        mConst( 0, 0 ) = tPre * ( 1.0 - tNu );
        mConst( 0, 1 ) = tPre * tNu;
        mConst( 1, 0 ) = tPre * tNu;
        mConst( 1, 1 ) = tPre * ( 1.0 - tNu );
        mConst( 3, 3 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;

        // axisymmetric contribution
        mConst( 0, 2 ) = tPre * tNu;
        mConst( 1, 2 ) = tPre * tNu;
        mConst( 2, 0 ) = tPre * tNu;
        mConst( 2, 1 ) = tPre * tNu;
        mConst( 2, 2 ) = tPre * ( 1.0 - tNu );
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::deviatoric_axisymmetric( std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        real tPre = tEmod / ( 3.0 * ( 1.0 + tNu ) );

        mConst( 0, 0 ) = tPre * 4.0;
        mConst( 0, 1 ) = tPre;
        mConst( 1, 0 ) = tPre;
        mConst( 1, 1 ) = tPre * 4.0;
        mConst( 3, 3 ) = tPre * 3.0 / 2.0;

        // axisymmetric contribution
        mConst( 0, 2 ) = tPre;
        mConst( 1, 2 ) = tPre;
        mConst( 2, 0 ) = tPre;
        mConst( 2, 1 ) = tPre;
        mConst( 2, 2 ) = tPre * 4.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::full_3d( std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        const real tPre = tEmod / ( 1.0 + tNu ) / ( 1.0 - 2.0 * tNu );

        mConst( 0, 0 ) = tPre * ( 1.0 - tNu );
        mConst( 0, 1 ) = tPre * tNu;
        mConst( 0, 2 ) = tPre * tNu;
        mConst( 1, 0 ) = tPre * tNu;
        mConst( 1, 1 ) = tPre * ( 1.0 - tNu );
        mConst( 1, 2 ) = tPre * tNu;
        mConst( 2, 0 ) = tPre * tNu;
        mConst( 2, 1 ) = tPre * tNu;
        mConst( 2, 2 ) = tPre * ( 1.0 - tNu );
        mConst( 3, 3 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
        mConst( 4, 4 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
        mConst( 5, 5 ) = tPre * ( 1.0 - 2.0 * tNu ) / 2.0;
    }

    //--------------------------------------------------------------------------------------------------------------

    void
    CM_Struc_Linear_Isotropic::deviatoric_3d(
            std::initializer_list< const real >&& tParams )
    {
        const real& tEmod = tParams.begin()[ 0 ];
        const real& tNu   = tParams.begin()[ 1 ];

        const real tPre = tEmod / ( 3.0 * ( 1.0 + tNu ) );

        mConst( 0, 0 ) = tPre * 4.0;
        mConst( 0, 1 ) = tPre;
        mConst( 0, 2 ) = tPre;
        mConst( 1, 0 ) = tPre;
        mConst( 1, 1 ) = tPre * 4.0;
        mConst( 1, 2 ) = tPre;
        mConst( 2, 0 ) = tPre;
        mConst( 2, 1 ) = tPre;
        mConst( 2, 2 ) = tPre * 4.0;
        mConst( 3, 3 ) = tPre * 3.0 / 2.0;
        mConst( 4, 4 ) = tPre * 3.0 / 2.0;
        mConst( 5, 5 ) = tPre * 3.0 / 2.0;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::fem
