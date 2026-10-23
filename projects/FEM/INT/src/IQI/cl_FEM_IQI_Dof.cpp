/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_FEM_IQI_Dof.cpp
 *
 */

#include "cl_FEM_Set.hpp"
#include "cl_FEM_Field_Interpolator_Manager.hpp"
#include "cl_FEM_IQI_Dof.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------

    IQI_Dof::IQI_Dof() {}

    //------------------------------------------------------------------------------

    void
    IQI_Dof::initialize()
    {
        if ( !mIsInitialized )
        {
            // size of parameter list
            uint tParamSize = mParameters.size();

            // extract spatial derivative order
            if ( tParamSize > 0 )
            {
                MORIS_ERROR( mParameters( 0 ).numel() == 2,
                        "IQI_Dof::initialize - Spatial gradient definition requires exactly two coefficients.\n" );

                mSpatialDerivativeDirection = mParameters( 0 )( 0 );
                mSpatialDerivativeOrder     = mParameters( 0 )( 1 );
            }

            // extract time derivative order
            if ( tParamSize > 1 )
            {
                MORIS_ERROR( mSpatialDerivativeOrder == 0,
                        "IQI_Dof::initialize - Time gradient can only be computed if spatial gradient order is zero.\n" );

                MORIS_ERROR( mParameters( 1 ).numel() == 1,
                        "IQI_Dof::initialize - Time gradient definition requires exactly one coefficient.\n" );

                mTimeDerivativeOrder = mParameters( 1 )( 0 );
            }

            // set initialize flag to true
            mIsInitialized = true;
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof::compute_QI( real aWStar )
    {
        // initialize if needed
        this->initialize();

        // get index for QI
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        Matrix< DDRMat > tMat;

        this->evaluate_QI( tMat );

        // a vector field needs a component index for a scalar QI
        MORIS_ERROR( tMat.numel() == 1,
                "IQI_Dof::compute_QI - IQI %s: vectorial_field_index not set for a vector field.",
                mName.c_str() );

        mSet->get_QI()( tQIIndex ) += aWStar * tMat;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof::compute_QI( Matrix< DDRMat >& aQI )
    {
        // initialize if needed
        this->initialize();

        // evaluate QI
        this->evaluate_QI( aQI );
    }

    //------------------------------------------------------------------------------

    sint
    IQI_Dof::get_dof_component_index( Field_Interpolator* aFI ) const
    {
        // number of field components
        uint tNumComps = aFI->get_number_of_fields();

        // a set index is used as given
        if ( mIQITypeIndex != -1 )
        {
            return this->get_component_index( tNumComps );
        }

        // scalar field
        if ( tNumComps == 1 )
        {
            return 0;
        }

        // single dof type of a vector group: its position in the group
        if ( mQuantityDofType.size() == 1 )
        {
            const Vector< MSI::Dof_Type >& tGroup = aFI->get_dof_type();

            for ( uint iComp = 0; iComp < tGroup.size(); iComp++ )
            {
                if ( tGroup( iComp ) == mQuantityDofType( 0 ) )
                {
                    return iComp;
                }
            }
        }

        // vector field without an index
        return -1;
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof::evaluate_QI( Matrix< DDRMat >& aMat )
    {
        // get field interpolator for a given dof type
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // check that field interpolator exists
        MORIS_ASSERT( tFI != nullptr,
                "IQI_Dof::evaluate_QI - field interpolator does not exist." );

        // get the field component
        sint tIndex = this->get_dof_component_index( tFI );

        // evaluate spatial derivative of dof
        if ( mSpatialDerivativeOrder > 0 )
        {
            MORIS_ERROR( tIndex != -1,
                    "IQI_Dof::evaluate_QI - IQI %s: vectorial_field_index not set for a vector field.",
                    mName.c_str() );

            const Matrix< DDRMat >& tSpatialGradient = tFI->gradx( mSpatialDerivativeOrder );

            aMat = { tSpatialGradient( mSpatialDerivativeDirection, tIndex ) };
        }

        // evaluate time derivative of dof
        else if ( mTimeDerivativeOrder > 0 )
        {
            MORIS_ERROR( tIndex != -1,
                    "IQI_Dof::evaluate_QI - IQI %s: vectorial_field_index not set for a vector field.",
                    mName.c_str() );

            const Matrix< DDRMat >& tTemporalGradient = tFI->gradt( mTimeDerivativeOrder );

            aMat = { tTemporalGradient( tIndex ) };
        }
        else if ( tIndex != -1 )
        {
            // evaluate DOF value
            aMat = { tFI->val()( tIndex ) };
        }
        // DO NOT DELETE THIS FUNCTIONALITY AGAIN
        else
        {
            aMat = tFI->val();
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof::compute_dQIdu( real aWStar )
    {
        // get field interpolator for a given dof type
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // get the column index to assemble in residual
        sint tQIIndex = mSet->get_QI_assembly_index( mName );

        // get the field component; derivatives need one
        sint tIndex = this->get_dof_component_index( tFI );

        MORIS_ERROR( tIndex != -1,
                "IQI_Dof::compute_dQIdu - IQI %s: vectorial_field_index not set for a vector field.",
                mName.c_str() );

        // get the number of leader dof type dependencies
        uint tNumDofDependencies = mRequestedLeaderGlobalDofTypes.size();

        // compute dQIdu for indirect dof dependencies
        for ( uint iDof = 0; iDof < tNumDofDependencies; iDof++ )
        {
            // get the treated dof type
            Vector< MSI::Dof_Type >& tDofType = mRequestedLeaderGlobalDofTypes( iDof );

            // get leader index for residual dof type, indices for assembly
            uint tLeaderDofIndex      = mSet->get_dof_index_for_type( tDofType( 0 ), mtk::Leader_Follower::LEADER );
            uint tLeaderDepStartIndex = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 0 );
            uint tLeaderDepStopIndex  = mSet->get_res_dof_assembly_map()( tLeaderDofIndex )( 0, 1 );

            // if derivative dof type is the group of the quantity dof type
            if ( tDofType( 0 ) == tFI->get_dof_type()( 0 ) )
            {
                // evaluate spatial derivative of dof
                if ( mSpatialDerivativeOrder > 0 )
                {
                    // Fixme: the following should be provided directly by field interpolator
                    // get number of field and number of bases
                    uint tNumVecFieldComps = tFI->val().numel();
                    uint tNumBasis         = tFI->dnNdxn( mSpatialDerivativeOrder ).n_cols();

                    Matrix< DDRMat > tSpatialDerivativeShapeFunction( 1, tNumVecFieldComps * tNumBasis, 0.0 );

                    const Matrix< DDRMat >& tdSpatialGradientdu = tFI->dnNdxn( mSpatialDerivativeOrder );

                    tSpatialDerivativeShapeFunction( { 0, 0 },                            //
                            { tIndex * tNumBasis, ( tIndex + 1 ) * tNumBasis - 1 } ) =    //
                            tdSpatialGradientdu.get_row( mSpatialDerivativeDirection );

                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex }, { 0, 0 } ) +=    //
                            aWStar * trans( tSpatialDerivativeShapeFunction );
                }

                // evaluate time derivative of dof
                else if ( mTimeDerivativeOrder > 0 )
                {
                    // Fixme: the following should be provided directly by field interpolator
                    // get number of field and number of bases
                    uint tNumVecFieldComps = tFI->val().numel();
                    uint tNumBasis         = tFI->dnNdtn( mTimeDerivativeOrder ).n_cols();

                    Matrix< DDRMat > tTimeDerivativeShapeFunction( 1, tNumVecFieldComps * tNumBasis, 0.0 );

                    tTimeDerivativeShapeFunction( { 0, 0 },                               //
                            { tIndex * tNumBasis, ( tIndex + 1 ) * tNumBasis - 1 } ) =    //
                            tFI->dnNdtn( mTimeDerivativeOrder ).matrix_data();

                    // get dof derivative of time derivative
                    // assemble into residual vector
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex }, { 0, 0 } ) +=    //
                            aWStar * trans( tTimeDerivativeShapeFunction );
                }
                else
                {
                    // build selection matrix
                    uint tNumVecFieldComps = tFI->val().numel();

                    Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );

                    tSelect( tIndex, 0 ) = 1.0;

                    // assemble into residual vector
                    mSet->get_residual()( tQIIndex )(
                            { tLeaderDepStartIndex, tLeaderDepStopIndex }, { 0, 0 } ) +=
                            aWStar * tFI->N_trans() * tSelect;
                }
            }
        }
    }

    //------------------------------------------------------------------------------

    void
    IQI_Dof::compute_dQIdu(
            Vector< MSI::Dof_Type >& aDofType,
            Matrix< DDRMat >&        adQIdu )
    {
        // get field interpolator for a given dof type
        Field_Interpolator* tFI =
                mLeaderFIManager->get_field_interpolators_for_type( mQuantityDofType( 0 ) );

        // get the field component; derivatives need one
        sint tIndex = this->get_dof_component_index( tFI );

        MORIS_ERROR( tIndex != -1,
                "IQI_Dof::compute_dQIdu - IQI %s: vectorial_field_index not set for a vector field.",
                mName.c_str() );

        // if derivative dof type is the group of the quantity dof type
        if ( aDofType( 0 ) == tFI->get_dof_type()( 0 ) )
        {
            // evaluate spatial derivative of dof
            if ( mSpatialDerivativeOrder > 0 )
            {
                // Fixme: the following should be provided directly by field interpolator
                // get number of field and number of bases
                uint tNumVecFieldComps = tFI->val().numel();
                uint tNumBasis         = tFI->dnNdxn( mSpatialDerivativeOrder ).n_cols();

                Matrix< DDRMat > tSpatialDerivativeShapeFunction( 1, tNumVecFieldComps * tNumBasis, 0.0 );

                const Matrix< DDRMat >& tdSpatialGradientdu = tFI->dnNdxn( mSpatialDerivativeOrder );

                tSpatialDerivativeShapeFunction( { 0, 0 },                            //
                        { tIndex * tNumBasis, ( tIndex + 1 ) * tNumBasis - 1 } ) =    //
                        tdSpatialGradientdu.get_row( mSpatialDerivativeDirection );

                // assemble into dof derivative of IQI
                adQIdu = tSpatialDerivativeShapeFunction;
            }

            // evaluate time derivative of dof
            else if ( mTimeDerivativeOrder > 0 )
            {
                // Fixme: the following should be provided directly by field interpolator
                // get number of field and number of bases
                uint tNumVecFieldComps = tFI->val().numel();
                uint tNumBasis         = tFI->dnNdtn( mTimeDerivativeOrder ).n_cols();

                Matrix< DDRMat > tTimeDerivativeShapeFunction( 1, tNumVecFieldComps * tNumBasis, 0.0 );

                tTimeDerivativeShapeFunction( { 0, 0 },                               //
                        { tIndex * tNumBasis, ( tIndex + 1 ) * tNumBasis - 1 } ) =    //
                        tFI->dnNdtn( mTimeDerivativeOrder ).matrix_data();

                // assemble into dof derivative of IQI
                adQIdu = tTimeDerivativeShapeFunction;
            }
            else
            {
                // build selection matrix
                uint tNumVecFieldComps = tFI->val().numel();

                Matrix< DDRMat > tSelect( tNumVecFieldComps, 1, 0.0 );

                tSelect( tIndex, 0 ) = 1.0;

                // assemble into dof derivative of IQI
                adQIdu = tFI->N_trans() * tSelect;
            }
        }
    }

    //------------------------------------------------------------------------------
}    // namespace moris::fem
