/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Summed_Fields.cpp
 *
 */

#include "cl_GEN_Summed_Fields.hpp"

namespace moris::gen
{

    //--------------------------------------------------------------------------------------------------------------

    Summed_Fields::Summed_Fields(
            const Vector< std::shared_ptr< Field > >& aFields,
            std::string                               aName )
            : Field( Vector< ADV >(), std::move( aName ) )
            , mFields( aFields )
    {
        MORIS_ERROR( mFields.size() > 0, "A GEN Summed_Fields must be created with at least one field." );
    }

    //--------------------------------------------------------------------------------------------------------------

    void Summed_Fields::import_advs( sol::Dist_Vector* aOwnedADVs )
    {
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            iField->import_advs( aOwnedADVs );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Summed_Fields::reset_nodal_data( mtk::Interpolation_Mesh* aMesh )
    {
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            iField->reset_nodal_data( aMesh );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    bool Summed_Fields::has_advs()
    {
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            if ( iField->has_advs() )
            {
                return true;
            }
        }
        return false;
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Summed_Fields::get_coefficient_adv_ids()
    {
        Vector< sint > tResult;
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            Vector< sint > tChildIDs = iField->get_coefficient_adv_ids();
            for ( sint iID : tChildIDs )
            {
                tResult.push_back( iID );
            }
        }
        return tResult;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Summed_Fields::update_dependencies( const Vector< std::shared_ptr< Field > >& aUpdatedFields )
    {
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            iField->update_dependencies( aUpdatedFields );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    real Summed_Fields::get_field_value(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        real tResult = 0.0;
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            tResult += iField->get_field_value( aNodeIndex, aCoordinates );
        }
        return tResult;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Summed_Fields::get_field_value(
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        real tResult = 0.0;
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            tResult += iField->get_field_value( aDerivedNode, aNodeManager );
        }
        return tResult;
    }

    //--------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >& Summed_Fields::get_dfield_dadvs(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        // Copy each child's result before calling the next child, since a child may return a reference to its own buffer
        Vector< Matrix< DDRMat > > tChildSensitivities( mFields.size() );
        uint                       tTotal = 0;
        for ( uint iFieldIndex = 0; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            tChildSensitivities( iFieldIndex ) = mFields( iFieldIndex )->get_dfield_dadvs( aNodeIndex, aCoordinates );
            tTotal += tChildSensitivities( iFieldIndex ).numel();
        }

        // Concatenate in child order
        mSensitivities.set_size( 1, tTotal );
        uint tOffset = 0;
        for ( const Matrix< DDRMat >& iChild : tChildSensitivities )
        {
            for ( uint iEntry = 0; iEntry < iChild.numel(); iEntry++ )
            {
                mSensitivities( tOffset + iEntry ) = iChild( iEntry );
            }
            tOffset += iChild.numel();
        }
        return mSensitivities;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Summed_Fields::get_dfield_dadvs(
            Matrix< DDRMat >&   aSensitivities,
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        Vector< Matrix< DDRMat > > tChildSensitivities( mFields.size() );
        uint                       tTotal = 0;
        for ( uint iFieldIndex = 0; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            mFields( iFieldIndex )->get_dfield_dadvs( tChildSensitivities( iFieldIndex ), aDerivedNode, aNodeManager );
            tTotal += tChildSensitivities( iFieldIndex ).numel();
        }

        // Concatenate in child order
        aSensitivities.set_size( 1, tTotal );
        uint tOffset = 0;
        for ( const Matrix< DDRMat >& iChild : tChildSensitivities )
        {
            for ( uint iEntry = 0; iEntry < iChild.numel(); iEntry++ )
            {
                aSensitivities( tOffset + iEntry ) = iChild( iEntry );
            }
            tOffset += iChild.numel();
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Summed_Fields::get_determining_adv_ids(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates )
    {
        Vector< sint > tResult;
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            Vector< sint > tChildIDs = iField->get_determining_adv_ids( aNodeIndex, aCoordinates );
            for ( sint iID : tChildIDs )
            {
                tResult.push_back( iID );
            }
        }
        return tResult;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Summed_Fields::get_determining_adv_ids(
            Vector< sint >&     aDeterminingADVIDs,
            const Derived_Node& aDerivedNode,
            const Node_Manager& aNodeManager )
    {
        aDeterminingADVIDs.clear();
        for ( const std::shared_ptr< Field >& iField : mFields )
        {
            Vector< sint > tChildIDs;
            iField->get_determining_adv_ids( tChildIDs, aDerivedNode, aNodeManager );
            for ( sint iID : tChildIDs )
            {
                aDeterminingADVIDs.push_back( iID );
            }
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    void Summed_Fields::get_dfield_dcoordinates(
            uint                    aNodeIndex,
            const Matrix< DDRMat >& aCoordinates,
            Matrix< DDRMat >&       aSensitivities )
    {
        mFields( 0 )->get_dfield_dcoordinates( aNodeIndex, aCoordinates, aSensitivities );
        for ( uint iFieldIndex = 1; iFieldIndex < mFields.size(); iFieldIndex++ )
        {
            // Children write into a pre-sized buffer, so give each the same size as the caller's
            Matrix< DDRMat > tChild( aSensitivities.n_rows(), aSensitivities.n_cols(), 0.0 );
            mFields( iFieldIndex )->get_dfield_dcoordinates( aNodeIndex, aCoordinates, tChild );
            aSensitivities = aSensitivities + tChild;
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    std::shared_ptr< mtk::Field > Summed_Fields::get_mtk_field()
    {
        // Remeshing (remapping onto a new mesh via an MTK field) is unsupported for summed designs
        return nullptr;
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
