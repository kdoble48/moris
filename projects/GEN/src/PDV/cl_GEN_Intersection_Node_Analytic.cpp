/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Analytic.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Analytic.hpp"
#include "cl_GEN_Basis_Node.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "fn_eye.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_dot.hpp"
#include "fn_stringify_matrix.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Analytic::Intersection_Node_Analytic(
            uint                              aNodeIndex,
            const Vector< Background_Node* >& aBackgroundNodes,
            const Parent_Node&                aFirstParentNode,
            const Parent_Node&                aSecondParentNode,
            mtk::Geometry_Type                aBackgroundGeometryType,
            mtk::Interpolation_Order          aBackgroundInterpolationOrder,
            Level_Set_Geometry&               aInterfaceGeometry )
            : Intersection_Node_Level_Set(
                    aNodeIndex,
                    aBackgroundNodes,
                    aFirstParentNode,
                    aSecondParentNode,
                    aBackgroundGeometryType,
                    aBackgroundInterpolationOrder,
                    aInterfaceGeometry )
            , mLinearFallback( aInterfaceGeometry.last_analytic_intersection_fell_back() )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Intersection_Node_Analytic::get_field_basis_nodes() const
    {
        return this->get_locator_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    real Intersection_Node_Analytic::get_dresidual_dxi( Matrix< DDRMat >& aGradient ) const
    {
        const Basis_Node&       tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node&       tSecondParentNode = this->get_second_parent_node();
        const Matrix< DDRMat >& tFirstX           = tFirstParentNode.get_global_coordinates();
        const Matrix< DDRMat >& tSecondX          = tSecondParentNode.get_global_coordinates();
        Matrix< DDRMat >        tEdge             = tSecondX - tFirstX;
        real                    tStep             = 1e-6 * norm( tEdge );

        aGradient = mInterfaceGeometry.get_dfield_dcoordinates_at_point( this->get_global_coordinates(), tStep );

        real tDrDxi = 0.5 * dot( aGradient, tEdge );

        // dr/dxi = 0 (the field touches the isocontour tangentially along the edge, or is flat there): the root has no
        // derivative with respect to the ADVs or the parent coordinates, so no sensitivity is returned
        MORIS_ERROR( std::isfinite( tDrDxi ) and std::abs( tDrDxi ) > MORIS_REAL_MIN,
                "Intersection_Node_Analytic: d phi / d xi = %e at the analytic intersection %s on the edge %s - %s; the field is "
                "flat or tangent to the edge there and the intersection has no coordinate sensitivity.",
                tDrDxi,
                ios::stringify_log( this->get_global_coordinates() ).c_str(),
                ios::stringify_log( tFirstX ).c_str(),
                ios::stringify_log( tSecondX ).c_str() );
        return tDrDxi;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Intersection_Node_Analytic::get_parent_field_gradient( const Basis_Node& aParentNode ) const
    {
        // exact spatial gradient of the field at the parent (central difference, step 1e-6 of the edge length); the
        // parent value enters the linear crossing exactly, so the chord slope of the parent's own locators would not do
        // (MB-083)
        Matrix< DDRMat > tEdge     = this->get_second_parent_node().get_global_coordinates() - this->get_first_parent_node().get_global_coordinates();
        Matrix< DDRMat > tGradient = mInterfaceGeometry.get_dfield_dcoordinates_at_point( aParentNode.get_global_coordinates(), 1e-6 * norm( tEdge ) );
        for ( uint iDim = 0; iDim < tGradient.numel(); iDim++ )
        {
            MORIS_ERROR( std::isfinite( tGradient( iDim ) ),
                    "Intersection_Node_Analytic: the field gradient at the parent %s of a linear-fallback intersection is not finite.",
                    ios::stringify_log( aParentNode.get_global_coordinates() ).c_str() );
        }
        return tGradient;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Intersection_Node_Analytic::get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const
    {
        // d xi / d phi_i of the linear crossing xi = ( 2 c - phi_1 - phi_2 ) / ( phi_2 - phi_1 ), as Intersection_Node_Linear;
        // used for the linear-fallback sensitivities (the analytic ones take the field sensitivities at x*)
        real tThreshold = mInterfaceGeometry.get_isocontour_threshold();
        real tPhi0      = mInterfaceGeometry.get_field_value( this->get_first_parent_node().get_index(), this->get_first_parent_node().get_global_coordinates() );
        real tPhi1      = mInterfaceGeometry.get_field_value( this->get_second_parent_node().get_index(), this->get_second_parent_node().get_global_coordinates() );
        return 2 * ( ( tPhi0 - tThreshold ) * ( aAncestorIndex == 1 ) - ( tPhi1 - tThreshold ) * ( aAncestorIndex == 0 ) ) / std::pow( ( tPhi1 - tPhi0 ), 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Intersection_Node_Analytic::get_dxi_dcoordinate_first_parent() const
    {
        if ( mLinearFallback )
        {
            return this->get_dxi_dfield_from_ancestor( 0 ) * this->get_parent_field_gradient( this->get_first_parent_node() );
        }
        Matrix< DDRMat > tGradient;
        real             tDrDxi = this->get_dresidual_dxi( tGradient );
        return ( -0.5 * ( 1.0 - this->get_local_coordinate() ) / tDrDxi ) * tGradient;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Intersection_Node_Analytic::get_dxi_dcoordinate_second_parent() const
    {
        if ( mLinearFallback )
        {
            return this->get_dxi_dfield_from_ancestor( 1 ) * this->get_parent_field_gradient( this->get_second_parent_node() );
        }
        Matrix< DDRMat > tGradient;
        real             tDrDxi = this->get_dresidual_dxi( tGradient );
        return ( -0.5 * ( 1.0 + this->get_local_coordinate() ) / tDrDxi ) * tGradient;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Intersection_Node_Analytic::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        // linear fallback: the coordinate is the linear crossing of the exact end values; its sensitivities are the
        // linear node's, through the field sensitivities at the parents and the exact field gradient at the parents
        if ( mLinearFallback )
        {
            Intersection_Node_Level_Set::append_dcoordinate_dadv( aCoordinateSensitivities, aSensitivityFactor );
            return;
        }

        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();
        real              tXi               = this->get_local_coordinate();

        // parent vector (column)
        Matrix< DDRMat > tParentVector = trans( tSecondParentNode.get_global_coordinates() - tFirstParentNode.get_global_coordinates() );

        // spatial gradient and d r / d xi at the intersection point, once for all terms below
        Matrix< DDRMat > tGradient;
        real             tDrDxi = this->get_dresidual_dxi( tGradient );

        // d x* / d p = 0.5 ( x_2 - x_1 ) d xi* / d p,  d xi* / d p = - ( d phi / d p )( x* ) / ( d r / d xi )
        Matrix< DDRMat > tFieldSensitivities = mInterfaceGeometry.get_dfield_dadvs_at_point( this->get_global_coordinates() );
        Matrix< DDRMat > tSensitivitiesToAdd = ( -0.5 / tDrDxi ) * aSensitivityFactor * tParentVector * tFieldSensitivities;

        // join
        uint tJoinedSensitivityLength = aCoordinateSensitivities.n_cols();
        aCoordinateSensitivities.resize( tSensitivitiesToAdd.n_rows(), tJoinedSensitivityLength + tSensitivitiesToAdd.n_cols() );
        for ( uint iCoordinate = 0; iCoordinate < tSensitivitiesToAdd.n_rows(); iCoordinate++ )
        {
            for ( uint iAdded = 0; iAdded < tSensitivitiesToAdd.n_cols(); iAdded++ )
            {
                aCoordinateSensitivities( iCoordinate, tJoinedSensitivityLength + iAdded ) = tSensitivitiesToAdd( iCoordinate, iAdded );
            }
        }

        // parent coordinate sensitivities, chain factor A * J as in Intersection_Node_Level_Set (MB-022), with
        // d xi* / d x_1 = - 0.5 ( 1 - xi* ) grad phi( x* ) / ( d r / d xi ),  d xi* / d x_2 = - 0.5 ( 1 + xi* ) grad phi( x* ) / ( d r / d xi )
        if ( tFirstParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tDxiDx             = ( -0.5 * ( 1.0 - tXi ) / tDrDxi ) * tGradient;
            Matrix< DDRMat > tLocCoord          = ( 1.0 - tXi ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * aSensitivityFactor * ( tLocCoord + tParentVector * tDxiDx );
            tFirstParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }
        if ( tSecondParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tDxiDx             = ( -0.5 * ( 1.0 + tXi ) / tDrDxi ) * tGradient;
            Matrix< DDRMat > tLocCoord          = ( 1.0 + tXi ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * aSensitivityFactor * ( tLocCoord + tParentVector * tDxiDx );
            tSecondParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Intersection_Node_Analytic::get_coordinate_determining_adv_ids() const
    {
        if ( mLinearFallback )
        {
            return Intersection_Node_Level_Set::get_coordinate_determining_adv_ids();
        }

        // same order as append_dcoordinate_dadv: the field ADVs at the intersection point, then the parents
        Vector< sint > tIDs = mInterfaceGeometry.get_determining_adv_ids_at_point( this->get_global_coordinates() );

        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();
        if ( tFirstParentNode.depends_on_advs() )
        {
            Intersection_Node::join_adv_ids( tIDs, tFirstParentNode.get_coordinate_determining_adv_ids() );
        }
        if ( tSecondParentNode.depends_on_advs() )
        {
            Intersection_Node::join_adv_ids( tIDs, tSecondParentNode.get_coordinate_determining_adv_ids() );
        }
        return tIDs;
    }

    //--------------------------------------------------------------------------------------------------------------
}    // namespace moris::gen
