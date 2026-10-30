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
        const Matrix< DDRMat >& tFirstX  = this->get_first_parent_node().get_global_coordinates();
        const Matrix< DDRMat >& tSecondX = this->get_second_parent_node().get_global_coordinates();
        Matrix< DDRMat >        tEdge    = tSecondX - tFirstX;
        real                    tStep    = 1e-6 * norm( tEdge );

        aGradient = mInterfaceGeometry.get_dfield_dcoordinates_at_point( this->get_global_coordinates(), tStep );

        real tDrDxi = 0.5 * dot( aGradient, tEdge );
        MORIS_ASSERT( std::abs( tDrDxi ) > MORIS_REAL_MIN, "Intersection_Node_Analytic: field tangent to the parent edge at the intersection." );
        return tDrDxi;
    }

    //--------------------------------------------------------------------------------------------------------------

    real Intersection_Node_Analytic::get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const
    {
        // linear-interpolation derivative w.r.t. the end values; not used by this node's sensitivities (which are taken
        // at the intersection point), kept for the interface
        real tThreshold = mInterfaceGeometry.get_isocontour_threshold();
        real tPhi0      = mInterfaceGeometry.get_field_value( this->get_first_parent_node().get_index(), this->get_first_parent_node().get_global_coordinates() );
        real tPhi1      = mInterfaceGeometry.get_field_value( this->get_second_parent_node().get_index(), this->get_second_parent_node().get_global_coordinates() );
        return 2 * ( ( tPhi0 - tThreshold ) * ( aAncestorIndex == 1 ) - ( tPhi1 - tThreshold ) * ( aAncestorIndex == 0 ) ) / std::pow( ( tPhi1 - tPhi0 ), 2 );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Intersection_Node_Analytic::get_dxi_dcoordinate_first_parent() const
    {
        Matrix< DDRMat > tGradient;
        real             tDrDxi = this->get_dresidual_dxi( tGradient );
        return ( -0.5 * ( 1.0 - this->get_local_coordinate() ) / tDrDxi ) * tGradient;
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat > Intersection_Node_Analytic::get_dxi_dcoordinate_second_parent() const
    {
        Matrix< DDRMat > tGradient;
        real             tDrDxi = this->get_dresidual_dxi( tGradient );
        return ( -0.5 * ( 1.0 + this->get_local_coordinate() ) / tDrDxi ) * tGradient;
    }

    //--------------------------------------------------------------------------------------------------------------

    void Intersection_Node_Analytic::append_dcoordinate_dadv(
            Matrix< DDRMat >&       aCoordinateSensitivities,
            const Matrix< DDRMat >& aSensitivityFactor ) const
    {
        const Basis_Node& tFirstParentNode  = this->get_first_parent_node();
        const Basis_Node& tSecondParentNode = this->get_second_parent_node();

        // parent vector (column)
        Matrix< DDRMat > tParentVector = trans( tSecondParentNode.get_global_coordinates() - tFirstParentNode.get_global_coordinates() );

        // d x* / d p = 0.5 ( x_2 - x_1 ) d xi* / d p,  d xi* / d p = - ( d phi / d p )( x* ) / ( d r / d xi )
        Matrix< DDRMat > tGradient;
        real             tDrDxi             = this->get_dresidual_dxi( tGradient );
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

        // parent coordinate sensitivities (chain factor A * J, as in Intersection_Node_Level_Set, MB-022)
        if ( tFirstParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tLocCoord          = ( 1.0 - this->get_local_coordinate() ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * aSensitivityFactor * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_first_parent() );
            tFirstParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }
        if ( tSecondParentNode.depends_on_advs() )
        {
            Matrix< DDRMat > tLocCoord          = ( 1.0 + this->get_local_coordinate() ) * eye( tParentVector.n_rows(), tParentVector.n_rows() );
            Matrix< DDRMat > tSensitivityFactor = 0.5 * aSensitivityFactor * ( tLocCoord + tParentVector * this->get_dxi_dcoordinate_second_parent() );
            tSecondParentNode.append_dcoordinate_dadv( aCoordinateSensitivities, tSensitivityFactor );
        }
    }

    //--------------------------------------------------------------------------------------------------------------

    Vector< sint > Intersection_Node_Analytic::get_coordinate_determining_adv_ids() const
    {
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
