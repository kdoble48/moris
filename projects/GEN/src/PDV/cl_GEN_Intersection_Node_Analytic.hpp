/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Analytic.hpp
 *
 */

#pragma once

#include "cl_GEN_Intersection_Node_Level_Set.hpp"

namespace moris::gen
{
    /**
     * Intersection node on an analytic level set field (use_analytic_intersection): the local coordinate xi* is the
     * root of r( xi ) = phi( x( xi ) ) - c along the parent edge, x( xi ) = 0.5 ( 1 - xi ) x_1 + 0.5 ( 1 + xi ) x_2,
     * found on the analytic field itself (Level_Set_Geometry::compute_intersection_local_coordinate_analytic).
     *
     * Sensitivities by implicit differentiation of r( xi*; p, x_1, x_2 ) = 0 at the intersection point x*:
     *   d xi* / d p   = - ( d phi / d p )( x* ) / ( d r / d xi ),   d r / d xi = 0.5 grad phi( x* ) . ( x_2 - x_1 )
     *   d xi* / d x_1 = - 0.5 ( 1 - xi* ) grad phi( x* ) / ( d r / d xi ),   d xi* / d x_2 = - 0.5 ( 1 + xi* ) grad phi( x* ) / ( d r / d xi )
     * The spatial gradient is a central difference of the analytic field (step 1e-6 of the parent edge length).
     */
    class Intersection_Node_Analytic : public Intersection_Node_Level_Set
    {
      public:
        Intersection_Node_Analytic(
                uint                              aNodeIndex,
                const Vector< Background_Node* >& aBackgroundNodes,
                const Parent_Node&                aFirstParentNode,
                const Parent_Node&                aSecondParentNode,
                mtk::Geometry_Type                aBackgroundGeometryType,
                mtk::Interpolation_Order          aBackgroundInterpolationOrder,
                Level_Set_Geometry&               aInterfaceGeometry );

        void append_dcoordinate_dadv(
                Matrix< DDRMat >&       aCoordinateSensitivities,
                const Matrix< DDRMat >& aSensitivityFactor ) const override;

        Vector< sint > get_coordinate_determining_adv_ids() const override;

      private:
        const Vector< Basis_Node >& get_field_basis_nodes() const override;

        real get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const override;

        Matrix< DDRMat > get_dxi_dcoordinate_first_parent() const override;

        Matrix< DDRMat > get_dxi_dcoordinate_second_parent() const override;

        /**
         * d r / d xi at the intersection point; fills the spatial gradient of the field there (1 x dim).
         */
        real get_dresidual_dxi( Matrix< DDRMat >& aGradient ) const;
    };
}    // namespace moris::gen
