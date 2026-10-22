/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Summed_Fields.hpp
 *
 */

#pragma once

#include "cl_GEN_Field.hpp"

namespace moris::gen
{
    /**
     * A field whose value is the sum of its child fields. ADV sensitivities and determining ADV IDs are the
     * concatenation of the children's (in child order); coordinate sensitivities are summed.
     */
    class Summed_Fields : public Field
    {

      private:
        Vector< std::shared_ptr< Field > > mFields;

      public:
        /**
         * Summed field constructor
         *
         * @param aFields Fields to sum (at least one)
         * @param aName Name of this field
         */
        Summed_Fields(
                const Vector< std::shared_ptr< Field > >& aFields,
                std::string                               aName = "" );

        /**
         * Imports the local ADVs required from the full owned ADV distributed vector, for every child field.
         *
         * @param aOwnedADVs Full owned distributed ADV vector
         */
        void import_advs( sol::Dist_Vector* aOwnedADVs ) override;

        /**
         * Resets all nodal information of every child field.
         *
         * @param aMesh Interpolation mesh
         */
        void reset_nodal_data( mtk::Interpolation_Mesh* aMesh ) override;

        /**
         * Gets whether any child field has ADVs.
         *
         * @return true if at least one child has ADVs
         */
        bool has_advs() override;

        /**
         * Gets the coefficient ADV IDs of every child, concatenated in child order.
         *
         * @return Coefficient ADV IDs
         */
        Vector< sint > get_coefficient_adv_ids() override;

        /**
         * Forwards the updated dependencies to every child field.
         *
         * @param aUpdatedFields Other fields that the children may depend on
         */
        void update_dependencies( const Vector< std::shared_ptr< Field > >& aUpdatedFields ) override;

        /**
         * Gets the sum of the child field values at a node.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Summed field value
         */
        real get_field_value(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets the sum of the child field values at a derived node.
         *
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         * @return Summed field value
         */
        real get_field_value(
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) override;

        /**
         * Gets the child field sensitivities with respect to ADVs, concatenated in child order.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @return Vector of sensitivities
         */
        const Matrix< DDRMat >& get_dfield_dadvs(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets the child field sensitivities with respect to ADVs at a derived node, concatenated in child order.
         *
         * @param aSensitivities Sensitivities to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        void get_dfield_dadvs(
                Matrix< DDRMat >&   aSensitivities,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) override;

        /**
         * Gets the sum of the child field sensitivities with respect to the nodal coordinates.
         * NOTE: a B-spline (discrete) child does not support coordinate sensitivities.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Vector of coordinate values
         * @param aSensitivities Sensitivities to be filled with d(field value)/d(coordinate_j)
         */
        void get_dfield_dcoordinates(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates,
                Matrix< DDRMat >&       aSensitivities ) override;

        /**
         * Gets the IDs of ADVs which this field depends on, concatenated in child order.
         *
         * @param aNodeIndex Node index
         * @param aCoordinates Node coordinates
         * @return Determining ADV IDs at this node
         */
        Vector< sint > get_determining_adv_ids(
                uint                    aNodeIndex,
                const Matrix< DDRMat >& aCoordinates ) override;

        /**
         * Gets the IDs of ADVs that this field depends on at a derived node, concatenated in child order.
         *
         * @param aDeterminingADVIDs Determining ADV IDs to fill for the given derived node
         * @param aDerivedNode Derived node
         * @param aNodeManager Node manager
         */
        void get_determining_adv_ids(
                Vector< sint >&     aDeterminingADVIDs,
                const Derived_Node& aDerivedNode,
                const Node_Manager& aNodeManager ) override;

        /**
         * Remeshing of a summed design is unsupported, so this always returns nullptr.
         *
         * @return nullptr
         */
        std::shared_ptr< mtk::Field > get_mtk_field() override;
    };
}
