/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_GEN_Intersection_Node_Bilinear.cpp
 *
 */

#include "cl_GEN_Intersection_Node_Bilinear.hpp"
#include "cl_GEN_Parent_Node.hpp"
#include "cl_GEN_Level_Set_Geometry.hpp"
#include "cl_GEN_Interpolation.hpp"

#include "cl_MTK_Interpolation_Function_Base.hpp"
#include "cl_MTK_Interpolation_Function_Factory.hpp"
#include "cl_MTK_Enums.hpp"

namespace moris::gen
{
    //--------------------------------------------------------------------------------------------------------------

    Intersection_Node_Bilinear::Intersection_Node_Bilinear(
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
            , mParametricParentVector( aSecondParentNode.get_parametric_coordinates() - aFirstParentNode.get_parametric_coordinates() )
    {
    }

    //--------------------------------------------------------------------------------------------------------------

    const Vector< Basis_Node >& Intersection_Node_Bilinear::get_field_basis_nodes() const
    {
        return this->get_background_nodes();
    }

    //--------------------------------------------------------------------------------------------------------------

    mtk::Interpolation_Function_Base*
    Intersection_Node_Bilinear::create_cell_interpolation( uint aNumDims, uint aNumBases )
    {
        mtk::Geometry_Type       tGeometryType = aNumDims == 2 ? mtk::Geometry_Type::QUAD : mtk::Geometry_Type::HEX;
        mtk::Interpolation_Order tOrder        = mtk::Interpolation_Order::UNDEFINED;

        if ( ( aNumDims == 2 and aNumBases == 4 ) or ( aNumDims == 3 and aNumBases == 8 ) )
        {
            tOrder = mtk::Interpolation_Order::LINEAR;
        }
        else if ( ( aNumDims == 2 and aNumBases == 9 ) or ( aNumDims == 3 and aNumBases == 27 ) )
        {
            tOrder = mtk::Interpolation_Order::QUADRATIC;
        }
        else if ( ( aNumDims == 2 and aNumBases == 16 ) or ( aNumDims == 3 and aNumBases == 64 ) )
        {
            tOrder = mtk::Interpolation_Order::CUBIC;
        }

        MORIS_ERROR( tOrder != mtk::Interpolation_Order::UNDEFINED,
                "Intersection_Node_Bilinear::create_cell_interpolation - no Lagrange QUAD/HEX cell with %u nodes in %uD.",
                aNumBases,
                aNumDims );

        mtk::Interpolation_Function_Factory tFactory;
        return tFactory.create_interpolation_function( tGeometryType, mtk::Interpolation_Type::LAGRANGE, tOrder );
    }

    //--------------------------------------------------------------------------------------------------------------

    real
    Intersection_Node_Bilinear::get_dxi_dfield_from_ancestor( uint aAncestorIndex ) const
    {
        // Cell-order intersection: the root of the cell's own-order interpolant phi(xi) = sum_k N_k(xi) phi_k along the
        // edge, so d(edge coordinate)/d(phi_k) = -N_k(xi*) / (d phi / d edge coordinate) for every cell node k
        if ( mInterfaceGeometry.use_cell_order_intersection() )
        {
            const Vector< Basis_Node >& tCellNodes = this->get_background_nodes();
            uint                        tNumBases  = tCellNodes.size();

            mtk::Interpolation_Function_Base* tCellInterpolation =
                    create_cell_interpolation( this->get_global_coordinates().length(), tNumBases );

            Matrix< DDRMat > tPhiCellNodes( tNumBases, 1 );
            for ( uint iNode = 0; iNode < tNumBases; ++iNode )
            {
                tPhiCellNodes( iNode ) = mInterfaceGeometry.get_field_value(
                        tCellNodes( iNode ).get_index(),
                        tCellNodes( iNode ).get_global_coordinates() );
            }

            const Matrix< DDRMat >& tCellCoordinate = this->get_parametric_coordinates();

            Matrix< DDRMat > tCellBasis;
            tCellInterpolation->eval_N( tCellCoordinate, tCellBasis );

            Matrix< DDRMat > tCellDBasisDxi;
            tCellInterpolation->eval_dNdXi( tCellCoordinate, tCellDBasisDxi );
            Matrix< DDRMat > tCellJacobian = 0.5 * mParametricParentVector * tCellDBasisDxi * tPhiCellNodes;

            delete tCellInterpolation;

            return -1.0 * tCellBasis( aAncestorIndex ) / ( tCellJacobian( 0 ) + MORIS_REAL_EPS );
        }

        // number of nodes to be used for interpolation
        uint tNumBases = 0;

        // determine number of basis for bi or tri-linear interpolation based on spatial dimensions
        switch ( this->get_global_coordinates().length() )
        {
            case 2:
            {
                tNumBases = 4;
                break;
            }
            case 3:
            {
                tNumBases = 8;
                break;
            }
            default:
            {
                MORIS_ERROR( false,
                        "Intersection_Node_Bilinear::compute_intersection_derivative - Improper spatial dimensions." );
            }
        }

        // check that aAncestorIndex <= number of basis
        // note: here only bi and tri-linear interpolation used irrespective of interpolation of background cell; thus only
        // corner nodes values are used; level set values of other nodes do not influence intersection position
        if ( aAncestorIndex >= tNumBases )
        {
            // return zero as only corner node level set values influence intersection
            return 0.0;
        }

        // build interpolator
        mtk::Interpolation_Function_Factory tFactory;
        mtk::Interpolation_Function_Base*   tInterpolation;

        // create interpolation function based on spatial dimension  of problem
        switch ( tNumBases )
        {
            case 4:
            {
                tInterpolation = tFactory.create_interpolation_function(
                        mtk::Geometry_Type::QUAD,
                        mtk::Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::LINEAR );
                break;
            }
            case 8:
            {
                tInterpolation = tFactory.create_interpolation_function(
                        mtk::Geometry_Type::HEX,
                        mtk::Interpolation_Type::LAGRANGE,
                        mtk::Interpolation_Order::LINEAR );
                break;
            }
            default:
            {
                MORIS_ERROR( false,
                        "Intersection_Node_Bilinear::compute_intersection_derivative - Interpolation type not implemented." );
            }
        }

        // allocate matrix for level set values at background cell nodes
        Matrix< DDRMat > tPhiBCNodes( tNumBases, 1 );

        // get level set values of corner nodes
        const Vector< Basis_Node >& tBackgroundNodes = this->get_background_nodes();
        for ( uint iBackgroundNodeIndex = 0; iBackgroundNodeIndex < tNumBases; ++iBackgroundNodeIndex )
        {
            tPhiBCNodes( iBackgroundNodeIndex ) = mInterfaceGeometry.get_field_value(
                    tBackgroundNodes( iBackgroundNodeIndex ).get_index(),
                    tBackgroundNodes( iBackgroundNodeIndex ).get_global_coordinates() );
        }

        // compute local coordinate in background cell CS
        const Matrix< DDRMat >& tCellCoordinate = this->get_parametric_coordinates();

        // compute basis function
        Matrix< DDRMat > tBasis;
        tInterpolation->eval_N( tCellCoordinate, tBasis );

        // compute derivative of residual
        real tDResidualDPhi = tBasis( aAncestorIndex );

        // compute Jacobian
        Matrix< DDRMat > tDBasisDxi;
        tInterpolation->eval_dNdXi( tCellCoordinate, tDBasisDxi );
        Matrix< DDRMat > tJac = 0.5 * mParametricParentVector * tDBasisDxi * tPhiBCNodes;

        // delete interpolator
        delete tInterpolation;

        // compute derivative of edge coordinate wrt. ancestor level set value
        return -1.0 * tDResidualDPhi / ( tJac( 0 ) + MORIS_REAL_EPS );
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Bilinear::get_dxi_dcoordinate_first_parent() const
    {
        MORIS_ERROR( false, "Intersections on intersections not implemented yet for bilinear case." );
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

    Matrix< DDRMat >
    Intersection_Node_Bilinear::get_dxi_dcoordinate_second_parent() const
    {
        MORIS_ERROR( false, "Intersections on intersections not implemented yet for bilinear case." );
        return { {} };
    }

    //--------------------------------------------------------------------------------------------------------------

}    // namespace moris::gen
