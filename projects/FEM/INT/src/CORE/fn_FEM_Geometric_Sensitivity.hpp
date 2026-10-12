/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Geometric_Sensitivity.hpp
 *
 * Helpers for analytical geometric sensitivities, i.e. the derivatives of residuals and
 * quantities of interest wrt the coordinates of the integration (IG) cell vertices.
 *
 * Moving IG vertex v in direction i by dX changes, at a fixed IG integration point eta:
 *   (1) the integration measure:  d(detJ) = detJ * m(v,i) dX                    (see eval_ig_measure_derivatives)
 *   (2) the evaluation point:     d(x_gp) = N_v(eta) e_i dX, so every field quantity
 *       f(x_gp) changes by        df      = N_v(eta) * df/dx_i dX               (affine IP element map)
 *   (3) for side sets the normal: dn      = dn/dX(v,i) dX                       (see eval_ig_normal_derivatives)
 * which is what the finite difference geometry path (IWG::select_dRdp_FD_geometry_bulk/sideset) perturbs.
 */

#pragma once

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"

namespace moris::fem
{
    class Geometry_Interpolator;
    class Field_Interpolator;

    //------------------------------------------------------------------------------
    /**
     * relative derivative of the IG integration measure wrt the IG vertex coordinates
     * @param[ in ]  aIGGI     IG geometry interpolator, evaluation point set
     * @param[ in ]  aIsSide   true if the IG cell is a side (line in 2D, triangle/quad facet in 3D)
     * @param[ out ] adMeasure ( number of IG vertices x number of space dimensions ),
     *                         adMeasure( v, i ) = ( d detJ / d X( v, i ) ) / detJ
     */
    void eval_ig_measure_derivatives(
            Geometry_Interpolator* aIGGI,
            bool                   aIsSide,
            Matrix< DDRMat >&      adMeasure );

    //------------------------------------------------------------------------------
    /**
     * derivative of the unit outward normal of an IG side wrt the IG side vertex coordinates
     * @param[ in ]  aIGGI    IG geometry interpolator of the side, evaluation point set
     * @param[ in ]  aNormal  unit normal used by the IWG/IQI (orientation reference)
     * @param[ out ] adNormal one matrix per space direction i ( number of space dimensions x number of IG vertices ),
     *                        column v of adNormal( i ) = d n / d X( v, i )
     */
    void eval_ig_normal_derivatives(
            Geometry_Interpolator*      aIGGI,
            const Matrix< DDRMat >&     aNormal,
            Vector< Matrix< DDRMat > >& adNormal );

    //------------------------------------------------------------------------------
    /**
     * row of the second spatial derivative matrix ( Field_Interpolator::dnNdxn( 2 ) ) that holds d2/dx_j dx_k
     * ordering 2D: xx, yy, xy;  3D: xx, yy, zz, yz, xz, xy
     */
    uint second_derivative_row(
            uint aJ,
            uint aK,
            uint aSpaceDim );

    //------------------------------------------------------------------------------
    /**
     * derivative of the first spatial derivatives of the shape functions wrt x_k
     * @param[ in ]  aFI         field interpolator, evaluation point set
     * @param[ in ]  aDirection  spatial direction k
     * @param[ out ] adNdxdxk    ( number of space dimensions x number of bases ), row j = d2N / dx_j dx_k
     */
    void eval_dNdx_spatial_derivative(
            Field_Interpolator* aFI,
            uint                aDirection,
            Matrix< DDRMat >&   adNdxdxk );

    //------------------------------------------------------------------------------
    /**
     * build the vector-field shape function matrix of a scalar shape-function row, in the layout of
     * Field_Interpolator::N() ( number of fields x number of fields * number of bases )
     * @param[ in ]  aNRow     ( 1 x number of bases ) shape function values or derivatives
     * @param[ in ]  aNumFields number of fields
     * @param[ out ] aN        block diagonal matrix
     */
    void build_vector_shape_function_matrix(
            const Matrix< DDRMat >& aNRow,
            uint                    aNumFields,
            Matrix< DDRMat >&       aN );

    //------------------------------------------------------------------------------
    /**
     * vectorize the field coefficients ( number of bases x number of fields ) column by column,
     * consistent with the layout of Field_Interpolator::N() and the structural test strain
     */
    void vectorize_coefficients(
            const Matrix< DDRMat >& aCoeff,
            Matrix< DDRMat >&       aVec );

    //------------------------------------------------------------------------------
    /**
     * true if at least one entry of the pdv local assembly map is active ( != -1 )
     */
    bool has_active_geometry_pdv( const Matrix< DDSMat >& aGeoLocalAssembly );

    //------------------------------------------------------------------------------
}    // namespace moris::fem
