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
 *
 * Properties enter (2) through their spatial variation dP/dx_i. The analytical path is used only if that variation is
 * known for every property involved (see classify_property_spatial_variation); otherwise the IWG/IQI reports no
 * analytical support and the FD geometry path is used (MB-062).
 */

#pragma once

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"

#include <memory>

namespace moris::fem
{
    class Geometry_Interpolator;
    class Field_Interpolator;
    class Property;

    //------------------------------------------------------------------------------
    /**
     * spatial variation of a property at the evaluation point, as seen by the analytical geometric sensitivities
     *   CONSTANT       - no dof, dv or field dependency, no space derivative functions, and either the default constant
     *                    value function or a value function whose value does not change when the evaluation point is
     *                    translated by +/- a fraction of the IG cell size in each direction (parameter-only functions)
     *   DIFFERENTIABLE - no dof, dv or field dependency, space derivative functions are set ( Property::dPropdx( 1 ) )
     *   UNKNOWN        - anything else (dof, dv or field dependent, or space dependent without a space derivative):
     *                    the analytical geometric sensitivities are not available, use the FD geometry path
     */
    enum class Property_Spatial_Variation
    {
        CONSTANT,
        DIFFERENTIABLE,
        UNKNOWN
    };

    //------------------------------------------------------------------------------
    /**
     * classify the spatial variation of a property at the current evaluation point (see Property_Spatial_Variation)
     * @param[ in ] aProperty property, the field interpolator manager is set at the evaluation point; nullptr is CONSTANT
     */
    Property_Spatial_Variation classify_property_spatial_variation( const std::shared_ptr< Property >& aProperty );

    //------------------------------------------------------------------------------
    /**
     * true if the spatial variation of the property is known (CONSTANT or DIFFERENTIABLE), or the property is nullptr
     */
    bool has_known_spatial_variation( const std::shared_ptr< Property >& aProperty );

    //------------------------------------------------------------------------------
    /**
     * true if the property is nullptr or CONSTANT
     */
    bool is_spatially_constant( const std::shared_ptr< Property >& aProperty );

    //------------------------------------------------------------------------------
    /**
     * true if the property is nullptr, CONSTANT, or DIFFERENTIABLE with a nonzero value:
     * a thickness that scales the whole contribution, whose derivative is added as R * ( dt/dx ) / t
     */
    bool has_known_spatial_variation_nonzero( const std::shared_ptr< Property >& aProperty );

    //------------------------------------------------------------------------------
    /**
     * derivative of the property value wrt the evaluation point location x_k
     * precondition: has_known_spatial_variation( aProperty ) (checked by the supports_analytical_* hooks)
     * @param[ in ]  aProperty  property, nullptr or of known spatial variation
     * @param[ in ]  aDirection spatial direction k
     * @param[ out ] adPropdxk  ( number of entries of the property value x 1 ), zero if the property is constant
     *                          (or nullptr: 1 x 1 zero). The space derivative function returns, like Field_Interpolator::gradx,
     *                          a ( number of space dimensions x number of entries ) matrix; for a scalar property a vector
     *                          of length number of space dimensions is accepted.
     */
    void eval_property_space_derivative(
            const std::shared_ptr< Property >& aProperty,
            uint                               aDirection,
            Matrix< DDRMat >&                  adPropdxk );

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
