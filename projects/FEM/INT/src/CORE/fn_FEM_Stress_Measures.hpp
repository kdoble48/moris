/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Stress_Measures.hpp
 *
 * Scalar stress measures of a constitutive-model flux and their derivatives wrt the flux components; the single
 * implementation shared by the stress IQIs (IQI_Stress and its children) and the stress IWGs (IWG_Struc_Stress).
 *
 * The flux is the stress vector of an elasticity CM in its own layout:
 *     3 components: 2D plane stress            ( s_xx, s_yy, s_xy )
 *     4 components: 2D plane strain / axisym.  ( s_xx, s_yy, s_zz, s_xy )
 *     6 components: 3D                         ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy )
 * Internally every measure works on the standardized 6-component vector ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy ),
 * missing components zero. Derivatives are returned in the flux layout, so the chain rule to the dofs is
 *     ds/du = trans( dFluxdDOF ) * ds/dflux   ( ndof x 1 ).
 */

#pragma once

#include "moris_typedefs.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_FEM_Enums.hpp"

namespace moris::fem
{
    //------------------------------------------------------------------------------
    /**
     * sort a flux into the standardized 6-component layout ( s_xx, s_yy, s_zz, s_yz, s_xz, s_xy )
     * @param[ in ] aFlux CM flux with 3, 4 or 6 components
     * @return standardized stress vector ( 6 x 1 )
     */
    Matrix< DDRMat > standardize_stress_vector( const Matrix< DDRMat >& aFlux );

    //------------------------------------------------------------------------------
    /**
     * map a vector in the standardized layout back to the flux layout (inverse of standardize_stress_vector;
     * components that the flux does not carry are dropped)
     * @param[ in ] aStandardized       vector in the standardized layout ( 6 x 1 )
     * @param[ in ] aNumFluxComponents  number of flux components ( 3, 4 or 6 )
     */
    Matrix< DDRMat > stress_vector_to_flux_layout(
            const Matrix< DDRMat >& aStandardized,
            uint                    aNumFluxComponents );

    //------------------------------------------------------------------------------
    /**
     * von Mises stress
     *     sigma_vm = sqrt( 1/2 [ ( s_xx - s_yy )^2 + ( s_yy - s_zz )^2 + ( s_zz - s_xx )^2 ] + 3 ( s_yz^2 + s_xz^2 + s_xy^2 ) )
     */
    real eval_von_mises_stress( const Matrix< DDRMat >& aFlux );

    //------------------------------------------------------------------------------
    /**
     * derivative of the von Mises stress wrt the flux components, in the flux layout
     *     d sigma_vm / d s_ii = ( 2 s_ii - s_jj - s_kk ) / ( 2 sigma_vm ),   d sigma_vm / d s_ij = 3 s_ij / sigma_vm
     * sigma_vm is not differentiable at a stress-free point; there ( sigma_vm <= MORIS_REAL_EPS ) the derivative is zero.
     */
    Matrix< DDRMat > eval_dvon_mises_stress_dflux( const Matrix< DDRMat >& aFlux );

    //------------------------------------------------------------------------------
    /**
     * i-th principal stress, sorted descending ( 1 = largest )
     * @param[ in ] aFlux                  CM flux
     * @param[ in ] aPrincipalStressIndex  1 or 2 in 2D, 1, 2 or 3 in 3D
     * @param[ in ] aNumSpaceDims          number of space dimensions ( 2: in-plane principal stresses )
     */
    real eval_principal_stress(
            const Matrix< DDRMat >& aFlux,
            uint                    aPrincipalStressIndex,
            uint                    aNumSpaceDims );

    //------------------------------------------------------------------------------
    /**
     * i-th normal stress ( 1, 2, 3 = s_xx, s_yy, s_zz )
     */
    real eval_normal_stress( const Matrix< DDRMat >& aFlux, uint aStressIndex );

    //------------------------------------------------------------------------------
    /**
     * i-th shear stress ( 1, 2, 3 = s_yz, s_xz, s_xy; in 2D only 3 exists, and 1 or 2 is an error )
     */
    real eval_shear_stress( const Matrix< DDRMat >& aFlux, uint aStressIndex );

    //------------------------------------------------------------------------------
    /**
     * scalar stress measure selected by a stress type
     * @param[ in ] aStressType     VON_MISES_STRESS, PRINCIPAL_STRESS, NORMAL_STRESS or SHEAR_STRESS
     * @param[ in ] aFlux           CM flux
     * @param[ in ] aIndex          1-based component / principal index (ignored for the von Mises stress)
     * @param[ in ] aNumSpaceDims   number of space dimensions
     */
    real eval_stress_measure(
            Stress_Type             aStressType,
            const Matrix< DDRMat >& aFlux,
            uint                    aIndex,
            uint                    aNumSpaceDims );

    //------------------------------------------------------------------------------
    /**
     * derivative of the scalar stress measure wrt the flux components, in the flux layout
     * (implemented for the von Mises stress; errors for the other stress types)
     */
    Matrix< DDRMat > eval_dstress_measure_dflux(
            Stress_Type             aStressType,
            const Matrix< DDRMat >& aFlux );

    //------------------------------------------------------------------------------
}    // namespace moris::fem
