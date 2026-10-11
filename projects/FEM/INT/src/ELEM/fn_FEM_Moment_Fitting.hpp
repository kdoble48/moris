/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_FEM_Moment_Fitting.hpp
 *
 * Moment-fitted cut-cell quadrature (Phase 1).
 *
 * For a cut cell cluster the material region is tessellated into simplices (XTK:
 * TETs in 3D, TRIs in 2D). This module fits a single non-negative quadrature rule
 * on the cluster's interpolation (parent) cell that reproduces the total-degree-d
 * moments of the material region:
 *
 *      A w = b ,   w >= 0
 *
 * where the candidate points (columns of A) are the cluster's EXISTING
 * tessellation Gauss points mapped to parent coordinates, A holds the
 * total-degree-d Legendre-product basis evaluated at the candidates, and b holds
 * the exact moments of the material region computed on the tessellation with a
 * degree-d-exact simplex rule. The non-negative least-squares (NNLS) active set
 * prunes the candidate set; a per-cluster fallback to the tessellated rule is
 * triggered when the relative moment residual exceeds a tolerance.
 *
 */

#ifndef SRC_FEM_FN_FEM_MOMENT_FITTING_HPP_
#define SRC_FEM_FN_FEM_MOMENT_FITTING_HPP_

#include "moris_typedefs.hpp"    // MRS/COR/src
#include "cl_Matrix.hpp"         // LINALG/src
#include "linalg_typedefs.hpp"
#include "cl_Vector.hpp"    // MRS/CNT/src

namespace moris::fem::moment_fitting
{
    //------------------------------------------------------------------------------
    /**
     * number of basis functions of the total-degree-d polynomial space
     * @param[ in ] aDim    spatial dimension ( 2 or 3 )
     * @param[ in ] aDegree total polynomial degree d
     * @return 2D: ( d + 1 )( d + 2 ) / 2 ; 3D: ( d + 1 )( d + 2 )( d + 3 ) / 6
     */
    uint number_of_basis_functions( uint aDim, uint aDegree );

    //------------------------------------------------------------------------------
    /**
     * evaluate the total-degree-d Legendre-product basis (2D or 3D from the row
     * count of aPoints)
     * ordering: by total degree, then lexicographic in the exponents,
     * 3D: for tot = 0...d : for i = 0...tot : for j = 0...tot-i : k = tot-i-j
     * 2D: for tot = 0...d : for i = 0...tot : j = tot-i
     * @param[ in ]  aPoints points ( dim x n ) in [-1,1]^dim (parent coordinates)
     * @param[ in ]  aDegree total polynomial degree d
     * @param[ out ] aBasis  basis values ( m x n ) with m = number_of_basis_functions( dim, d )
     */
    void evaluate_legendre_basis(
            const Matrix< DDRMat >& aPoints,
            uint                    aDegree,
            Matrix< DDRMat >&       aBasis );

    //------------------------------------------------------------------------------
    /**
     * non-negative least squares by the Lawson-Hanson active-set algorithm
     * (small dense systems; the passive-set solves use modified Gram-Schmidt QR)
     * @param[ in ]  aA       matrix ( m x n )
     * @param[ in ]  aB       right-hand side ( m x 1 )
     * @param[ out ] aX       solution ( n x 1 ), all entries >= 0
     * @param[ in ]  aMaxIter maximum number of active-set iterations ( 0: 10 * n )
     * @return || A x - b ||_2
     */
    real nnls(
            const Matrix< DDRMat >& aA,
            const Matrix< DDRMat >& aB,
            Matrix< DDRMat >&       aX,
            uint                    aMaxIter = 0 );

    //------------------------------------------------------------------------------
    /**
     * affine map of a linear simplex given by its vertex coordinates in some
     * target frame, using the MORIS parametrizations
     * TET ( 4 x 3 ): N = [ z1, z2, 1-z1-z2-z3, z3 ], vertices at
     *                (1,0,0), (0,1,0), (0,0,0), (0,0,1)
     * TRI ( 3 x 2 ): N = [ z1, z2, 1-z1-z2 ], vertices at (1,0), (0,1), (0,0)
     *      xi( zeta ) = aOrigin + aM * zeta
     * @param[ in ]  aCellCoords vertex coordinates ( nv x dim ), rows are vertices
     * @param[ out ] aM          map matrix ( dim x dim )
     * @param[ out ] aOrigin     map origin ( dim x 1 ) = vertex 2
     * @param[ out ] aDetJ       det( aM ) / dim! (MORIS simplex space_det_J
     *                           convention: weights summing to one integrate the
     *                           simplex measure)
     */
    void simplex_affine_map(
            const Matrix< DDRMat >& aCellCoords,
            Matrix< DDRMat >&       aM,
            Matrix< DDRMat >&       aOrigin,
            real&                   aDetJ );

    //------------------------------------------------------------------------------
    /**
     * input for the per-cluster fit
     */
    struct Fit_Input
    {
        // per simplex: vertex coordinates in the parent (IP) parametric frame
        // ( nv x dim : 4 x 3 for TET, 3 x 2 for TRI )
        Vector< Matrix< DDRMat > > mCellParamCoords;

        // candidate rule on the reference simplex: the set's existing space rule
        // ( dim x nGP, 1 x nGP )
        Matrix< DDRMat > mCandSpacePoints;
        Matrix< DDRMat > mCandSpaceWeights;

        // moment rule on the reference simplex: degree-d exact ( dim x nGM, 1 x nGM )
        Matrix< DDRMat > mMomSpacePoints;
        Matrix< DDRMat > mMomSpaceWeights;

        // total polynomial degree of the fitting basis
        uint mDegree = 4;

        // relative moment residual tolerance for accepting the fit
        real mRelTol = 1e-10;
    };

    //------------------------------------------------------------------------------
    /**
     * result of the per-cluster fit
     */
    struct Fit_Result
    {
        // retained (pruned) quadrature points in parent coordinates ( dim x nRet )
        Matrix< DDRMat > mPoints;

        // weights carrying the material measure in parent coordinates ( 1 x nRet )
        Matrix< DDRMat > mWeights;

        // relative moment residual || A w - b || / || b ||
        real mResidual = -1.0;

        // number of candidate points offered to the fit
        uint mNumCandidates = 0;

        // material measure in parent coordinates (zeroth moment)
        real mMaterialVolume = 0.0;

        // true if the fit was accepted
        bool mSuccess = false;
    };

    //------------------------------------------------------------------------------
    /**
     * fit a single non-negative quadrature rule for one cut cell cluster
     * (see file header); fails (aResult.mSuccess = false) if a simplex is
     * inverted, the material measure vanishes, or the relative moment residual
     * exceeds the tolerance - the caller then keeps the tessellated rule
     */
    Fit_Result fit_cluster_rule( const Fit_Input& aInput );

    //------------------------------------------------------------------------------
}    // namespace moris::fem::moment_fitting

#endif /* SRC_FEM_FN_FEM_MOMENT_FITTING_HPP_ */
