/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * fn_linsolve.hpp
 *
 */

#ifndef PROJECTS_LINALG_SRC_FN_LINSOLVE_HPP_
#define PROJECTS_LINALG_SRC_FN_LINSOLVE_HPP_

// MORIS library header files.
#include "cl_Matrix.hpp"

#ifdef MORIS_USE_EIGEN
#include "Eigen_Impl/fn_linsolve_Eigen.hpp"
#endif

#ifdef MORIS_USE_ARMA
#include "fn_linsolve_Arma.hpp"
#endif

namespace moris
{
    /**
     * @brief Solve for a linear set of equations Ax = B.
     *
     * @param[in] A The LHS Matrix
     * @param[in] B The RHS Vector
     * @param[in] aSolver Optional solver type. NOTE: currently IGNORED by both
     *            backends (kept for API compatibility). Both debug and opt
     *            builds use the backend's default (accuracy-checked) algorithm.
     *
     * @return The vector of solutions, x. Similar to B/A in Matlab
     *
     * @note If A is square, solve() is faster and more accurate than using X = inv(A)*B .
     * If A is non-square, solve() will try to provide approximate solutions to under-determined
     * as well as over-determined systems.
     * Eigen provides various options for decomposition of matrices to facilitate a linear solve.
     * We use the default option of QR decomposition with column pivoting.
     *
     */
    template< typename Matrix_Type >
    auto
    solve( Matrix< Matrix_Type > const &  aA,
            Matrix< Matrix_Type > const & aB,
            std::string const &           aSolver = "default" )
            -> decltype( solve( aA.matrix_data(), aB.matrix_data() ) )
    {
        MORIS_ASSERT( aA.n_rows() > 10, "For matrices smaller than 10x10 use inv() instead of solve().\n" );

        return solve( aA.matrix_data(), aB.matrix_data() );
    }

    /**
     * @brief Dense least-squares solve into a caller-provided output:
     *        aX = argmin || aA aX - aB ||_2 for a tall ( m >= k )
     *        full-column-rank system, via the backend's compiled LAPACK
     *        path (dgels for the Armadillo backend).
     *
     * Intended for REPEATED small/medium dense solves: writes into aX
     * (resized only when needed - reuse the same output across calls to
     * avoid per-call allocation of the result), and never falls back to an
     * approximate (pseudo-inverse) solution.
     *
     * @param[in]  aA tall dense matrix ( m x k, m >= k )
     * @param[in]  aB right-hand side ( m x nrhs )
     * @param[out] aX solution ( k x nrhs ), written in place
     *
     * @return true on success; false if the backend flags the system as
     *         (numerically) rank deficient or the solve fails - the caller
     *         must handle that case (no approximate solution is returned).
     *
     * @note Only implemented for the Armadillo backend; the Eigen build
     *       returns false so callers exercise their fallback path.
     */
    inline bool
    solve_least_squares(
            Matrix< DDRMat > const & aA,
            Matrix< DDRMat > const & aB,
            Matrix< DDRMat >&        aX )
    {
#ifdef MORIS_USE_ARMA
        return arma::solve(
                aX.matrix_data(),
                aA.matrix_data(),
                aB.matrix_data(),
                arma::solve_opts::fast + arma::solve_opts::no_approx );
#else
        return false;
#endif
    }

    /**
     * @brief Workspace variant of solve_least_squares: solves on the LEADING
     *        aNumCols columns of aA (backend submatrix view - no copy of aA).
     *
     * Lets a caller keep one geometrically-grown gather buffer and re-solve
     * subproblems of varying column count without any per-call allocation of
     * the system matrix (active-set methods, column-selection loops).
     *
     * @param[in]  aA       workspace matrix ( m x >= aNumCols )
     * @param[in]  aNumCols number of leading columns forming the system
     * @param[in]  aB       right-hand side ( m x nrhs )
     * @param[out] aX       solution ( aNumCols x nrhs ), written in place
     *
     * @return true on success (see solve_least_squares)
     */
    inline bool
    solve_least_squares(
            Matrix< DDRMat > const & aA,
            uint                     aNumCols,
            Matrix< DDRMat > const & aB,
            Matrix< DDRMat >&        aX )
    {
        MORIS_ASSERT( aNumCols > 0 && aNumCols <= aA.n_cols(),
                "solve_least_squares - invalid leading column count." );

#ifdef MORIS_USE_ARMA
        return arma::solve(
                aX.matrix_data(),
                aA.matrix_data().head_cols( aNumCols ),
                aB.matrix_data(),
                arma::solve_opts::fast + arma::solve_opts::no_approx );
#else
        return false;
#endif
    }
}    // namespace moris

#endif /* PROJECTS_LINALG_SRC_FN_LINSOLVE_HPP_ */
