/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_MMA.hpp
 *
 * Native in-tree port of Svanberg's Globally Convergent Method of Moving Asymptotes
 * (GCMMA, 2002). Ported from the canonical reference implementation github.com/kkmaute/gcmma
 * (files mma.cpp/solve.cpp/asymp.cpp/gcmmasub.cpp/raaupdate.cpp/subsolve.cpp/kktcheck.cpp).
 *
 * Owning the outer loop in-tree (instead of linking the external TPL) lets MORIS update the
 * design vector mid-iteration -- e.g. an in-place SDF redistance each step without a GCMMA
 * restart. Variable names mirror the reference for line-by-line validation against the TPL.
 */

#ifndef MORIS_CL_OPT_ALGORITHM_MMA_HPP_
#define MORIS_CL_OPT_ALGORITHM_MMA_HPP_

#include <vector>

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_OPT_Problem.hpp"

namespace moris::opt
{
    class Algorithm_MMA : public Algorithm
    {
      private:
        // ---- user parameters (mirror the GCMMA parameter list) ----
        bool mPrint        = false;
        sint mMaxInnerIter = 0;      ///< max inner conservative cycles (maxinnerit)
        real mNormDrop     = 1e-4;   ///< KKT relative-residual accuracy (acc)
        real mAsympAdapt0  = 0.5;    ///< initial asymptote spread factor (saFac)
        real mAsympShrink  = 0.7;    ///< oscillation shrink factor (sbFac)
        real mAsympExpand  = 1.2;    ///< monotone expand factor (scFac)
        real mStepSize     = 0.01;   ///< move limit fraction (dstep)
        real mPenalty      = 100.0;  ///< constraint penalty c_j (penal; <0 => |f0|)

        // ---- problem dimensions ----
        int mNumVar = 0;   ///< number of design variables (numvar)
        int mNumCon = 0;   ///< number of constraints (numcon)

        // ---- algorithm scalars (reference names) ----
        real epsimin = 0.0;
        real raa0eps = 1.0e-7;
        real a0      = 1.0;
        real raa0    = 0.0;
        real r0      = 0.0;
        real f0val   = 0.0;
        real f0valnew = 0.0;
        real f0app   = 0.0;
        real z = 0.0, zet = 0.0, zmma = 0.0, zetmma = 0.0;

        // ---- 1D state, size numvar ----
        std::vector< double > xval, xold1, xold2, xmax, xmin;
        std::vector< double > low, upp, alfa, beta;
        std::vector< double > p0, q0, df0dx;
        // subsolve working arrays (size numvar, dx is numvar+1)
        std::vector< double > x, xold, delx, diagx, diagxinv, xmma, dx;
        std::vector< double > xsi, xsiold, dxsi, eta, etaold, deta;

        // ---- 1D state, size numcon ----
        std::vector< double > fval, fvalnew, a, b, c, d, r, raa, fapp, gvec;
        std::vector< double > y, yold, dy, dely, diagy, diagyinv;
        std::vector< double > lam, lamold, dlam, dellam, dellamyi;          // dlam is numcon+1
        std::vector< double > diaglam, diaglamyi, diaglamyiinv;
        std::vector< double > mu, muold, dmu, s, sold, ds;

        // ---- 2D state, numcon x numvar ----
        std::vector< std::vector< double > > P, Q, GG, dfdx;

        // ---- dense Newton system (size min(numcon,numvar)+1) ----
        std::vector< double > AA;   ///< flat row-major k x k
        std::vector< double > bb;   ///< rhs, size k

        // ---- optional in-place reinitialization hook ----
        bool mReinitHook = false;   ///< if true, call Problem reinit on xval each outer iteration
        std::vector< double > mLastRedistanced;   ///< redistanced ADVs from the most recent func() eval
                                                  ///< (kept separate so xmma stays the raw subproblem
                                                  ///< solution for raaupdate/conservatism); committed to xval

        // ---- NaN guard ----
        real mMoveScale = 1.0;      ///< move-limit multiplier; shrunk while backtracking off a bad step

        /** true if the current trial criteria (f0valnew, fvalnew) are all finite */
        bool criteria_finite() const;

        // ================= helpers (ports of the reference functions) =================

        /** allocate/size all arrays and set epsimin, a0, a, d (port of initialize.cpp) */
        void mma_initialize();

        /**
         * Evaluate objective+constraints at x via the inherited criteria call (port of func_wrap).
         * aX is mutable: when the in-place reinit hook is active, any redistance the workflow
         * applies to the ADV vector during criteria evaluation is written back into aX, so the
         * redistanced design persists in the optimizer state (no GCMMA restart needed).
         */
        void func( double* aX, real& aF0, std::vector< double >& aF );

        /** evaluate gradients at x via inherited criteria-gradient call (port of grad_wrap) */
        void grad( const double* aX );

        /** update asymptotes low/upp and shrink raa (port of asymp.cpp) */
        void asymp( int aIter );

        /** build conservative convex approximation, solve subproblem, eval f0app/fapp (gcmmasub.cpp) */
        void gcmmasub( int aIter );

        /** increase raa0/raa until conservative (raaupdate.cpp) */
        void raaupdate();

        /** primal-dual interior-point solve of the MMA subproblem (subsolve.cpp) */
        void subsolve();

        /** KKT residual convergence test; returns 1 if converged (kktcheck.cpp) */
        int kktcheck( real& aKKTrefres );

        /** small dense LU solve with partial pivoting: solves A(k x k) x = b in place of dgesv */
        static void solve_dense( std::vector< double >& aA, std::vector< double >& aX,
                const std::vector< double >& aB, int aN );

        inline double mx( double a, double b ) { return a > b ? a : b; }
        inline double mn( double a, double b ) { return a < b ? a : b; }

      public:
        /** Constructor */
        explicit Algorithm_MMA( const moris::Parameter_List& aParameterList );

        /** Destructor */
        ~Algorithm_MMA() override;

        /**
         * @brief Solve the optimization problem with the native GCMMA
         */
        uint solve(
                uint                                   aCurrentOptAlgInd,
                std::shared_ptr< moris::opt::Problem > aOptProb ) override;

        /** Run the GCMMA outer loop on proc 0 (port of solve.cpp) */
        int mma_solve();
    };
}    // namespace moris::opt

#endif /* MORIS_CL_OPT_ALGORITHM_MMA_HPP_ */
