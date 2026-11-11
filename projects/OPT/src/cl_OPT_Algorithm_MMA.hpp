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

#include <string>
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
        bool mShiftHistory = true;   ///< reinit mode: shift xold1/xold2/low/upp by the redistance (warm restart)
        bool mReinitOnAcceptOnly = false;   ///< reinit mode + inner loop: trials see the raw design, only accepted iterates
                                            ///< are redistanced (one extra evaluation per outer iteration)
        bool mDebugInnerTrials   = false;   ///< log each trial's approximation vs true criteria and the criteria at the
                                            ///< other design (raw psi vs R(psi)), by one extra diagnostic evaluation
        int  mNumLocked    = 0;      ///< variables whose asymptote spread is at the floor, from the last asymp()
                                                  ///< (kept separate so xmma stays the raw subproblem
                                                  ///< solution for raaupdate/conservatism); committed to xval

        // ---- NaN guard ----
        real mMoveScale = 1.0;      ///< move-limit multiplier; shrunk while backtracking off a bad step

        // ---- oscillation-adaptive per-ADV move limit (move_adapt) ----
        bool mMoveAdapt       = false;    ///< off => every ADV uses mStepSize (the original fixed box)
        real mMoveShrink      = 0.7;      ///< factor applied when the last two increments flip sign
        real mMoveGrow        = 1.2;      ///< factor applied when they agree
        real mMoveMin         = 1.0e-4;   ///< floor (fraction of range)
        real mMoveMax         = -1.0;     ///< cap (fraction of range); < 0 => mStepSize
        real mMoveInit        = -1.0;     ///< initial value; < 0 => mStepSize
        bool mMoveGate        = true;     ///< adapt only ADVs whose step hit the subproblem move bound
        std::vector< double > mMove;      ///< per-ADV move limit (fraction of range)
        std::vector< int >    mMoveHit;   ///< 1 if the accepted step of this ADV sat on its move bound
        std::vector< int >    mMoveHitPrev;   ///< same for the previous accepted step

        // ---- rise guard (move_guard_*) ----
        real mGuardRise   = -1.0;     ///< relative objective rise that rejects a step from a feasible design; < 0 off
        real mGuardFeas   = 2.0e-3;   ///< feasibility tolerance of the current design
        real mGuardShrink = 0.5;      ///< move-box factor per rejected trial
        sint mGuardMaxBt  = 4;        ///< backtracks before the last trial is accepted

        // ---- objective scaling (degenerate-evaluation safeguard) ----
        // GCMMA assumes an O(1-100) objective (Svanberg's own scaling recommendation). Degenerate
        // XFEM cuts can hand back gradients of 1e10+, which overflow the subproblem's interior-point
        // arithmetic into NaN iterates. When max|df0dx| exceeds the target, the objective and its
        // gradient are scaled UNIFORMLY for the current outer iteration -- this preserves the descent
        // direction exactly (unlike per-entry clipping, which distorts it) and returns the subproblem
        // to its proven operating range. The stored objective value is re-booked at each grad() call
        // so all quantities within an outer iteration share one consistent scale.
        real mObjScale             = 1.0;      ///< current uniform objective scale (1 = unscaled)
        static constexpr real sObjGradEngage = 1.0e9;    ///< scale only above this -- the solver demonstrably
                                                         ///< handles 1e8 raw and overflows ~1e10; engaging in
                                                         ///< proven regimes needlessly perturbs the trajectory
        static constexpr real sObjGradTarget = 1.0e8;    ///< decade-quantized scale brings max into (1e7, 1e8]

        real mKKTRef = -1.0;    ///< KKT reference residual (first kktcheck of the run; < 0 = not set yet)

      public:
        /**
         * Optimizer state carried across a restart (restart_keep_state). Widths and move steps are in ADV units, so
         * that they survive a change of the bounds; the centre is the design they belong to.
         */
        struct Restart_State
        {
            bool                  mValid = false;
            std::vector< double > mCentre;         ///< design the asymptote widths are centred on
            std::vector< double > mLowWidth;       ///< centre - low of the next asymptotes (> 0)
            std::vector< double > mUppWidth;       ///< upp - centre of the next asymptotes (> 0)
            std::vector< double > mMoveStep;       ///< per-ADV move limit times the ADV range
            std::vector< double > mLowerBounds;    ///< bounds the state was saved with (diagnostic)
            std::vector< double > mUpperBounds;
            std::vector< double > mRaa;        ///< constraint conservatism weights raa_j
            std::vector< double > mPenalty;    ///< constraint penalties c_j
            real                  mRaa0     = 0.0;
            real                  mObjScale = 1.0;
            real                  mKKTRef   = -1.0;
            Matrix< IdMat >       mIjklIds;          ///< ADV IDs in the order of the arrays (empty: the ADVs carry none)
            sint                  mNumProcs = -1;    ///< processor count the order belongs to (-1: unknown, old files)
        };

        /** what apply_restart_state() did per ADV */
        struct Restart_Counts
        {
            int mKept  = 0;    ///< state kept (re-centred)
            int mMoved = 0;    ///< mapped, but moved by more than the moved fraction of the range: fresh
            int mFresh = 0;    ///< mapped to -1 (row map): fresh
        };

        /**
         * Re-centre a kept state on the restart design aX (bounds aXMin, aXMax): for ADV i with old index
         * k = aRowMap[ i ] (identity if aRowMap is empty), low_i = x_i - w_low_k, upp_i = x_i + w_upp_k and
         * move_i = step_k / range_i, clamped to [aMoveMin, aMoveMax]. ADV i starts fresh (low/upp = x_i -+ aAsympAdapt0
         * * range_i, move_i = aFreshMove) when k = -1, or when aMovedFraction > 0 and |x_i - centre_k| > aMovedFraction
         * * range_i (the restart moved it out of the region the kept state describes). Widths are clamped to the
         * asymptote limits [0.01, 10] * range_i of asymp(). aMove is left unchanged if aKeepMove is false (fixed box).
         */
        static Restart_Counts apply_restart_state(
                const Restart_State&         aState,
                const std::vector< sint >&   aRowMap,
                const std::vector< double >& aX,
                const std::vector< double >& aXMin,
                const std::vector< double >& aXMax,
                real                         aMovedFraction,
                real                         aAsympAdapt0,
                bool                         aKeepMove,
                real                         aFreshMove,
                real                         aMoveMin,
                real                         aMoveMax,
                std::vector< double >&       aLow,
                std::vector< double >&       aUpp,
                std::vector< double >&       aMove );

        /**
         * Map from the current ADV order to the order of a state read from a file: entry i is the position in aState
         * of the ADV at position i now. Empty if the two orders are the same.
         *
         * The ADV order depends on the processor count, so the state identifies its ADVs like the ADV restart file
         * does, by their IDs (IjklIds). With IDs on both sides the state is mapped by ID (an ADV without an ID, gNoID,
         * keeps its position). Without them (a state file written before the IDs were stored, or ADVs that carry
         * none), the order can only be the same if the processor count is: that of the state, or for an old file that
         * of the ADV restart file read with it (aDesignFileNumProcs, -1 if none was read). Anything else is an error.
         *
         * @param aState state as read from the file
         * @param aIjklIds current ADV IDs (any other size than aNumVar: the ADVs carry no IDs)
         * @param aNumVar current number of ADVs (must be that of the state)
         * @param aNumProcs current processor count
         * @param aDesignFileNumProcs processor count of the ADV restart file of this run, -1 if there is none
         * @param aFileName state file name, for the messages
         */
        static std::vector< sint > restart_state_id_map(
                const Restart_State&   aState,
                const Matrix< IdMat >& aIjklIds,
                uint                   aNumVar,
                sint                   aNumProcs,
                sint                   aDesignFileNumProcs,
                const std::string&     aFileName );

        /**
         * write aState to an HDF5 file (datasets Centre, LowWidth, UppWidth, MoveStep, bounds, Raa, Penalty, scalars,
         * NumProcs, and IjklIds if the state has them)
         */
        static void write_restart_state( const std::string& aFileName, const Restart_State& aState, uint aOptIter );

        /** read a state written by write_restart_state(); a file without IjklIds / NumProcs gives none / -1 */
        static Restart_State read_restart_state( const std::string& aFileName );

        /** file name of the state saved with the committed design of optimization iteration aOptIter */
        static std::string restart_state_file_name( uint aAlgIndex, uint aOptIter );

        /** the kept state: of the last completed iteration, or the one applied at the start of the current solve */
        const Restart_State&
        get_kept_state() const
        {
            return mKeptState;
        }

        /** the state applied at the start of the last warm solve (re-centred), and what was kept of it */
        const Restart_State&
        get_warm_start_state() const
        {
            return mWarmStartState;
        }

        const Restart_Counts&
        get_warm_start_counts() const
        {
            return mWarmStartCounts;
        }

        /** number of solves that started from a kept state */
        int
        get_num_warm_solves() const
        {
            return mNumWarmSolves;
        }

      private:
        // ---- warm restart (restart_keep_state) ----
        bool                mKeepState     = false;    ///< keep the optimizer state across restarts
        real                mMovedFraction = 0.25;     ///< restart_keep_state_moved_fraction
        std::string         mStateFile;                ///< restart_state_file (read once, at the first solve)
        std::vector< sint > mStateRowMap;              ///< restart_state_row_map(_file)
        bool                mStateFileRead   = false;
        bool                mWarmSolve       = false;    ///< the current solve started from a kept state
        uint                mLastFuncOptIter = 0;        ///< optimization iteration of the last func() evaluation
        Restart_State       mKeptState;
        Restart_State       mWarmStartState;
        Restart_Counts      mWarmStartCounts;
        int                 mNumWarmSolves = 0;

        /** snapshot of the state at the end of a completed iteration (next asymptotes, move limits, raa, scale) */
        Restart_State capture_restart_state();

        /** record the order of the per-ADV arrays of aState: the problem's ADV IDs and the processor count */
        void set_restart_state_order( Restart_State& aState ) const;

        /** apply the kept (or file) state at the start of a solve; returns true if this solve is warm */
        bool begin_warm_solve();

        /** asymptote history update of asymp() for one variable (iteration > 2), into aLowOut / aUppOut */
        static void history_asymptote( double aX, double aXOld1, double aXOld2, double aLow, double aUpp, double aXMin, double aXMax, double aShrink, double aExpand, double& aLowOut, double& aUppOut );

        /** true if the objective and all constraints are finite */
        bool all_finite( real aF0, const std::vector< double >& aF ) const;

        /** true if the current trial criteria (f0valnew, fvalnew) are all finite */
        bool criteria_finite() const;

        /** true if the approximations f0app, fapp bound the trial criteria from above (GCMMA conservatism) */
        bool is_conservative() const;

        /** outcome of the NaN-guard backtrack on a trial design */
        enum class Trial_Status
        {
            finite,     ///< the trial criteria are finite (possibly after shrinking the move limit)
            restart,    ///< the workflow requested an optimization restart during an evaluation
            failed      ///< still non-finite after the maximum number of backtracks
        };

        /**
         * NaN guard for the trial point xmma: while its criteria are non-finite, halve the move limit, re-solve
         * the subproblem and re-evaluate. Used for the first trial AND every conservatism inner trial, so a failed
         * evaluation is never tested for conservatism, committed, or differentiated (MB-023).
         */
        Trial_Status backtrack_non_finite( int aIter, int& aNumFunc );

        /**
         * Re-evaluate the accepted design xval after a rejected iteration. Every failed trial evaluation has
         * overwritten the shared Problem (ADVs, geometry, analysis, criteria) that the Manager reads; this restores
         * a state consistent with the design the optimizer holds (MB-024). Returns false if the restore evaluation
         * requested a restart; a non-finite restore is a MORIS_ERROR (the evaluation is not reproducible).
         */
        bool restore_accepted_state( int aIter, int& aNumFunc );

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

        /** true if inner trials are evaluated at the raw design (reinit_in_place + reinit_on_accept_only + inner loop) */
        bool accept_only_redistance() const;

        /**
         * Evaluate the trial point xmma into f0valnew/fvalnew (redistanced unless accept_only_redistance()). With
         * debug_inner_trials, first evaluates the other design (R(psi) or psi) for the diagnostic log line.
         */
        void trial_func( int aIter, int aInner, int& aNumFunc );

        /** debug_inner_trials: one log line comparing the approximation with the true trial criteria */
        void log_trial( int aIter, int aInner, real aAltF0, const std::vector< double >& aAltF, real aDeltaR, bool aHaveAlt ) const;

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

        static double mx( double a, double b ) { return a > b ? a : b; }
        static double mn( double a, double b ) { return a < b ? a : b; }

      public:
        /**
         * Oscillation-adaptive move-limit update for one accepted outer iteration (move_adapt).
         * For each ADV i with increments d1 = aX - aXOld1, d0 = aXOld1 - aXOld2 and range r_i, if aHit[i] != 0
         * (the step sat on its move bound; pass all ones for the ungated rule):
         *   d1 * d0 < 0 (sign flip)  => aMove[i] = aShrink * min( aMove[i], |d1| / r_i )
         *   d1 * d0 > 0 (agreement)  => aMove[i] = aGrow * aMove[i]
         * otherwise unchanged; then clamped to [aMin, aMax].
         * @return number of ADVs shrunk
         */
        static int adapt_move_limits(
                std::vector< double >&       aMove,
                const std::vector< double >& aX,
                const std::vector< double >& aXOld1,
                const std::vector< double >& aXOld2,
                const std::vector< double >& aRange,
                const std::vector< int >&    aHit,
                real                         aShrink,
                real                         aGrow,
                real                         aMin,
                real                         aMax );

        /**
         * Rise-guard test: true if the current design is feasible (all aCon <= aFeas) and the trial objective
         * exceeds the current one by more than aRise * |aF0|.
         */
        static bool guard_rejects( real aF0, real aF0New, const std::vector< double >& aCon, real aRise, real aFeas );

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

        /**
         * Number of design variables whose asymptote spread (upp - low) / range was below 0.05 in the
         * last asymp() call ("locked": their move limit is choked to ~1% of the range). Diagnostic.
         */
        int
        get_num_locked_asymptotes() const
        {
            return mNumLocked;
        }

        /**
         * Reinit mode only: whether to shift the history and asymptotes by the redistance (default true, the
         * warm restart). false reproduces the pre-2026-10-10 behaviour; used as the unit test's negative control.
         */
        void
        set_shift_history( bool aShift )
        {
            mShiftHistory = aShift;
        }
    };
}    // namespace moris::opt

#endif /* MORIS_CL_OPT_ALGORITHM_MMA_HPP_ */
