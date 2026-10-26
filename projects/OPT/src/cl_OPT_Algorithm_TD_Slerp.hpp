/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_TD_Slerp.hpp
 *
 */

#pragma once

#include <string>
#include <vector>

#include "cl_OPT_Algorithm.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_OPT_Problem.hpp"

namespace moris::opt
{
    /**
     * @brief Topological-derivative level-set evolution on the unit sphere ("td_slerp"; Amstutz and Andrae 2006,
     * Teixeira, Krenn, Gangl and Marussig 2025, Algorithm 1).
     *
     * The design phi (the level-set B-spline coefficients, every ADV) is kept on the unit sphere of the inner product
     * <a, b>_M = a^T M b (M the design mass matrix of the criteria interface, the L2(D) inner product of the level-set
     * functions). With g~ the filtered topological derivative of the last evaluation (criteria interface
     * get_design_direction, MORISGENERAL TD_OPERATOR) and theta = acos( <phi, g~>_M / ( |phi|_M |g~|_M ) ):
     *
     *   phi_t = [ sin( (1 - k) theta ) phi / |phi|_M + sin( k theta ) g~ / |g~|_M ] / sin( theta ),   phi_t /= |phi_t|_M
     *
     * with k = kappa_init, halved until the trial objective satisfies J_t <= J_i ( 1 + accept_rel_tol ) (non-strict:
     * TKGM25 accept a trial that leaves J unchanged, e.g. their iteration 1 at k = 1/4). Exits, checked before each
     * update and logged: THETA (theta < theta_tol, converged), THETA_PI (theta > 180 - theta_pi_tol: g~ ~ -s phi, check
     * the sign convention), N_MAX (max_its accepted iterations), LINE_SEARCH (k < kappa_min: phi_i is restored and
     * written as the last restart file), G_ZERO (|g~|_M = 0).
     *
     * One evaluation = compute_design_criteria; the restart file ADV_Alg_<alg>_Iter_<n>.hdf5 is written before every
     * evaluation (index captured before the evaluation advances it, MB-016, MB-021). g~ is fetched only right after an
     * accepted evaluation (the interface holds the last evaluation, a rejected trial's g~ is never used). No gradients
     * and no constraints are used. There is no bound projection: a coefficient at or beyond a bound is a fatal error
     * (raise the discretization bounds), and an evaluation that changes the design (in-place reinitialization) is an
     * error, since it would leave the sphere.
     *
     * The start (the problem's ADVs, the restart ADVs on a restart) is evaluated as given, since M exists only after an
     * evaluation; with normalize_start it is then divided by |phi_0|_M and evaluated again, unless |phi_0|_M = 1 to 1e-12.
     *
     * Every ADV must be a coefficient of the TD_OPERATOR design field: the update rescales the whole ADV vector.
     * restart_index offsets the accepted-iteration count (max_its counts in total); the evaluation (restart file)
     * numbering of the logger is not changed.
     */
    class Algorithm_TD_Slerp : public Algorithm
    {
      public:
        /// exit reasons
        enum class Exit
        {
            NONE,
            THETA,
            THETA_PI,
            N_MAX,
            LINE_SEARCH,
            G_ZERO
        };

        /// name of an exit reason
        static std::string get_exit_name( Exit aExit );

      private:
        real mThetaTol;          // exit THETA below this angle (degrees)
        real mThetaPiTol;        // exit THETA_PI above 180 - this angle (degrees)
        real mKappaInit;         // first trial step fraction of an iteration
        real mKappaMin;          // line search fails below this step fraction
        real mAcceptRelTol;      // accept J_t <= J_i ( 1 + tol )
        bool mKappaCarry;        // kappa_policy "carry": start from the last accepted kappa
        bool mNormalizeStart;    // divide the start by its M norm

        // results of the last solve (rank 0)
        Exit                mExit           = Exit::NONE;
        uint                mNumAccepted    = 0;       // accepted-iteration count at the exit (incl. restart_index)
        uint                mNumEvaluations = 0;       // criteria evaluations of this solve
        std::vector< real > mObjectiveHistory;         // J_i per accepted iterate (index 0: the start)
        std::vector< real > mThetaHistory;             // theta_i (degrees, M inner product) per accepted iterate
        std::vector< real > mThetaEuclideanHistory;    // theta_i with the Euclidean inner product (diagnostic)
        std::vector< real > mKappaHistory;             // accepted kappa per accepted iteration

        // one evaluation: restart file, criteria, objective
        real evaluate( Vector< real >& aADVs );

        // <a, b>_M and |a|_M through the problem's criteria interface
        real mass_inner( const Vector< real >& aFirst, const Vector< real >& aSecond );
        real mass_norm( const Vector< real >& aVector );

        // errors if a coefficient is not finite or at/beyond a bound
        void check_bounds( const Vector< real >& aADVs );

        // the iteration on rank 0
        void td_slerp_solve();

      public:
        explicit Algorithm_TD_Slerp( const Parameter_List& aParameterList );

        ~Algorithm_TD_Slerp() override;

        uint solve( uint aCurrentOptAlgInd, std::shared_ptr< Problem > aOptProb ) override;

        //------------------------------------------------------------------------------
        // results of the last solve (rank 0)

        Exit get_exit() const { return mExit; }
        uint get_num_accepted() const { return mNumAccepted; }
        uint get_num_evaluations() const { return mNumEvaluations; }

        const std::vector< real >& get_objective_history() const { return mObjectiveHistory; }
        const std::vector< real >& get_theta_history() const { return mThetaHistory; }
        const std::vector< real >& get_theta_euclidean_history() const { return mThetaEuclideanHistory; }
        const std::vector< real >& get_kappa_history() const { return mKappaHistory; }
    };
}    // namespace moris::opt
