/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_ROL.hpp
 *
 */

#ifndef MORIS_CL_OPT_ALGORITHM_ROL_HPP_
#define MORIS_CL_OPT_ALGORITHM_ROL_HPP_

#include <string>
#include <vector>

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_OPT_Problem.hpp"

namespace moris::opt
{
    /**
     * @brief Trilinos ROL (Rapid Optimization Library) optimizer wrapper.
     *
     * Wraps ROL's trust-region / Lin-More method (the algorithm family Sandia's Plato uses),
     * with an optional Augmented-Lagrangian outer loop for general constraints. ROL is templated
     * C++ and consumes the MORIS Problem through the same value/gradient callbacks as the other
     * algorithms; see cl_OPT_Algorithm_ROL.cpp for the ROL::Objective / ROL::Constraint adapters.
     *
     * Header stays free of ROL/Teuchos types so it can be included without the Trilinos headers;
     * everything ROL-typed lives behind #ifdef MORIS_HAVE_ROL in the .cpp.
     */
    class Algorithm_ROL : public Algorithm
    {
      private:
        // --- curated settings mapped onto the ROL ParameterList (see build in rol_solve) ---
        std::string mStepType;                 // "Trust Region" (bound) / "Augmented Lagrangian" (constrained)
        std::string mSubproblemModel;          // trust-region subproblem model, e.g. "Lin-More"
        std::string mSubproblemSolver;         // trust-region subproblem solver, e.g. "Truncated CG"
        std::string mSecantType;               // Hessian approximation, e.g. "Limited-Memory BFGS"
        sint        mSecantStorage;            // secant history length
        real        mInitialTRRadius;          // initial trust-region radius
        real        mMaxTRRadius;              // maximum trust-region radius
        bool        mUseAugmentedLagrangian;   // wrap constraints in an AL outer loop
        sint        mSubproblemIterationLimit; // inner (bound-constrained TR) subproblem iteration cap in the AL loop
        real        mALInitialPenalty;         // AL initial penalty parameter
        real        mALPenaltyGrowth;          // AL penalty growth factor
        real        mGradientTol;              // status-test gradient tolerance
        real        mConstraintTol;            // status-test constraint tolerance
        real        mStepTol;                  // status-test step tolerance
        std::string mXmlFile;                  // optional full ROL ParameterList XML (overrides curated keys)

        // --- shared evaluation cache -------------------------------------------------------
        // ROL calls the objective and (M2) constraint adapters at the same design point; both a
        // forward criteria solve and an adjoint gradient solve are expensive, so we cache the ADV
        // vector each was last evaluated at and skip redundant solves. compute_design_criteria and
        // compute_design_criteria_gradients populate BOTH objective and constraint quantities in a
        // single solve, so one entry per stage suffices.
        std::vector< real > mCriteriaADVs;      // ADVs of the last forward (criteria) solve
        std::vector< real > mGradientADVs;      // ADVs of the last adjoint (gradient) solve
        bool                mHaveCriteria  = false;
        bool                mHaveGradients = false;

      public:
        /**
         * Constructor
         *
         * @param aParameterList algorithm parameter list (see prm::create_rol_parameter_list)
         */
        explicit Algorithm_ROL( const Parameter_List& aParameterList );

        /**
         * Destructor
         */
        ~Algorithm_ROL() override;

        /**
         * @brief MORIS interface for solving the optimization problem using ROL
         *
         * @param[in] aCurrentOptAlgInd index of optimization algorithm
         * @param[in] aOptProb          Problem holding ADVs, objective and constraints
         */
        uint solve(
                uint                       aCurrentOptAlgInd,
                std::shared_ptr< Problem > aOptProb ) override;

        /**
         * @brief Run the ROL solver (rank 0 only). Builds the ROL problem + parameter list and
         *        drives ROL::Solver; no-op MORIS_ERROR when compiled without ROL.
         */
        void rol_solve();

        /**
         * @brief Record an accepted design iterate (writes the restart file). Called by the ROL
         *        objective adapter's update() only on ROL::UpdateType::Accept, so restart files are
         *        written once per accepted step rather than once per (trial) forward evaluation.
         */
        void notify_accepted_step( const std::vector< real >& aADVs );

        /**
         * @brief Ensure design criteria (objective + constraints) are current at aADVs.
         *        Triggers a forward criteria solve only if aADVs differ from the cached point.
         *        Called by the ROL::Objective / ROL::Constraint adapters.
         */
        void ensure_criteria( const std::vector< real >& aADVs );

        /**
         * @brief Ensure design-criteria gradients are current at aADVs (implies ensure_criteria).
         *        Triggers an adjoint gradient solve only if aADVs differ from the cached point.
         */
        void ensure_gradients( const std::vector< real >& aADVs );
    };
}    // namespace moris::opt

#endif /* MORIS_CL_OPT_ALGORITHM_ROL_HPP_ */
