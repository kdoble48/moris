/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_Adam_Projected.hpp
 *
 */

#ifndef MORIS_CL_OPT_ALGORITHM_ADAM_PROJECTED_HPP_
#define MORIS_CL_OPT_ALGORITHM_ADAM_PROJECTED_HPP_

#include <string>
#include <vector>

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_OPT_Problem.hpp"

namespace moris::opt
{
    /**
     * @brief Adam with an exact one-constraint projection ("adam_proj"; 3DGS-style first-order method).
     *
     * For problems min f(x) s.t. g(x) <= 0 (or = 0) with ONE constraint, x within bounds. Each iteration:
     *   1. Adam step on the objective gradient (optionally its component tangent to the constraint), with a
     *      per-ADV learning rate lr_k * pattern_i, lr_k decaying exponentially from learning_rate to
     *      learning_rate * lr_decay_final over max_its iterations; clipped to the bounds.
     *   2. Restore the constraint: x <- x + t d along one fixed direction d (projection_direction "uniform": the
     *      projection_mask pattern itself, e.g. every alpha ADV = a global opacity shift; "gradient": the masked
     *      constraint gradient), t from a 1D root find on the TRUE constraint (linearized predictor from the
     *      current gradient, then secant corrections, each a forward solve) until |g| <= projection_tol.
     *   3. Gradients at the projected point; every iterate is accepted (restart file written).
     * There is no constraint in the step itself, so no active-constraint bounce; the decaying rate makes the tail
     * converge by construction.
     */
    class Algorithm_Adam_Projected : public Algorithm
    {
      private:
        real                mLearningRate;          // initial base learning rate (ADV units)
        std::vector< real > mRatePattern;           // per-ADV multipliers, cycled over the ADVs
        real                mDecayFinal;            // lr at the last iteration / initial lr
        real                mBeta1;
        real                mBeta2;
        real                mEpsilon;
        std::vector< real > mMaskPattern;           // projection mask, cycled over the ADVs (empty: all ones)
        std::string         mDirection;             // "uniform" | "gradient"
        std::string         mMode;                  // "equality" | "inequality"
        real                mProjectionTol;
        sint                mProjectionMaxIts;
        bool                mTangentMoments;        // feed Adam the objective gradient tangent to the constraint
        sint                mNumForwardSolves = 0;

        static std::vector< real > parse_pattern( const std::string& aString );

        // forward solve at aX: returns objective and constraint (first constraint, 0 if none)
        void evaluate( const std::vector< real >& aX, real& aF, real& aG );

      public:
        explicit Algorithm_Adam_Projected( const Parameter_List& aParameterList );

        ~Algorithm_Adam_Projected() override;

        uint solve( uint aCurrentOptAlgInd, std::shared_ptr< Problem > aOptProb ) override;

        void adam_solve();
    };
}    // namespace moris::opt

#endif /* MORIS_CL_OPT_ALGORITHM_ADAM_PROJECTED_HPP_ */
