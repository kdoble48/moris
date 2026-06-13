/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Algorithm_GCMMA.hpp
 *
 */

#ifndef MORIS_CL_OPT_ALGORITHM_GCMMA_HPP_
#define MORIS_CL_OPT_ALGORITHM_GCMMA_HPP_

#include "core.hpp"
#include "cl_OPT_Algorithm.hpp"
#include "cl_Parameter_List.hpp"
#include "cl_OPT_Problem.hpp"

using namespace moris;

class OptAlgGCMMA : public moris::opt::Algorithm
{

private:
    bool mPrint;

    moris::uint mResFlag = 0;         // flag from GCMMA describing result of optimization algorithm
    moris::sint mMaxInnerIterations;  // maximum inner iterations per every optimization iteration
    moris::real mNormDrop;            // relative change in objective convergence criteria
    moris::real mAsympAdapt0;         // initial asymptote adaptation factor
    moris::real mAsympShrink;         // shrinking asymptote adaptation factor
    moris::real mAsympExpand;         // expanding asymptote adaptation factor
    moris::real mStepSize;            // GCMMA step size
    moris::real mPenalty;             // GCMMA constraint penalty
    moris::sint mIvers;               // GCMMA version

    // Uniform objective scaling (degenerate-evaluation safeguard): the TPL solver's internal
    // arithmetic overflows to NaN when degenerate XFEM cuts push the objective gradient past
    // ~1e10, while gradients up to ~1e8 are demonstrably handled raw. When max|dObj/dADV|
    // exceeds the engage threshold, objective and gradient are scaled uniformly (preserves
    // the descent direction exactly) by a DECADE-QUANTIZED factor that brings the max into
    // [target/10, target]. Quantization matters: the TPL's stored objective value cannot be
    // re-booked when the scale changes, so every change costs one inconsistent
    // conservativeness check -- a continuously re-derived scale corrupts every iteration,
    // while a quantized one only transitions when the gradient magnitude crosses a decade.
    moris::real                  mObjScale = 1.0;
    static constexpr moris::real sObjGradEngage = 1.0e9;    ///< scale only above this (true pathology)
    static constexpr moris::real sObjGradTarget = 1.0e8;    ///< scaled max lands in (1e7, 1e8]

    /**
     * @brief External function call for computing objective and constraints, to
     *        interface with gcmma library
     */
    friend void opt_alg_gcmma_func_wrap(
            OptAlgGCMMA* aOptAlgGCMMA,
            int&         aIter,
            double*      aAdv,
            double&      aObjval,
            double*      aConval );

    /**
     * @brief External function call for computing sensitivities, to interface
     *        with gcmma library
     */
    friend void opt_alg_gcmma_grad_wrap(
            OptAlgGCMMA* aOptAlgGCMMA,
            double*      aAdv,
            double*      aD_Obj,
            double**     aD_Con,
            int*         aActive );

public:

    /**
     * Constructor
     */
  OptAlgGCMMA( const moris::Parameter_List& aParameterList );

  /**
   * Destructor
   */
  ~OptAlgGCMMA() override;

  /**
   * @brief MORIS interface for solving of optimization problem using
   *        GCMMA
   *
   * @param[in] aCurrentOptAlgInd index of optimization algorithm
   * @param[in] aOptProb Object of type Problem containing relevant
   *            data regarding ADVs, the objective and constraints
   */
  uint solve(
          uint                                   aCurrentOptAlgInd,
          std::shared_ptr< moris::opt::Problem > aOptProb ) override;

  /**
   *@brief Run GCMMA algorithm
   */
  void gcmma_solve();

  /**
   *@brief Prints result of the GCMMA algorithm based on mStopFlag
   */
  void printresult();
};

#endif /* MORIS_CL_OPT_ALGORITHM_GCMMA_HPP_ */

