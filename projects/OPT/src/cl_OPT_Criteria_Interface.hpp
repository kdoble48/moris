/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Criteria_Interface.hpp
 *
 */

#pragma once

#include "cl_Matrix.hpp"
#include "cl_Parameter_List.hpp"

namespace moris::opt
{
    class Criteria_Interface
    {
      private:
        bool             mEvaluated = false;
        Matrix< DDRMat > mSensitivities;

      protected:
        bool mInitializeOptimizationRestart = false;

        /// set by the optimizer for an evaluation whose design must not be redistanced in place (an inner GCMMA trial
        /// in reinit_on_accept_only mode); read on rank 0 by the workflow and broadcast there
        bool mReinitSuppressed = false;

        uint mReinitializeIterInterval = MORIS_UINT_MAX;
        uint mReinitializeIterFirst    = MORIS_UINT_MAX;

        uint mIter = 0;

      public:
        /**
         * Constructor
         */
        Criteria_Interface()
        {
        }

        /**
         * Destructor
         */
        virtual ~Criteria_Interface()
        {
        }

        /**
         * Initializes the vectors of ADV values, lower bounds, and upper bounds
         *
         * @param aADVs Initial ADVs to be filled.
         * @param aLowerBounds Lower ADV bounds to be filled.
         * @param aUpperBounds Upper ADV bounds to be filled.
         */
        virtual void initialize(
                Vector< real >& aADVs,
                Vector< real >& aLowerBounds,
                Vector< real >& aUpperBounds,
                Matrix< IdMat >&  aIjklIds ) = 0;

        /**
         * Gets the criteria values given a new set of ADVs
         *
         * @return vector of criteria
         */
        Vector< real > get_criteria( Vector< real >& aNewADVs );

        /**
         * Gets the derivatives of the criteria with respect to the advs, and computes if not already done
         *
         * @return matrix d(criteria)_i/d(adv)_j
         */
        Matrix< DDRMat > get_dcriteria_dadv();

        /**
         * indicates a request to restart the optimization
         *
         */
        bool
        get_restart_optimization()
        {
            return mInitializeOptimizationRestart;
        }

        /**
         * requests a restart of the optimization (e.g. when an optimization algorithm
         * produces a non-finite design and needs to be re-initialized from the current
         * design); the flag is cleared by the workflow on the restart's initialize
         */
        void
        request_restart_optimization()
        {
            mInitializeOptimizationRestart = true;
        }

        /**
         * Suppress (or allow) the in-place SDF redistance for the following criteria evaluations. Used by the native
         * MMA in reinit_on_accept_only mode: inner conservatism trials are evaluated at the raw subproblem solution,
         * only the accepted iterate is redistanced. Virtual so that an interface manager can forward it.
         */
        virtual void
        set_reinit_suppressed( bool aSuppressed )
        {
            mReinitSuppressed = aSuppressed;
        }

        bool
        is_reinit_suppressed() const
        {
            return mReinitSuppressed;
        }

        void
        set_reinitialize_iter( uint aReinitializeIter )
        {
            mReinitializeIterInterval = aReinitializeIter;
        }

        void
        set_first_reinitialize_iter( uint aReinitializeIter )
        {
            mReinitializeIterFirst = aReinitializeIter;
        }

      private:
        /**
         * Gets the criteria values given a new set of ADVs
         *
         * @return vector of criteria
         */
        virtual Vector< real > perform( Vector< real >& aNewADVs ) = 0;

        /**
         * Computes the derivatives of the criteria with respect to the advs
         *
         * @return matrix d(criteria)_i/d(adv)_j
         */
        virtual Matrix< DDRMat > compute_dcriteria_dadv() = 0;
    };
}
