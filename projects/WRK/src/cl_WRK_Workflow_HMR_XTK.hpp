/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Workflow_HMR_XTK.hpp
 *
 */

#pragma once

#include "cl_WRK_Workflow.hpp"
#include "cl_WRK_Reinitialize_Performer.hpp"
#include <memory>
#include "cl_OPT_Criteria_Interface.hpp"
#include "moris_typedefs.hpp"    //MRS/COR/src
#include "cl_Vector.hpp"

#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_MTK_Enums.hpp"
#include "cl_Communication_Tools.hpp"

namespace xtk
{
    class Model;
}

namespace moris
{
    class Library_IO;
    //------------------------------------------------------------------------------
    namespace hmr
    {
        class HMR;
    }
    namespace mtk
    {
        class Mesh_Manager;
    }
    namespace ge
    {
        class GEN_Geometry_Engine;
    }
    namespace mdl
    {
        class Model;
    }
    namespace opt
    {
        class Manager;
    }

    namespace wrk
    {
        class Performer_Manager;
        //------------------------------------------------------------------------------
        // Naming convention here means the background mesh is constructed by HMR and the integration mesh is constructed by XTK
        class Workflow_HMR_XTK : public Workflow
        {
          private:
            /// Flag to track that SDF reinit has triggered NANs and coefficients are pending
            bool mSDFReinitPending = false;
            
            /// Flag to skip SDF reinit on the first iteration after restart
            bool mSkipNextSDFReinit = false;

            /// Flag to track that the one-time start-of-run SDF initialization has fired
            bool mSDFInitialized = false;

            /// Flag to permanently disable SDF reinit once it has weaned off (design converged)
            bool mSDFWeanedOff = false;

            /// Whether the MDL fields (e.g. an error indicator) were computed on the current mesh: set by a forward
            /// solve, cleared by a remesh. A remesh requested before any solve on the new mesh (XTK failed) must not
            /// refine on them (MB-033).
            bool mMDLFieldsOnCurrentMesh = false;

            /// Design variables (ADVs) captured at the previous scheduled reinit, used to measure
            /// the design-variable change between reinits for the wean-off decision
            Vector< real > mPrevReinitADVs;

            /// in-place taper over accepted iterates (created from the reinit performer's parameters on first use)
            std::unique_ptr< SDF_Wean_Monitor > mWeanMonitor;

          public:
            //------------------------------------------------------------------------------
            /**
             * constructor
             */
            Workflow_HMR_XTK( wrk::Performer_Manager* aPerformerManager );

            //------------------------------------------------------------------------------
            /**
             * destructor
             */
            ~Workflow_HMR_XTK() override{};

            //------------------------------------------------------------------------------
            /**
             * Initializes the vectors of ADV values, lower bounds, and upper bounds
             */
            void initialize(
                    Vector< real >& aADVs,
                    Vector< real >& aLowerBounds,
                    Vector< real >& aUpperBounds,
                    Matrix< IdMat >& aIjklIDs ) override;

            //------------------------------------------------------------------------------
            /**
             * Gets the criteria values given a new set of ADVs
             *
             * @return vector of criteria
             */
            Vector< real > perform( Vector< real >& aNewADVs ) override;

            //------------------------------------------------------------------------------
            /**
             * SDF taper (in-place every-step mode): feeds the relative step ||aAccepted - aPrevious|| / ||aPrevious|| of
             * an accepted optimizer iterate to the wean monitor; weans the redistance off for the rest of the run once it
             * fires. Rank 0 (the optimizer's rank); perform() broadcasts the result. An empty aAccepted marks an iterate
             * accepted at the inner-loop cap (not a convergence sample: resets the run of small changes).
             */
            void notify_accepted_iterate( const Vector< real >& aPrevious, const Vector< real >& aAccepted ) override;

            /** false once the in-place redistance has been tapered off */
            bool is_reinit_active() const override
            {
                return !mSDFWeanedOff;
            }

            //------------------------------------------------------------------------------
            /**
             * Gets the derivative of the criteria with respect to the advs
             *
             * @return matrix d(criteria)_i/d(adv)_j
             */
            Matrix< DDRMat > compute_dcriteria_dadv() override;

            //------------------------------------------------------------------------------
            /**
             * @brief Performs T-Matrix and MPC outputs (used for project work) if requested
             * by user.
             *
             * @param aMTKPerformer pointer to MTK Mesh_Manager
             * @param aXTKPerformer pointer to XTK_Model
             */
            bool
            output_T_matrices(
                    const std::shared_ptr< mtk::Mesh_Manager >& aMTKPerformer,
                    xtk::Model* const &                         aXTKPerformer );
        };
        //------------------------------------------------------------------------------
    }    // namespace wrk
} /* namespace moris */
