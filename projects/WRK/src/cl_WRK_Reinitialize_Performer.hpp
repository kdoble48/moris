/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_Reinitialize_Performer.hpp
 *
 */

#ifndef SRC_cl_WRK_Reinitialize_Performer
#define SRC_cl_WRK_Reinitialize_Performer

#include <memory>
#include "cl_Parameter_List.hpp"
#include "cl_Matrix.hpp"
#include "cl_Library_IO.hpp"

namespace moris
{
    namespace hmr
    {
        class HMR;
        class Mesh;
    }    // namespace hmr

    namespace mtk
    {
        class Field;
        class Mesh_Manager;
        class Mesh;
    }    // namespace mtk

    namespace MSI
    {
        enum class Dof_Type;
    }

    namespace gen
    {
        class Geometry_Engine;
    }

    namespace mdl
    {
        class Model;
    }

    namespace xtk
    {
        class Cut_Integration_Mesh;
    }

    namespace wrk
    {
        /**
         * @brief Mode for SDF reinitialization
         */
        enum class ReinitMode
        {
            FIELD_REMAP,     ///< Original: remap solution field to target mesh
            INTERFACE_SDF    ///< New: compute exact SDF from XTK interface facets
        };

        class Reinitialize_Performer
        {
          private:
            std::shared_ptr< Library_IO > mLibrary = nullptr;

            // Solution Dof types
            Vector< moris::MSI::Dof_Type > mDofTypes;

            // adv filed name
            std::string mADVFiledName;

            // discretization mesh index for adof
            moris_index mAdofMeshIndex;

            // coefficents that will be populated with target and source field coeff
            moris::Matrix< moris::DDRMat > mCoefficients;

            // reinitialization frequency
            sint mReinitializationFrequency;

            // mtk fields that will be used to overwrite gen fields
            Vector< std::shared_ptr< mtk::Field > > mMTKFields;

            // mesh outputting params
            std::string mOutputMeshFile;
            moris::real mTimeOffset;

            // SDF reinitialization mode and parameters
            ReinitMode  mReinitMode     = ReinitMode::FIELD_REMAP;
            moris_index mMaterialPhase  = 0;    ///< Bulk phase index for "inside" (negative SDF)
            Vector< moris_index > mExcludePhases;  ///< Bulk phases to exclude from interface (e.g., external zone)
            uint mSmoothingIterations = 0;  ///< Number of Laplacian smoothing iterations (0 = disabled)
            real mSmoothingLambda = 0.5;  ///< Laplacian smoothing strength (0.0-1.0)
            bool mInitializeAtStart = false;  ///< Reinit once at the first iteration to start from a true SDF
            real mWeanTol = 0.0;              ///< Wean-off tolerance on the reinit's relative change (0 = disabled)
            bool mLastReinitBelowWeanTol = false;  ///< Set when the most recent reinit's change fell below mWeanTol
            bool mReinitEveryStep = false;    ///< Redistance every iteration in-place (no restart); smooth every reinit_freq

          public:
            //------------------------------------------------------------------------------

            /**
             * @brief Construct a new Reinitialize_Performer object
             *
             * @param aParameterlist
             * @param aLibrary
             */

            Reinitialize_Performer( const std::shared_ptr< Library_IO >& aLibrary );

            //------------------------------------------------------------------------------

            /**
             * @brief Destroy the Reinitialize_Performer object
             *
             */

            ~Reinitialize_Performer(){};

            //------------------------------------------------------------------------------

            /**
             * @brief
             *
             * @param aHMRPerformers
             * @param aGENPerformer
             * @param aMTKPerformer
             * @param mMDLPerformer
             */

            void perform(
                    Vector< std::shared_ptr< hmr::HMR > >&            aHMRPerformers,
                    Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
                    Vector< std::shared_ptr< mtk::Mesh_Manager > >&   aMTKPerformer,
                    Vector< std::shared_ptr< mdl::Model > >           mMDLPerformer );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the reinitialization frequency object
             *
             * @return moris::uint
             */

            moris::sint
            get_reinitialization_frequency() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the coefficients object
             *
             * @return Matrix <DDRMat> const&
             */

            Matrix< DDRMat > const &
            get_coefficients() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Get the name of the ADV (level-set) field this performer reinitializes.
             *        Used to target that design's ADV slice when writing reinitialized coefficients.
             *
             * @return std::string const&
             */
            std::string const &
            get_adv_field_name() const { return mADVFiledName; }

            //------------------------------------------------------------------------------

            /**
             * @brief clip values of advs
             *
             * @param aGENPerformer
             */

            void
            impose_upper_lower_bound( Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer, mtk::Field* aField );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the mtk fields object get the replaced fields to use in GEN
             *
             * @return Vector< std::shared_ptr< mtk::Field > >
             */
            Vector< std::shared_ptr< mtk::Field > >
            get_mtk_fields() const;

            //------------------------------------------------------------------------------

            /**
             * @brief  output the mapped and original field on the same mesh
             *
             * @param aTarget
             * @param aSource
             */
            void
            output_fields( mtk::Field* aTarget, mtk::Field* aSource, std::string aExoFileName ) const;

            //------------------------------------------------------------------------------

            /**
             * @brief Set the reinitialization mode
             *
             * @param aMode FIELD_REMAP (original) or INTERFACE_SDF (new)
             */
            void
            set_reinit_mode( ReinitMode aMode );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the current reinitialization mode
             */
            ReinitMode
            get_reinit_mode() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Whether to initialize the design as an SDF on the first iteration
             *
             * @return true if a one-time reinit should fire at the first optimization iteration
             */
            bool
            get_initialize_at_start() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Whether the most recent reinit's relative change fell below the wean tolerance
             *
             * @return true if the reinit has stopped doing meaningful work and should be weaned off
             */
            bool
            last_reinit_below_wean_tol() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Wean-off tolerance on the relative design-variable change (0 = disabled)
             *
             * @return the configured sdf_reinit_wean_tol
             */
            real
            get_wean_tol() const { return mWeanTol; }

            //------------------------------------------------------------------------------

            /**
             * @brief Whether to redistance the SDF every iteration in-place (no GCMMA restart)
             *
             * @return true if every-step redistancing is enabled
             */
            bool
            get_reinit_every_step() const { return mReinitEveryStep; }

            //------------------------------------------------------------------------------

            /**
             * @brief Set the material phase index (for sign determination)
             *
             * @param aPhase Bulk phase index that represents "inside" (negative SDF)
             */
            void
            set_material_phase( moris_index aPhase );

            //------------------------------------------------------------------------------

            /**
             * @brief Get the excluded phases for SDF interface computation
             *
             * @return Vector of bulk phase indices to exclude from interface facets
             */
            Vector< moris_index > const&
            get_exclude_phases() const;

            //------------------------------------------------------------------------------

            /**
             * @brief Set phases to exclude from SDF interface computation
             *
             * @param aExcludePhases Vector of bulk phase indices to exclude
             */
            void
            set_exclude_phases( Vector< moris_index > const& aExcludePhases );

            //------------------------------------------------------------------------------

            /**
             * @brief Get number of Laplacian smoothing iterations for SDF reinit
             */
            uint
            get_smoothing_iterations() const { return mSmoothingIterations; }

            /**
             * @brief Set number of Laplacian smoothing iterations for SDF reinit
             * @param aNumIter Number of iterations (0 = disabled, 1-5 recommended)
             */
            void
            set_smoothing_iterations( uint aNumIter ) { mSmoothingIterations = aNumIter; }

            //------------------------------------------------------------------------------

            /**
             * @brief Compute SDF from XTK interface facets
             *
             * Uses the triangulated interface from XTK to compute exact signed
             * distances at all mesh nodes. This is an alternative to field remapping
             * that maintains the exact SDF property (|∇φ| = 1).
             *
             * @param aCutMesh XTK cut integration mesh containing interface facets
             * @param aInterpMesh Interpolation mesh (Lagrange nodes)
             * @param aGENPerformer Geometry engine for ADV updates
             * @param aMTKPerformer Mesh manager
             * @param aSDFField Output: MTK field with reinitialized SDF values
             * @param aApplySmoothing Whether to apply Laplacian smoothing this call (decoupled from
             *        redistancing in every-step mode, where smoothing runs only every reinit_freq)
             */
            void
            compute_sdf_from_interface(
                    xtk::Cut_Integration_Mesh&                       aCutMesh,
                    mtk::Mesh*                                       aInterpMesh,
                    Vector< std::shared_ptr< gen::Geometry_Engine > >& aGENPerformer,
                    Vector< std::shared_ptr< mtk::Mesh_Manager > >&  aMTKPerformer,
                    std::shared_ptr< mtk::Field >&                   aSDFField,
                    bool                                             aApplySmoothing = true );

            //------------------------------------------------------------------------------

            /**
             * @brief Apply Laplacian smoothing to the ADV field
             *
             * This method smooths the ADV field by averaging each node's value
             * with its neighbors. Used for failure recovery to remove thin features
             * and small islands that cause numerical instabilities.
             *
             * @param aCoefficients ADV values to smooth (modified in place)
             * @param aInterpMesh   Interpolation mesh (for connectivity)
             * @param aNumIterations Number of smoothing iterations (default: 3)
             * @param aLambda       Smoothing factor in [0,1], higher = more smoothing (default: 0.5)
             */
            static void
            apply_laplacian_smoothing(
                    Matrix< DDRMat >&   aCoefficients,
                    mtk::Mesh*          aInterpMesh,
                    uint                aNumIterations = 3,
                    real                aLambda = 0.5 );
        };
    }    // namespace wrk
}    // namespace moris

#endif /* cl_WRK_Reinitailize_Performer.hpp */
