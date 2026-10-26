/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_WRK_TD_Design_Operator.hpp
 *
 */

#pragma once

#include <memory>
#include <string>

#include "moris_typedefs.hpp"
#include "cl_Vector.hpp"
#include "cl_Matrix.hpp"
#include "linalg_typedefs.hpp"
#include "cl_Parameter_List.hpp"

namespace moris
{
    namespace mtk
    {
        class Field;
        class Interpolation_Mesh;
    }    // namespace mtk

    namespace gen
    {
        class Geometry_Engine;
    }

    namespace wrk
    {
        //------------------------------------------------------------------------------

        /**
         * Topological-derivative design operator (MORISGENERAL sublist TD_OPERATOR).
         *
         * Maps an element-wise topological derivative g_e to the design B-spline space of a level-set field and filters it
         * (Galerkin route of the TD-driven level-set evolution, Teixeira et al. 2025 with Amstutz and Andrae 2006):
         *
         *   b_A = sum_e g_e int_e B_A                (load vector, exact in sum_A b_A = int_D g)
         *   ( gamma K + M ) g~ = b                     (Helmholtz filter, natural boundary condition)
         *
         * with M_AB = int_D B_A B_B and K_AB = int_D grad B_A . grad B_B assembled over the owned cells of the Lagrange
         * interpolation mesh that carries the design B-spline mesh, B_A = sum_n T_nA N_n (the vertex T-matrices of the
         * design discretization index), Gauss order p_L + 1 per direction. M also serves the inner product <a, b>_M of the
         * update. Only coefficients that appear in a T-matrix of the domain are unknowns; the others (zero T columns, e.g.
         * functions supported in the padding only) get g~ = 0.
         *
         * The values g_e are read per background cell from ELEMENTAL FEM fields with field_elemental_normalization =
         * ip_cell: every enriched interpolation cell row is added to its base cell, so the per-phase fields of a cut cell sum
         * to g_e = rho_e <g_in> + ( 1 - rho_e ) g_out.
         *
         * All vectors at the public interface are full ADV vectors (the level-set coefficients at their positions in the
         * ADV vector, every other entry ignored on input and zero on output). Serial only (asserted).
         */
        class TD_Design_Operator
        {
          private:
            // parameters
            std::string           mADVFieldName;
            Vector< std::string > mTDFieldNames;
            std::string           mCutRatioFieldName;
            real                  mVoidTDValue = 0.0;
            real                  mFilterGamma = 0.0;
            std::string           mOutputFile;

            // the mesh the operators were built on and its design discretization index
            mtk::Interpolation_Mesh* mMesh                = nullptr;
            uint                     mDiscretizationIndex = 0;
            uint                     mNumBuilds           = 0;

            // index maps: unknown -> position in the ADV vector, number of ADVs, cells of the mesh
            uint           mNumADVs  = 0;
            uint           mNumCells = 0;
            Vector< uint > mUnknownToADV;
            uint           mNumDroppedCoefficients = 0;

            // CSR (all rows of the symmetric matrices): M, K
            Vector< uint > mRowOffsets;
            Vector< uint > mColumns;
            Vector< real > mMassValues;
            Vector< real > mStiffnessValues;

            // per cell: unknowns and int_e B_A (CSR by cell index)
            Vector< uint > mCellOffsets;
            Vector< uint > mCellUnknowns;
            Vector< real > mCellIntegrals;
            Vector< real > mCellMeasures;

            // envelope Cholesky factor of gamma K + M in the reverse Cuthill-McKee order
            Vector< uint > mPermutation;        // new -> old
            Vector< uint > mEnvelopeFirst;      // first column of each row (new order)
            Vector< uint > mEnvelopeOffsets;    // start of each row in mFactor
            Vector< real > mFactor;             // L( i, first( i ) .. i )

            // diagnostics of the build
            real mDomainMeasure = 0.0;
            real mMassSum       = 0.0;
            real mMinDiagonal   = 0.0;
            real mMinPivot      = 0.0;

            // state of the last evaluation
            uint           mNumEvaluations = 0;
            Vector< real > mEvaluatedADVs;
            Vector< real > mCellValues;
            Vector< real > mLoad;
            Vector< real > mDirection;

            //------------------------------------------------------------------------------

            void factorize();

            Vector< real > to_unknowns( const Vector< real >& aADVVector ) const;

            Vector< real > to_adv_vector( const Vector< real >& aUnknownVector ) const;

            Vector< real > csr_product( const Vector< real >& aValues, const Vector< real >& aUnknownVector ) const;

            void write_output() const;

          public:
            //------------------------------------------------------------------------------

            /**
             * constructor
             * @param[ in ] aParameterList TD_OPERATOR sublist (prm::create_td_operator_parameter_list)
             */
            explicit TD_Design_Operator( const Parameter_List& aParameterList );

            //------------------------------------------------------------------------------

            ~TD_Design_Operator() = default;

            //------------------------------------------------------------------------------

            /**
             * builds M, K, the cell integrals and the factor of gamma K + M on a Lagrange interpolation mesh, and runs the
             * sanity checks (sum M = measure of the domain, row sums of M = int B_A, K 1 = 0, no zero rows, positive
             * Cholesky pivots); errors if one fails
             *
             * @param[ in ] aMesh                  Lagrange interpolation mesh with the T-matrices of the design index
             * @param[ in ] aDiscretizationIndex   design B-spline mesh index on aMesh
             * @param[ in ] aCoefficientADVPositions position in the ADV vector of every coefficient (local index), -1 if none
             * @param[ in ] aNumADVs               length of the ADV vector
             */
            void build_operators(
                    mtk::Interpolation_Mesh* aMesh,
                    uint                     aDiscretizationIndex,
                    const Vector< sint >&    aCoefficientADVPositions,
                    uint                     aNumADVs );

            /**
             * builds the operators for the GEN design named by adv_field (mesh, discretization index and ADV positions
             * from the geometry engine)
             */
            void build_operators( gen::Geometry_Engine& aGeometryEngine, uint aNumADVs );

            //------------------------------------------------------------------------------

            /**
             * per background cell (cell index of the mesh), the sum of the rows of the td fields over the enriched
             * interpolation cells of that cell, plus ( 1 - rho_e ) void_td_value if a cut ratio field is set
             */
            Vector< real > collect_cell_values( const Vector< std::shared_ptr< mtk::Field > >& aFields ) const;

            /** b_A = sum_e g_e int_e B_A, full ADV vector */
            Vector< real > assemble_load( const Vector< real >& aCellValues ) const;

            /** solves ( gamma K + M ) g~ = b, full ADV vectors */
            Vector< real > filter( const Vector< real >& aLoad ) const;

            /** a^T M b, full ADV vectors */
            real mass_product( const Vector< real >& aFirst, const Vector< real >& aSecond ) const;

            /** M a and K a, full ADV vectors */
            Vector< real > mass_times( const Vector< real >& aVector ) const;
            Vector< real > stiffness_times( const Vector< real >& aVector ) const;

            //------------------------------------------------------------------------------

            /**
             * one evaluation (called by Workflow_HMR_XTK::perform after the forward solve): builds the operators on the
             * first call and after a change of the mesh, then collects g_e from the FEM fields, assembles b, filters, and
             * caches the direction for the ADVs of this evaluation
             */
            void perform(
                    gen::Geometry_Engine&                          aGeometryEngine,
                    const Vector< std::shared_ptr< mtk::Field > >& aFields,
                    const Vector< real >&                          aADVs );

            //------------------------------------------------------------------------------
            // results of the last evaluation and diagnostics

            const Vector< real >& get_design_direction() const { return mDirection; }
            const Vector< real >& get_evaluated_advs() const { return mEvaluatedADVs; }
            const Vector< real >& get_cell_values() const { return mCellValues; }
            const Vector< real >& get_load() const { return mLoad; }
            uint                  get_num_evaluations() const { return mNumEvaluations; }
            uint                  get_num_builds() const { return mNumBuilds; }
            uint                  get_num_unknowns() const { return mUnknownToADV.size(); }
            uint                  get_num_dropped_coefficients() const { return mNumDroppedCoefficients; }
            real                  get_domain_measure() const { return mDomainMeasure; }
            real                  get_mass_sum() const { return mMassSum; }
            real                  get_min_diagonal() const { return mMinDiagonal; }
            real                  get_min_pivot() const { return mMinPivot; }
            real                  get_filter_gamma() const { return mFilterGamma; }
            const std::string&    get_adv_field_name() const { return mADVFieldName; }

            /** int_e B_A summed over the cells (= M 1), full ADV vector */
            Vector< real > get_coefficient_integrals() const;
        };

        //------------------------------------------------------------------------------
    }    // namespace wrk
}    // namespace moris
