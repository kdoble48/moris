/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_OPT_Problem.cpp
 *
 */

#include "cl_OPT_Problem.hpp"
#include "op_plus.hpp"
#include "fn_norm.hpp"
#include "fn_trans.hpp"
#include "fn_Parsing_Tools.hpp"
#include "cl_Communication_Tools.hpp"
#include "HDF5_Tools.hpp"

// Logger package
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"

#include "fn_stringify_matrix.hpp"

#include <algorithm>
#include <vector>
#include <cmath>

namespace moris::opt
{

    // -------------------------------------------------------------------------------------------------------------
    // Public functions
    // -------------------------------------------------------------------------------------------------------------

    Problem::Problem(
            Parameter_List&                        aParameterList,
            std::shared_ptr< Criteria_Interface >& aInterface )
    {
        // Set interface
        mInterface = aInterface;

        // Parameters: restart file name
        mRestartFile = aParameterList.get< std::string >( "restart_file" );

        // Parameters: gradient-explosion clip factor; <= 0 disables the clip
        if ( aParameterList.exists( "grad_clip_factor" ) )
        {
            mGradClipFactor = aParameterList.get< real >( "grad_clip_factor" );
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    Problem::~Problem()
    {
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Problem::initialize()
    {
        // Trace optimization problem
        Tracer tTracer( "OPT", "Problem", "Initialize" );

        // Initialize ADVs
        mInterface->initialize( mADVs, mLowerBounds, mUpperBounds, mIjklIds );

        // Override interface ADVs
        this->override_advs();

        // Log number of optimization variables
        MORIS_LOG_SPEC( "Number of optimization variables", mADVs.size() );

        // Diagnostic: non-finite entries handed back by the interface poison the optimizer
        // (NaN ADVs make every downstream consistency check fail). Detect them at the source.
        if ( par_rank() == 0 )
        {
            uint tNanADVs    = 0;
            uint tNanBounds  = 0;
            for ( uint iADV = 0; iADV < mADVs.size(); iADV++ )
            {
                if ( !std::isfinite( mADVs( iADV ) ) ) tNanADVs++;
                if ( !std::isfinite( mLowerBounds( iADV ) ) || !std::isfinite( mUpperBounds( iADV ) ) ) tNanBounds++;
            }
            if ( tNanADVs > 0 || tNanBounds > 0 )
            {
                MORIS_LOG_WARNING( "Problem::initialize - non-finite entries after interface initialize: %u ADVs, %u bounds (of %zu).",
                        tNanADVs, tNanBounds, mADVs.size() );
            }
        }

        MORIS_ERROR( mADVs.size() == mLowerBounds.size() and mADVs.size() == mUpperBounds.size(),
                "ADVs and its lower and upper bound vectors need to have same length.\n" );

        // Read ADVs from restart file
        if ( par_rank() == 0 && !mRestartFile.empty() )
        {
            this->read_restart_file();
        }

        // Initialize constraint types
        mConstraintTypes = this->get_constraint_types();

        MORIS_ERROR( this->check_constraint_order(), "The constraints are not ordered properly (Eq, Ineq)" );    // TODO move call to alg

        // Set number of objectives and constraints
        mNumObjectives  = 1;
        mNumConstraints = mConstraintTypes.numel();

        // Log number of objectives and constraints
        MORIS_LOG_SPEC( "Number of objectives ", mNumObjectives );
        MORIS_LOG_SPEC( "Number of constraints", mNumConstraints );
    }

    // -------------------------------------------------------------------------------------------------------------

    bool
    Problem::check_constraint_order()
    {
        // Checks that the constraint order is  1. equality constraints   (typ=0)
        //                                      2. inequality constraints (typ=1)

        int tSwitches = 1;
        for ( uint tConstraintIndex = 0; tConstraintIndex < mConstraintTypes.numel(); tConstraintIndex++ )
        {
            if ( mConstraintTypes( tConstraintIndex ) == tSwitches )
            {
                tSwitches -= 1;
            }
        }
        return ( tSwitches >= 0 );
    }

    // -------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Problem::get_objectives()
    {
        MORIS_ERROR( par_rank() == 0,
                "Problem::get_objectives - called by another processor but zero.\n" );

        // log objective
        MORIS_LOG_SPEC( "Objective", ios::stringify_log( mObjectives ) );

        return mObjectives;
    }

    // -------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Problem::get_constraints()
    {
        MORIS_ERROR( par_rank() == 0,
                "Problem::get_constraints - called by another processor but zero.\n" );

        // log constraints
        MORIS_LOG_SPEC( "Constraints", ios::stringify_log( mConstraints ) );

        return mConstraints;
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Problem::set_objectives_and_constraints(
            const Matrix< DDRMat >& aObjectives,
            const Matrix< DDRMat >& aConstraints )
    {
        // check for proper size of input matrices
        MORIS_ERROR( aObjectives.n_rows() == mObjectives.n_rows() and aObjectives.n_rows() == mObjectives.n_rows(),
                "Problem::set_objectives_and_constraints - input objective vector has incorrect size\n." );

        MORIS_ERROR( aConstraints.n_rows() == mConstraints.n_rows() and aConstraints.n_rows() == mConstraints.n_rows(),
                "Problem::set_objectives_and_constraints - input objective vector has incorrect size\n." );

        // copy objectives and constraints
        mObjectives  = aObjectives;
        mConstraints = aConstraints;
    }

    // -------------------------------------------------------------------------------------------------------------

    uint
    Problem::get_num_equality_constraints()
    {
        uint tNumEqualityConstraints = 0;

        // Loop over the constraint types matrix
        for ( uint tConstraintIndex = 0; tConstraintIndex < mConstraintTypes.numel(); tConstraintIndex++ )
        {
            if ( mConstraintTypes( tConstraintIndex ) == 0 )
            {
                tNumEqualityConstraints++;
            }
        }

        return tNumEqualityConstraints;
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Problem::compute_design_criteria( Vector< real >& aNewADVs )
    {
        // check for proper size of input ADV vector (on processor 0 only as all other processors
        // have zero-size ADV vector
        if ( par_rank() == 0 )
        {
            // Check that input ADV vector has correct size
            MORIS_ERROR( mADVs.size() == aNewADVs.size(),
                    "Problem::compute_design_criteria - size of ADV vectors does not match.\n" );
        }

        // compute design criteria
        mCriteria = mInterface->get_criteria( aNewADVs );

        // compute objective and constraints (on processor 0 only)
        if ( par_rank() == 0 )
        {
            // save the new ADV vector ( this is done after the criteria is computed since get-criteria might make some changes
            //  in the adv vector e.g reinitialization of the ADV vector)
            mADVs = aNewADVs;

            // log criteria and ADVs
            MORIS_LOG_SPEC( "Criteria", ios::stringify_log( mCriteria ) );

            if ( mADVs.size() > 0 )
            {
                MORIS_LOG_SPEC( "MinADV", mADVs.min() );
                MORIS_LOG_SPEC( "MaxADV", mADVs.max() );
            }

            // compute objective
            mObjectives = this->compute_objectives();

            MORIS_ASSERT( mObjectives.numel() == mNumObjectives,
                    "compute_design_criteria - Number of objectives is incorrect (%zu vs %i).\n",
                    mObjectives.numel(),
                    mNumObjectives );

            // compute constraints
            mConstraints = this->compute_constraints();

            MORIS_ASSERT( mConstraints.numel() == mNumConstraints,
                    "compute_design_criteria - Number of constraints is incorrect (%zu vs %i).\n",
                    mConstraints.numel(),
                    mNumConstraints );
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Problem::compute_design_criteria_gradients( const Vector< real >& aNewADVs )
    {
        // Check that input ADV vector is identical to the one stored in problem
        // as otherwise new criteria evaluation would be needed; this is done only
        // on processor 0 only all other processors have zero-size ADV vector
        if ( par_rank() == 0 )
        {
            for ( uint iADVIndex = 0; iADVIndex < mADVs.size(); iADVIndex++ )
            {
                MORIS_ERROR( std::abs( mADVs( iADVIndex ) - aNewADVs( iADVIndex ) ) <= mADVNormTolerance,
                        "Problem::compute_design_criteria_gradients - given ADV vector does not match ADVs used to compute design criteria.\n" );
            }
        }

        // compute derivatives of design criteria
        mInterface->get_dcriteria_dadv();

        // compute objective and constraints (on processor 0 only)
        if ( par_rank() == 0 )
        {
            // compute gradients of objective
            mObjectiveGradient =
                    this->compute_dobjective_dadv() +    //
                    this->compute_dobjective_dcriteria() * mInterface->get_dcriteria_dadv();

            MORIS_ASSERT( mObjectiveGradient.n_rows() == mNumObjectives and mObjectiveGradient.n_cols() == mADVs.size(),
                    "Problem::compute_design_criteria_gradients  - gradient of objective matrix has incorrect size.\n." );

            // compute gradients of constraints
            mConstraintGradient =
                    this->compute_dconstraint_dadv() +    //
                    this->compute_dconstraint_dcriteria() * mInterface->get_dcriteria_dadv();

            MORIS_ASSERT( mConstraintGradient.n_rows() == mNumConstraints and mConstraintGradient.n_cols() == mADVs.size(),
                    "Problem::compute_design_criteria_gradients  - gradient of constraint matrix has incorrect size.\n." );

            // ----- Gradient explosion detection and per-ADV clipping -----
            real tObjGradNorm = norm( mObjectiveGradient );
            real tConGradNorm = norm( mConstraintGradient );

            bool tGradClipped = false;

            // ----- Self-calibrating gradient-explosion clip (opt-in via grad_clip_factor) -----
            // Small XFEM cut cells make a few dIQI/dADV entries explode by many orders of magnitude
            // (dIQI/dPDV ~ 1/cell_volume as the interface nears a background node). The design
            // gradient is bimodal: most B-spline coeffs are EXACTLY 0 (inactive, no interface in
            // their support), a band is healthy active sensitivity, and a few explode. Clip each
            // entry to mGradClipFactor x the MEDIAN of the NONZERO entries -- the active-band scale.
            // This is (a) self-calibrating from the current gradient, so it needs no "previous
            // healthy" baseline -- dense seedings explode from iteration 1 and never establish one,
            // which is why a relative-jump test (ratio > threshold) never fired and the raw 1e8
            // gradient froze the design; and (b) immune to the inactive-zero floor that made a
            // P90-of-ALL-entries reference collapse to ~0 and crush the whole active gradient.
            // The cap also flattens legitimately large sensitivities, so the clip is disabled
            // unless explicitly requested (grad_clip_factor > 0).
            if ( mGradClipFactor > 0.0 )
            {
                auto tActiveScale = []( const Matrix< DDRMat >& aGrad ) -> real {
                    std::vector< real > tNz;
                    tNz.reserve( aGrad.n_rows() * aGrad.n_cols() );
                    for ( uint ii = 0; ii < aGrad.n_rows(); ++ii )
                        for ( uint jj = 0; jj < aGrad.n_cols(); ++jj )
                        {
                            real tv = std::abs( aGrad( ii, jj ) );
                            if ( tv > 0.0 ) tNz.push_back( tv );
                        }
                    if ( tNz.empty() ) return 0.0;
                    uint tMid = tNz.size() / 2;
                    std::nth_element( tNz.begin(), tNz.begin() + tMid, tNz.end() );
                    return tNz[ tMid ];
                };

                real tObjClipVal = mGradClipFactor * tActiveScale( mObjectiveGradient );
                real tConClipVal = mGradClipFactor * tActiveScale( mConstraintGradient );

                // degenerate (all-zero) gradient: disable clipping for it
                if ( !( tObjClipVal > 0.0 ) ) tObjClipVal = 1.0e30;
                if ( !( tConClipVal > 0.0 ) ) tConClipVal = 1.0e30;

                uint tObjClipCount = 0;
                uint tConClipCount = 0;

                for ( uint j = 0; j < mADVs.size(); j++ )
                {
                    for ( uint i = 0; i < mObjectiveGradient.n_rows(); i++ )
                    {
                        if ( std::abs( mObjectiveGradient( i, j ) ) > tObjClipVal )
                        {
                            real tSign                  = ( mObjectiveGradient( i, j ) > 0.0 ) ? 1.0 : -1.0;
                            mObjectiveGradient( i, j )  = tSign * tObjClipVal;
                            tObjClipCount++;
                        }
                    }
                    for ( uint i = 0; i < mConstraintGradient.n_rows(); i++ )
                    {
                        if ( std::abs( mConstraintGradient( i, j ) ) > tConClipVal )
                        {
                            real tSign                  = ( mConstraintGradient( i, j ) > 0.0 ) ? 1.0 : -1.0;
                            mConstraintGradient( i, j ) = tSign * tConClipVal;
                            tConClipCount++;
                        }
                    }
                }

                tGradClipped = ( tObjClipCount + tConClipCount ) > 0;

                if ( tGradClipped )
                {
                    MORIS_LOG_WARNING(
                            "Gradient clip: %u obj + %u con entries clipped (clip_val = %.4e / %.4e, "
                            "raw obj_norm = %.4e).",
                            tObjClipCount,
                            tConClipCount,
                            tObjClipVal,
                            tConClipVal,
                            tObjGradNorm );
                }
            }

            // ----- Gradient diagnostics logging -----
            real tObjGradMax = 0.0;
            real tConGradMax = 0.0;
            {
                uint tAtLower = 0;
                uint tAtUpper = 0;

                for ( uint j = 0; j < mADVs.size(); j++ )
                {
                    // Max absolute gradient entry
                    for ( uint i = 0; i < mObjectiveGradient.n_rows(); i++ )
                    {
                        tObjGradMax = std::max( tObjGradMax, std::abs( mObjectiveGradient( i, j ) ) );
                    }
                    for ( uint i = 0; i < mConstraintGradient.n_rows(); i++ )
                    {
                        tConGradMax = std::max( tConGradMax, std::abs( mConstraintGradient( i, j ) ) );
                    }

                    // Bound saturation (ADV within 1e-10 of bound)
                    if ( std::abs( mADVs( j ) - mLowerBounds( j ) ) < 1.0e-10 )
                    {
                        tAtLower++;
                    }
                    if ( std::abs( mADVs( j ) - mUpperBounds( j ) ) < 1.0e-10 )
                    {
                        tAtUpper++;
                    }
                }

                MORIS_LOG_INFO(
                    "GradDebug: obj_norm=%.4e  con_norm=%.4e  obj_max=%.4e  con_max=%.4e  at_lower=%u  at_upper=%u  n_advs=%zu",
                    tObjGradNorm,
                    tConGradNorm,
                    tObjGradMax,
                    tConGradMax,
                    tAtLower,
                    tAtUpper,
                    mADVs.size() );
            }

            // Update the reference gradient for next-iteration explosion detection. CRITICAL:
            // when clipping fired this iteration, do NOT adopt the (artificially shrunk) clipped
            // norm as the reference. Doing so makes the clip threshold (mPrevObjGradNorm/sqrt(N))
            // spiral toward zero over successive explosions, eventually clipping the *healthy*
            // optimization gradient to ~0 and freezing the design at its initial topology. Keep the
            // last healthy reference instead, so the clip caps explosive (small-cut-cell) ADVs at
            // the healthy gradient scale while letting the legitimate gradient drive the design.
            if ( !tGradClipped )
            {
                mPrevObjectiveGradient  = mObjectiveGradient;
                mPrevConstraintGradient = mConstraintGradient;
                mPrevObjGradNorm        = norm( mObjectiveGradient );
                mPrevConGradNorm        = norm( mConstraintGradient );
                mPrevObjGradMax         = tObjGradMax;    // healthy per-entry scale for absolute clip
                mPrevConGradMax         = tConGradMax;
            }
            mHasPrevGradients = true;
        }
    }

    // -------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Problem::get_objective_gradients()
    {
        MORIS_ERROR( par_rank() == 0,
                "Problem::get_objective_gradients - called by another processor but zero.\n" );

        return mObjectiveGradient;
    }

    // -------------------------------------------------------------------------------------------------------------

    const Matrix< DDRMat >&
    Problem::get_constraint_gradients()
    {
        MORIS_ERROR( par_rank() == 0,
                "Problem::get_constraint_gradients - called by another processor but zero.\n" );

        return mConstraintGradient;
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Problem::scale_solution()
    {
        // TODO Need to decide on a framework for the scaling of the solution
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Problem::update_problem()
    {
        // TODO Need to decide on a framework for the update of the problem
    }

    // -------------------------------------------------------------------------------------------------------------

    void
    Problem::read_restart_file()
    {
        MORIS_LOG_INFO( "Reading ADVs from restart file: %s", mRestartFile.c_str() );

        // Open open restart file
        // Note: only processor 0 reads this file; therefore no parallel file name extension is used
        hid_t tFileID = open_hdf5_file( mRestartFile, false );

        // Define matrix in which to read restart ADVs
        Vector< real >  tRestartADVs;
        Vector< real >  tRestartUpperBounds;
        Vector< real >  tRestartLowerBounds;
        Matrix< IdMat > tIjklIdsFile;
        sint            tParSizeFile = 0;

        // Read ADVS from restart file
        herr_t tStatus = 0;
        load_vector_from_hdf5_file( tFileID, "ADVs", tRestartADVs.data(), tStatus );
        load_vector_from_hdf5_file( tFileID, "UpperBounds", tRestartUpperBounds.data(), tStatus );
        load_vector_from_hdf5_file( tFileID, "LowerBounds", tRestartLowerBounds.data(), tStatus );

        load_scalar_from_hdf5_file( tFileID, "NumProcs", tParSizeFile, tStatus );

        // FIXME this function also should work if you want to restart with same number procs but new proc layout
        // in this case this if criteria should be deleted
        if ( par_size() != tParSizeFile )
        {
            load_matrix_from_hdf5_file( tFileID, "IjklIds", tIjklIdsFile, tStatus );

            MORIS_ERROR( tIjklIdsFile.numel() == tRestartADVs.size(), "restart adv vector and ijkl ID vector not of same size" );
            MORIS_ERROR( tIjklIdsFile.numel() == mIjklIds.numel(), "restart ijkl ID vector and new ijkl ID vector not of same size" );
        }

        // Close restart file
        close_hdf5_file( tFileID );

        // Check for matching sizes
        MORIS_ERROR( tRestartADVs.size() == mADVs.size(),
                "Number of restart ADVS does not match.\n" );
        MORIS_ERROR( tRestartUpperBounds.size() == mUpperBounds.size(),
                "Number of restart upper bounds of ADVs does not match.\n" );
        MORIS_ERROR( tRestartLowerBounds.size() == mLowerBounds.size(),
                "Number of restart lower bounds of ADVs does not match.\n" );

        // reorder adv for new proc layout. FIXME see above
        if ( par_size() != tParSizeFile )
        {
            moris::map< moris_id, sint > tIjklIdToPosMap;
            for ( uint Ik = 0; Ik < mIjklIds.numel(); Ik++ )
            {
                if ( mIjklIds( Ik ) != gNoID )
                {
                    tIjklIdToPosMap[ mIjklIds( Ik ) ] = Ik;
                }
            }

            for ( uint Ik = 0; Ik < tIjklIdsFile.numel(); Ik++ )
            {
                if ( tIjklIdToPosMap.key_exists( tIjklIdsFile( Ik ) ) )
                {
                    sint tIndxed            = tIjklIdToPosMap.find( tIjklIdsFile( Ik ) );
                    mADVs( tIndxed )        = tRestartADVs( Ik );
                    mUpperBounds( tIndxed ) = tRestartUpperBounds( Ik );
                    mLowerBounds( tIndxed ) = tRestartLowerBounds( Ik );
                }
                else
                {
                    mADVs( Ik )        = tRestartADVs( Ik );
                    mUpperBounds( Ik ) = tRestartUpperBounds( Ik );
                    mLowerBounds( Ik ) = tRestartLowerBounds( Ik );
                }
            }

            MORIS_ERROR( mADVs.max() != MORIS_REAL_MAX, "Restart ADV is MORIS_REAL_MAX" );
            MORIS_ERROR( mUpperBounds.max() != MORIS_REAL_MAX, "Restart upper bound is MORIS_REAL_MAX" );
            MORIS_ERROR( mLowerBounds.max() != MORIS_REAL_MAX, "Restart lower bound is MORIS_REAL_MAX" );
        }

        // Copy restart vectors on member variables
        mADVs        = tRestartADVs;
        mUpperBounds = tRestartUpperBounds;
        mLowerBounds = tRestartLowerBounds;

        mRestartFile = "";
    }
}    // namespace moris::opt
