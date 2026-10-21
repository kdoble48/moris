/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 * LevelSet_ADV_Filter: sensitivity sweep of a plate with a B-spline level-set hole, without (case 0) and with (case 1)
 * the Kreissl-Maute filter of the B-spline coefficients. Each case is its own TEST_CASE and ctest process
 * (tags [EXA_LevelSet_ADV_Filter_0], [EXA_LevelSet_ADV_Filter_1]).
 *
 * Checks per case:
 *   - analytic dObjective/dADV and dConstraint/dADV agree with central finite differences for every ADV
 *     (relative to the largest analytic entry);
 *   - the criteria at the initial design match reference values; case 1 differs from case 0 (the fields use F psi).
 */

#include <catch.hpp>
#include <iomanip>

#include "cl_Logger.hpp"    // MRS/IOS/src
#include "HDF5_Tools.hpp"
#include "fn_norm.hpp"

#include "EXA_Globals.hpp"    // shared deck-visible globals (dlopen ABI; see doc/internal/EXA_RUNNER_RFC.md)

using namespace moris;

//---------------------------------------------------------------

// defined at global scope in WRK - must be declared outside the example namespace
int fn_WRK_Workflow_Main_Interface( int argc, char* argv[] );

//---------------------------------------------------------------

namespace exa_levelset_adv_filter
{
    // reference criteria at the initial design: { strain energy, volume } for case 0 (no filter) and case 1 (r = 1.5 h)
    // (recorded 2026-10-20, feature/ls-adv-filter; a failing check prints the new value)
    static const real tReferenceObjective[ 2 ]  = { 23.7730481965, 22.7173375688 };
    static const real tReferenceConstraint[ 2 ] = { 1.7243383991, 1.7568370074 };

    //---------------------------------------------------------------

    static void
    check_results( uint aCaseIndex )
    {
        Matrix< DDRMat > tObjective;
        Matrix< DDRMat > tConstraint;
        Matrix< DDRMat > tObjectiveAnalytic;
        Matrix< DDRMat > tConstraintAnalytic;
        Matrix< DDRMat > tObjectiveFD;
        Matrix< DDRMat > tConstraintFD;

        hid_t  tFileID = open_hdf5_file( "levelset_adv_filter.hdf5" );
        herr_t tStatus = 0;
        load_matrix_from_hdf5_file( tFileID, "objectives eval_1-1", tObjective, tStatus );
        load_matrix_from_hdf5_file( tFileID, "constraints eval_1-1", tConstraint, tStatus );
        load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 analytical", tObjectiveAnalytic, tStatus );
        load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 analytical", tConstraintAnalytic, tStatus );
        load_matrix_from_hdf5_file( tFileID, "objective_gradients eval_1-1 epsilon_1-1 fd_central", tObjectiveFD, tStatus );
        load_matrix_from_hdf5_file( tFileID, "constraint_gradients eval_1-1 epsilon_1-1 fd_central", tConstraintFD, tStatus );
        close_hdf5_file( tFileID );

        if ( gPrintReferenceValues )
        {
            std::cout << "LevelSet_ADV_Filter case " << aCaseIndex << ": objective " << std::scientific << std::setprecision( 15 )
                      << tObjective( 0 ) << ", constraint " << tConstraint( 0 ) << '\n';
        }

        // criteria at the initial design
        CHECK( tObjective( 0 ) == Approx( tReferenceObjective[ aCaseIndex ] ).epsilon( 1e-8 ) );
        CHECK( tConstraint( 0 ) == Approx( tReferenceConstraint[ aCaseIndex ] ).epsilon( 1e-8 ) );

        // gradients: every ADV, compared to the largest analytic entry of the criterion
        REQUIRE( tObjectiveAnalytic.numel() == tObjectiveFD.numel() );
        REQUIRE( tConstraintAnalytic.numel() == tConstraintFD.numel() );
        REQUIRE( tObjectiveAnalytic.numel() > 0 );

        real tObjectiveScale  = std::max( norm( tObjectiveAnalytic ) / std::sqrt( (real)tObjectiveAnalytic.numel() ), MORIS_REAL_EPS );
        real tConstraintScale = std::max( norm( tConstraintAnalytic ) / std::sqrt( (real)tConstraintAnalytic.numel() ), MORIS_REAL_EPS );
        for ( uint iADV = 0; iADV < tObjectiveAnalytic.numel(); iADV++ )
        {
            tObjectiveScale  = std::max( tObjectiveScale, std::abs( tObjectiveAnalytic( iADV ) ) );
            tConstraintScale = std::max( tConstraintScale, std::abs( tConstraintAnalytic( iADV ) ) );
        }

        real tMaxObjectiveError  = 0.0;
        real tMaxConstraintError = 0.0;
        uint tNumNonzero         = 0;
        for ( uint iADV = 0; iADV < tObjectiveAnalytic.numel(); iADV++ )
        {
            real tObjectiveError  = std::abs( tObjectiveAnalytic( iADV ) - tObjectiveFD( iADV ) ) / tObjectiveScale;
            real tConstraintError = std::abs( tConstraintAnalytic( iADV ) - tConstraintFD( iADV ) ) / tConstraintScale;
            tMaxObjectiveError    = std::max( tMaxObjectiveError, tObjectiveError );
            tMaxConstraintError   = std::max( tMaxConstraintError, tConstraintError );
            tNumNonzero += std::abs( tConstraintAnalytic( iADV ) ) > 1e-12 ? 1 : 0;
            CHECK( tObjectiveError < 1e-5 );
            CHECK( tConstraintError < 1e-5 );
        }
        MORIS_LOG_INFO(
                "LevelSet_ADV_Filter case %u: %u ADVs, %u with a nonzero volume gradient; max |analytic - FD| / max|analytic|: "
                "objective %.2e, constraint %.2e",
                aCaseIndex,
                (uint)tObjectiveAnalytic.numel(),
                tNumNonzero,
                tMaxObjectiveError,
                tMaxConstraintError );
    }

    //---------------------------------------------------------------

    static void
    run_case( uint aCaseIndex )
    {
        gCaseIndex = aCaseIndex;

        if ( par_rank() == 0 )
        {
            MORIS_ERROR( std::system( "rm -f levelset_adv_filter.hdf5" ) == 0, "LevelSet_ADV_Filter: removing the sweep file failed" );
        }

        int   argc       = 2;
        char  tString1[] = "";
        char  tString2[] = "LevelSet_ADV_Filter.so";
        char* argv[ 2 ]  = { tString1, tString2 };

        if ( par_size() == 1 )
        {
            int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );
            REQUIRE( tRet == 0 );
            check_results( aCaseIndex );
        }
    }

    //---------------------------------------------------------------

    TEST_CASE( "LevelSet_ADV_Filter_0", "[moris],[example],[optimization],[sweep],[EXA_LevelSet_ADV_Filter_0]" )
    {
        run_case( 0 );
    }

    TEST_CASE( "LevelSet_ADV_Filter_1", "[moris],[example],[optimization],[sweep],[EXA_LevelSet_ADV_Filter_1]" )
    {
        run_case( 1 );
    }

}    // namespace exa_levelset_adv_filter
