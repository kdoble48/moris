/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * example_test_case.cpp
 *
 * LevelSet_ADV_Filter: B-spline level-set designs without and with the Kreissl-Maute filter of the B-spline
 * coefficients. Each case is its own TEST_CASE and ctest process (tags [EXA_LevelSet_ADV_Filter_<0..3>]); the
 * optimization cases 2 and 3 also carry [EXA_LevelSet_ADV_Filter], the tag of the gallery's showcase run.
 *
 * Cases 0 (no filter) and 1 (filter): sensitivity sweep of a plate with one hole
 *   - analytic dObjective/dADV and dConstraint/dADV agree with central finite differences for every ADV
 *     (relative to the largest analytic entry);
 *   - the criteria at the initial design match reference values; case 1 differs from case 0 (the fields use F psi).
 * Cases 2 (no filter) and 3 (filter): 100 MMA iterations of a cantilever, min strain energy s.t. volume <= 1.0
 *   - the final design meets the volume limit and its strain energy and volume match reference values.
 */

#include <catch.hpp>
#include <iomanip>

#include "cl_Logger.hpp"    // MRS/IOS/src
#include "HDF5_Tools.hpp"
#include "fn_norm.hpp"
#include "cl_MTK_Exodus_IO_Helper.hpp"    // MTK/src
#include <fstream>
#include <cstdio>

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

    // reference final strain energy and volume of the optimization cases 2 (no filter) and 3 (filter)
    static const real tReferenceFinalSE[ 2 ]     = { 31.913990251, 32.070175458 };
    static const real tReferenceFinalVolume[ 2 ] = { 0.99999978200, 0.99964353239 };

    //---------------------------------------------------------------

    // output file stem of a case (the deck's output_name())
    static std::string
    output_name( uint aCaseIndex )
    {
        return "levelset_adv_filter_" + std::to_string( aCaseIndex );
    }

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

        hid_t  tFileID = open_hdf5_file( output_name( aCaseIndex ) + ".hdf5" );
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

    // VIS writes one file per evaluation: levelset_adv_filter_<case>.exo.e-s.0000, .0001, ...
    static std::string
    exodus_step_file(
            uint aCaseIndex,
            uint aStep )
    {
        char tSuffix[ 16 ];
        std::snprintf( tSuffix, sizeof( tSuffix ), "%04u", aStep );
        return output_name( aCaseIndex ) + ".exo.e-s." + std::string( tSuffix );
    }

    static void
    check_optimization_results( uint aCaseIndex )
    {
        // the last evaluation's file holds the final design (MMA accepts its last trial)
        uint tLastStep = 0;
        while ( std::ifstream( exodus_step_file( aCaseIndex, tLastStep + 1 ) ).good() )
        {
            tLastStep++;
        }
        std::string tLastFile = exodus_step_file( aCaseIndex, tLastStep );
        REQUIRE( tLastStep > 100 );

        moris::mtk::Exodus_IO_Helper tExoIO( tLastFile.c_str(), 0, true, false );
        real                         tStrainEnergy = 0.0;
        real                         tVolume       = 0.0;
        for ( uint iGlobal = 0; iGlobal < 2; iGlobal++ )
        {
            std::string tName                                     = tExoIO.get_global_variable_name( iGlobal );
            ( tName == "StrainEnergy" ? tStrainEnergy : tVolume ) = tExoIO.get_global_variable( iGlobal, 0 );
        }

        MORIS_LOG_INFO( "LevelSet_ADV_Filter case %u: %u evaluations; final strain energy %.10e, volume %.10e (limit 1.0)",
                aCaseIndex,
                tLastStep,
                tStrainEnergy,
                tVolume );

        CHECK( tVolume <= 1.0 + 1e-3 );
        CHECK( tStrainEnergy == Approx( tReferenceFinalSE[ aCaseIndex - 2 ] ).epsilon( 1e-6 ) );
        CHECK( tVolume == Approx( tReferenceFinalVolume[ aCaseIndex - 2 ] ).epsilon( 1e-6 ) );
    }

    //---------------------------------------------------------------

    static void
    run_case( uint aCaseIndex )
    {
        gCaseIndex = aCaseIndex;

        if ( par_rank() == 0 )
        {
            std::string tRemove = "rm -f " + output_name( aCaseIndex ) + ".hdf5 " + output_name( aCaseIndex ) + ".exo* ADV_Alg_0_Iter_*.hdf5";
            MORIS_ERROR( std::system( tRemove.c_str() ) == 0,
                    "LevelSet_ADV_Filter: removing the output of an earlier run failed" );
        }

        int   argc       = 2;
        char  tString1[] = "";
        char  tString2[] = "LevelSet_ADV_Filter.so";
        char* argv[ 2 ]  = { tString1, tString2 };

        if ( par_size() == 1 )
        {
            int tRet = fn_WRK_Workflow_Main_Interface( argc, argv );
            REQUIRE( tRet == 0 );
            if ( aCaseIndex < 2 )
            {
                check_results( aCaseIndex );
            }
            else
            {
                check_optimization_results( aCaseIndex );
            }
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

    TEST_CASE( "LevelSet_ADV_Filter_2", "[moris],[example],[optimization],[EXA_LevelSet_ADV_Filter_2],[EXA_LevelSet_ADV_Filter]" )
    {
        run_case( 2 );
    }

    TEST_CASE( "LevelSet_ADV_Filter_3", "[moris],[example],[optimization],[EXA_LevelSet_ADV_Filter_3],[EXA_LevelSet_ADV_Filter]" )
    {
        run_case( 3 );
    }

}    // namespace exa_levelset_adv_filter
