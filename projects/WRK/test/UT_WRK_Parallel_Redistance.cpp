/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * UT_WRK_Parallel_Redistance.cpp
 *
 * MB-113: on several processors the in-place interface redistance (sdf_reinit_mode = "interface",
 * sdf_reinit_every_step) measured each node's distance to the interface facets of its own processor only, wrote only
 * processor 0's coefficients into the ADV vector, and re-cut on processor 0 alone (MPI_ERR_TRUNCATE in the next
 * collective); sdf_preserve_interface stopped with "serial only".
 *
 * The workflow is called as the optimizer calls it (the ADVs on processor 0, an empty vector on the others) on the deck
 * WRK_Input_4.cpp: one evaluation with the redistance suppressed (the criteria of the design psi), one with the
 * redistance (it replaces the ADVs by phi = R(psi) and evaluates phi). The redistanced ADV vector and the criteria are
 * compared with the serial values (the 1-processor run of this test), for the plain redistance, the interface-preserving
 * redistance and the redistance with Laplacian smoothing. Run it under a timeout: a hang is a failure mode.
 *
 */

#include "catch.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <vector>

#include "cl_Communication_Tools.hpp"
#include "cl_Logger.hpp"
#include "cl_Tracer.hpp"
#include "paths.hpp"

#include "cl_WRK_Performer_Manager.hpp"
#include "cl_WRK_Workflow_HMR_XTK.hpp"
#include "cl_Library_Factory.hpp"

namespace moris::wrk
{
    namespace
    {
        // a signature of the ADV vector that does not depend on the order of the ADVs (the order follows the partition)
        struct ADV_Signature
        {
            uint mSize       = 0;
            real mSum        = 0.0;
            real mSumSquares = 0.0;
            real mSumCubes   = 0.0;
            real mRankSum    = 0.0;    // sum_k k * c_(k) over the sorted values
            real mMin        = 0.0;
            real mMax        = 0.0;
        };

        ADV_Signature
        signature( const Vector< real >& aADVs )
        {
            std::vector< real > tSorted( aADVs.size() );
            for ( uint iADV = 0; iADV < aADVs.size(); iADV++ )
            {
                tSorted[ iADV ] = aADVs( iADV );
            }
            std::sort( tSorted.begin(), tSorted.end() );

            ADV_Signature tSignature;
            tSignature.mSize = tSorted.size();
            for ( uint iADV = 0; iADV < tSorted.size(); iADV++ )
            {
                real tValue = tSorted[ iADV ];
                tSignature.mSum += tValue;
                tSignature.mSumSquares += tValue * tValue;
                tSignature.mSumCubes += tValue * tValue * tValue;
                tSignature.mRankSum += ( iADV + 1.0 ) * tValue;
            }
            if ( tSorted.size() > 0 )
            {
                tSignature.mMin = tSorted.front();
                tSignature.mMax = tSorted.back();
            }
            return tSignature;
        }

        struct Result
        {
            Vector< real > mCriteriaPsi;    // criteria of the design before the redistance
            Vector< real > mCriteriaPhi;    // criteria of the redistanced design
            ADV_Signature  mPsi;
            ADV_Signature  mPhi;
            real           mMaxChange = 0.0;    // max |phi - psi| over the ADVs
        };

        // the ADVs as the optimizer hands them to the workflow: complete on processor 0, empty on the other processors
        Result
        run( bool aPreserve, uint aSmoothing, const std::string& aLabel )
        {
            setenv( "WRK_INPUT_4_PRESERVE", aPreserve ? "1" : "0", 1 );
            setenv( "WRK_INPUT_4_SMOOTHING", std::to_string( aSmoothing ).c_str(), 1 );

            std::string tInputFilePath = moris::get_moris_bin_dir() + "/lib/WRK_Input_4.so";

            Library_Factory               tLibraryFactory;
            std::shared_ptr< Library_IO > tLibrary = tLibraryFactory.create_Library( Library_Type::STANDARD );
            tLibrary->load_parameter_list( tInputFilePath, File_Type::SO_FILE );
            tLibrary->finalize();

            Performer_Manager tPerformerManager( tLibrary );
            Workflow_HMR_XTK  tWorkflow( &tPerformerManager );

            Vector< real >  tADVs;
            Vector< real >  tLowerBounds;
            Vector< real >  tUpperBounds;
            Matrix< IdMat > tIjklIDs;
            tWorkflow.initialize( tADVs, tLowerBounds, tUpperBounds, tIjklIDs );

            // 17 x 17 linear B-spline coefficients on processor 0, none elsewhere
            REQUIRE( tADVs.size() == ( par_rank() == 0 ? 289 : 0 ) );

            // the workflow reads the optimization iteration from the optimizer's log entry; the redistance starts at 1
            Tracer tOptimizerTracer( "OPT", "Manager", "Perform" );
            gLogger.iterate();

            Result               tResult;
            const Vector< real > tPsi = tADVs;
            tResult.mPsi              = signature( tPsi );

            // psi: redistance suppressed (as for the optimizer's inner trials)
            tWorkflow.set_reinit_suppressed( true );
            tResult.mCriteriaPsi = tWorkflow.get_criteria( tADVs );
            REQUIRE( tADVs.size() == tPsi.size() );
            for ( uint iADV = 0; iADV < tPsi.size(); iADV++ )
            {
                REQUIRE( tADVs( iADV ) == tPsi( iADV ) );
            }

            // phi = R(psi): the workflow writes phi into the ADVs and evaluates it
            tWorkflow.set_reinit_suppressed( false );
            tResult.mCriteriaPhi = tWorkflow.get_criteria( tADVs );
            REQUIRE( tADVs.size() == tPsi.size() );
            tResult.mPhi = signature( tADVs );
            for ( uint iADV = 0; iADV < tPsi.size(); iADV++ )
            {
                tResult.mMaxChange = std::max( tResult.mMaxChange, std::abs( tADVs( iADV ) - tPsi( iADV ) ) );
            }

            // a NaN criterion would restart the optimizer; on processors > 0 it shows a processor that left the others
            REQUIRE( tResult.mCriteriaPsi.size() == 2 );
            REQUIRE( tResult.mCriteriaPhi.size() == 2 );
            CHECK( std::isfinite( tResult.mCriteriaPsi( 0 ) ) );
            CHECK( std::isfinite( tResult.mCriteriaPsi( 1 ) ) );
            CHECK( std::isfinite( tResult.mCriteriaPhi( 0 ) ) );
            CHECK( std::isfinite( tResult.mCriteriaPhi( 1 ) ) );

            if ( par_rank() == 0 )
            {
                std::cout << std::setprecision( 15 ) << "MB-113 " << aLabel << " on " << par_size() << " processor(s):\n"
                          << "  criteria psi " << tResult.mCriteriaPsi( 0 ) << " " << tResult.mCriteriaPsi( 1 ) << '\n'
                          << "  criteria phi " << tResult.mCriteriaPhi( 0 ) << " " << tResult.mCriteriaPhi( 1 ) << '\n'
                          << "  phi: sum " << tResult.mPhi.mSum << " squares " << tResult.mPhi.mSumSquares << " cubes "
                          << tResult.mPhi.mSumCubes << " rank sum " << tResult.mPhi.mRankSum << " min " << tResult.mPhi.mMin
                          << " max " << tResult.mPhi.mMax << " max|phi - psi| " << tResult.mMaxChange << '\n';
            }
            return tResult;
        }

        // the serial values (1-processor run of this test)
        struct Serial
        {
            real mCriteriaPsi[ 2 ];
            real mCriteriaPhi[ 2 ];
            real mSum;
            real mSumSquares;
            real mSumCubes;
            real mRankSum;
            real mMin;
            real mMax;
            real mMaxChange;
        };

        void
        check_serial_values( const Result& aResult, const Serial& aSerial )
        {
            if ( par_rank() != 0 )
            {
                return;
            }
            // the L2 projection and the forward problem are solved by direct solvers: round-off only
            const real tEps = 1.0e-9;
            CHECK( aResult.mCriteriaPsi( 0 ) == Approx( aSerial.mCriteriaPsi[ 0 ] ).epsilon( tEps ) );
            CHECK( aResult.mCriteriaPsi( 1 ) == Approx( aSerial.mCriteriaPsi[ 1 ] ).epsilon( tEps ) );
            CHECK( aResult.mCriteriaPhi( 0 ) == Approx( aSerial.mCriteriaPhi[ 0 ] ).epsilon( tEps ) );
            CHECK( aResult.mCriteriaPhi( 1 ) == Approx( aSerial.mCriteriaPhi[ 1 ] ).epsilon( tEps ) );
            CHECK( aResult.mPhi.mSize == 289 );
            CHECK( aResult.mPhi.mSum == Approx( aSerial.mSum ).epsilon( tEps ) );
            CHECK( aResult.mPhi.mSumSquares == Approx( aSerial.mSumSquares ).epsilon( tEps ) );
            CHECK( aResult.mPhi.mSumCubes == Approx( aSerial.mSumCubes ).epsilon( tEps ) );
            CHECK( aResult.mPhi.mRankSum == Approx( aSerial.mRankSum ).epsilon( tEps ) );
            CHECK( aResult.mPhi.mMin == Approx( aSerial.mMin ).epsilon( tEps ) );
            CHECK( aResult.mPhi.mMax == Approx( aSerial.mMax ).epsilon( tEps ) );
            CHECK( aResult.mMaxChange == Approx( aSerial.mMaxChange ).epsilon( tEps ) );
        }
    }    // namespace

    TEST_CASE( "MB-113 in-place interface redistance on several processors", "[WRK],[sdf_redistance_parallel],[MB-113]" )
    {
        SECTION( "plain redistance" )
        {
            Result tResult = run( false, 0, "plain redistance" );

            // the redistance changes the design, and with it the cut and the criteria
            CHECK( ( par_rank() != 0 || tResult.mMaxChange > 1.0e-2 ) );
            CHECK( ( par_rank() != 0 || std::abs( tResult.mCriteriaPhi( 1 ) - tResult.mCriteriaPsi( 1 ) ) > 1.0e-8 ) );

            check_serial_values( tResult,
                    { { 0.17528043994646, 0.220453383866188 },
                            { 0.175927647604536, 0.219508772795647 },
                            -38.6918906795042,
                            10.7115028001155,
                            -2.49209520834741,
                            -2346.81042834885,
                            -0.3,
                            0.207217815769316,
                            0.148235162093118 } );
        }

        SECTION( "interface-preserving redistance" )
        {
            Result tResult = run( true, 0, "interface-preserving redistance" );

            // the far field is redistanced, the cut and with it the criteria are those of psi (MB-027)
            if ( par_rank() == 0 )
            {
                CHECK( tResult.mMaxChange > 1.0e-2 );
                CHECK( tResult.mCriteriaPhi( 0 ) == Approx( tResult.mCriteriaPsi( 0 ) ).epsilon( 1.0e-12 ) );
                CHECK( tResult.mCriteriaPhi( 1 ) == Approx( tResult.mCriteriaPsi( 1 ) ).epsilon( 1.0e-12 ) );
            }

            check_serial_values( tResult,
                    { { 0.17528043994646, 0.220453383866188 },
                            { 0.17528043994646, 0.220453383866188 },
                            -39.4989792973614,
                            11.8706308997911,
                            -2.5396851304924,
                            -2319.54214810031,
                            -0.3,
                            0.207217815769316,
                            0.148235162093118 } );
        }

        SECTION( "redistance with Laplacian smoothing" )
        {
            Result tResult = run( false, 3, "redistance with 3 smoothing iterations" );

            check_serial_values( tResult,
                    { { 0.17528043994646, 0.220453383866188 },
                            { 0.185705459515691, 0.204003959284502 },
                            -38.0998253394164,
                            9.56783257787638,
                            -2.14312599580639,
                            -2561.36701211414,
                            -0.299650980752602,
                            0.160523515286414,
                            0.143310853007463 } );
        }
    }
}    // namespace moris::wrk
