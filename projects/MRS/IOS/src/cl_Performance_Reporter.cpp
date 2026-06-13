/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Performance_Reporter.cpp
 *
 */

#include "cl_Performance_Reporter.hpp"

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <vector>

#include <mpi.h>

#include "cl_Communication_Tools.hpp"    // COM/src   (par_rank, par_size)
#include "cl_Json_Object.hpp"            // MRS/IOS/src
#include "cl_Logger.hpp"                 // MRS/IOS/src (MORIS_LOG)

// global instance - defined here so any executable linking IOS gets it
moris::Performance_Reporter gPerfReporter;

namespace moris
{
    //------------------------------------------------------------------------------

    namespace
    {
        // collective reduction of a single scalar across all ranks
        real
        reduce_all( real aLocal, MPI_Op aOp )
        {
            real tGlobal = aLocal;
            MPI_Allreduce( &aLocal, &tGlobal, 1, MPI_DOUBLE, aOp, MPI_COMM_WORLD );
            return tGlobal;
        }

        // display row pairing an accumulated record with its cross-rank wall stats
        struct Display_Row
        {
            Module_Performance_Record mRecord;
            real                      mWallMax = 0.0;
            real                      mWallMin = 0.0;
        };
    }    // namespace

    //------------------------------------------------------------------------------

    void
    Performance_Reporter::initialize( sint aLevel, const std::string& aJsonPath )
    {
        mEnabled  = true;
        mLevel    = std::min< sint >( std::max< sint >( aLevel, 0 ), 3 );
        mJsonPath = aJsonPath;

        // mark the source of the memory numbers reported (see Logger::current_memory_mb)
#ifdef WITHGPERFTOOLS
        mMemSource = "tcmalloc";
#else
        mMemSource = "proc_rss";
#endif

        mRunWallStart = std::chrono::system_clock::now();
    }

    //------------------------------------------------------------------------------

    bool
    Performance_Reporter::include_in_report( uint aDepth ) const
    {
        switch ( mLevel )
        {
            case 0:
                return false;
            case 1:
                return aDepth <= 2;
            case 2:
                return aDepth <= 4;
            default:
                return true;
        }
    }

    //------------------------------------------------------------------------------

    void
    Performance_Reporter::record_region(
            const std::string& aEntity,
            const std::string& aType,
            const std::string& aAction,
            uint               aDepth,
            real               aWallSeconds,
            real               aMemDeltaMb )
    {
        // nothing to accumulate when disabled or only run totals are requested
        if ( !mEnabled || mLevel <= 0 )
        {
            return;
        }

        const std::string tKey = aEntity + "|" + aType + "|" + aAction;

        Module_Performance_Record& tRec = mRecords[ tKey ];

        // initialize identifying fields on first visit
        if ( tRec.mVisitCount == 0 )
        {
            tRec.mEntity = aEntity;
            tRec.mType   = aType;
            tRec.mAction = aAction;
            tRec.mDepth  = aDepth;
        }

        tRec.mVisitCount++;
        tRec.mWallTimeSum += aWallSeconds;
        tRec.mMemDeltaSum += aMemDeltaMb;

        if ( aMemDeltaMb > tRec.mMemDeltaPeak )
        {
            tRec.mMemDeltaPeak = aMemDeltaMb;
        }
    }

    //------------------------------------------------------------------------------

    void
    Performance_Reporter::update_peak_rss( real aMemMb )
    {
        if ( !mEnabled )
        {
            return;
        }

        if ( aMemMb > mPeakMemMb )
        {
            mPeakMemMb = aMemMb;
        }
    }

    //------------------------------------------------------------------------------

    void
    Performance_Reporter::finalize()
    {
        // no-op if a run was never started (e.g. unit-test executables)
        if ( !mEnabled )
        {
            return;
        }

        // ----- run-level aggregation (collective on all ranks) -----
        std::chrono::duration< double > tWall = std::chrono::system_clock::now() - mRunWallStart;

        real tRunWallLocal = tWall.count();
        real tRunWallMax   = reduce_all( tRunWallLocal, MPI_MAX );
        real tPeakMemMax   = reduce_all( mPeakMemMb, MPI_MAX );

        // ----- check that all ranks recorded the same number of regions -----
        real tCountLocal = (real)mRecords.size();
        bool tAligned    = ( reduce_all( tCountLocal, MPI_MAX ) == reduce_all( tCountLocal, MPI_MIN ) );

        // ----- per-region cross-rank wall-time stats (collective, in sorted-key lockstep) -----
        std::vector< Display_Row > tRows;
        tRows.reserve( mRecords.size() );

        for ( auto& tEntry : mRecords )
        {
            Display_Row tRow;
            tRow.mRecord = tEntry.second;

            if ( tAligned )
            {
                tRow.mWallMax = reduce_all( tEntry.second.mWallTimeSum, MPI_MAX );
                tRow.mWallMin = reduce_all( tEntry.second.mWallTimeSum, MPI_MIN );
            }
            else
            {
                // ranks diverged - fall back to local values, no cross-rank min/max
                tRow.mWallMax = tEntry.second.mWallTimeSum;
                tRow.mWallMin = tEntry.second.mWallTimeSum;
            }

            tRows.push_back( tRow );
        }

        // only rank 0 emits the report; all collective calls are above this guard
        if ( par_rank() != 0 )
        {
            return;
        }

        // drop regions filtered out by the current granularity level and sort by cost
        tRows.erase(
                std::remove_if( tRows.begin(), tRows.end(),
                        [ & ]( const Display_Row& aRow ) { return !include_in_report( aRow.mRecord.mDepth ); } ),
                tRows.end() );

        std::sort( tRows.begin(), tRows.end(),
                []( const Display_Row& aL, const Display_Row& aR ) {
                    return aL.mRecord.mWallTimeSum > aR.mRecord.mWallTimeSum;
                } );

        // ----- console summary table -----
        std::ostringstream tOut;
        tOut << "\n";
        tOut << "==================== MORIS Performance Report (level " << mLevel << ") ====================\n";

        if ( mLevel > 0 )
        {
            tOut << std::left << std::setw( 38 ) << "Module (entity - action)"
                 << std::right << std::setw( 7 ) << "Visits"
                 << std::setw( 12 ) << "Wall[s]"
                 << std::setw( 18 ) << "Wall max/min"
                 << std::setw( 13 ) << "MemD[MB]"
                 << std::setw( 13 ) << "MemDpeak\n";
            tOut << "-------------------------------------------------------------------------------------------------------\n";

            for ( const Display_Row& tRow : tRows )
            {
                const Module_Performance_Record& tRec = tRow.mRecord;

                std::string tLabel = tRec.mEntity + " - " + tRec.mAction;
                if ( tLabel.size() > 37 )
                {
                    tLabel = tLabel.substr( 0, 34 ) + "...";
                }

                std::ostringstream tMaxMin;
                tMaxMin << std::fixed << std::setprecision( 2 ) << tRow.mWallMax << "/" << tRow.mWallMin;

                // clamp the reported growth at zero (negative deltas are allocator noise)
                real tGrowth = tRec.mMemDeltaSum > 0.0 ? tRec.mMemDeltaSum : 0.0;

                tOut << std::left << std::setw( 38 ) << tLabel
                     << std::right << std::setw( 7 ) << tRec.mVisitCount
                     << std::setw( 12 ) << std::fixed << std::setprecision( 3 ) << tRec.mWallTimeSum
                     << std::setw( 18 ) << tMaxMin.str()
                     << std::setw( 13 ) << std::fixed << std::setprecision( 1 ) << tGrowth
                     << std::setw( 12 ) << std::fixed << std::setprecision( 1 ) << tRec.mMemDeltaPeak
                     << "\n";
            }

            tOut << "-------------------------------------------------------------------------------------------------------\n";
        }

        tOut << "Run total wall (max over ranks): " << std::fixed << std::setprecision( 3 ) << tRunWallMax << " s"
             << "   |   Peak memory (max over ranks): " << std::fixed << std::setprecision( 1 ) << tPeakMemMax
             << " MB [" << mMemSource << "]\n";
        tOut << "Report written to: " << mJsonPath << "\n";
        tOut << "=======================================================================================================\n";

        // emit to console and ASCII log
        MORIS_LOG( "%s", tOut.str().c_str() );

        // ----- JSON report -----
        Json tRoot;
        tRoot.put( "level", mLevel );
        tRoot.put( "run_total_wall_s", tRunWallMax );
        tRoot.put( "peak_mem_mb", tPeakMemMax );
        tRoot.put( "mem_source", mMemSource );
        tRoot.put( "ranks_aligned", tAligned );

        Json tModules;
        for ( const Display_Row& tRow : tRows )
        {
            const Module_Performance_Record& tRec = tRow.mRecord;

            Json tModule;
            tModule.put( "entity", tRec.mEntity );
            tModule.put( "type", tRec.mType );
            tModule.put( "action", tRec.mAction );
            tModule.put( "depth", tRec.mDepth );
            tModule.put( "visits", tRec.mVisitCount );
            tModule.put( "wall_s", tRec.mWallTimeSum );
            tModule.put( "wall_s_max", tRow.mWallMax );
            tModule.put( "wall_s_min", tRow.mWallMin );
            tModule.put( "mem_delta_mb", tRec.mMemDeltaSum );
            tModule.put( "mem_delta_peak_mb", tRec.mMemDeltaPeak );

            tModules.push_back( { "", tModule } );
        }
        tRoot.add_child( "modules", tModules );

        write_json( mJsonPath, tRoot );
    }

    //------------------------------------------------------------------------------

}    // namespace moris
