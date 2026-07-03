/*
 * Copyright (c) 2022 University of Colorado
 * Licensed under the MIT license. See LICENSE.txt file in the MORIS root for details.
 *
 *------------------------------------------------------------------------------------
 *
 * cl_Performance_Reporter.hpp
 *
 */

#ifndef MORIS_IOS_CL_PERFORMANCE_REPORTER_HPP_
#define MORIS_IOS_CL_PERFORMANCE_REPORTER_HPP_

#include <chrono>
#include <map>
#include <string>

#include "moris_typedefs.hpp"

namespace moris
{
    //------------------------------------------------------------------------------

    /**
     * One accumulated performance record for a single tracer region, keyed by the
     * (entity, type, action) triple it was signed in with. Values are summed over
     * every visit to that region across the whole run (i.e. across optimization
     * iterations, since each iteration re-enters the same named regions).
     */
    struct Module_Performance_Record
    {
        std::string mEntity;
        std::string mType;
        std::string mAction;

        uint mDepth        = 0;      // indentation level when the region was first seen
        uint mVisitCount   = 0;      // number of sign-out events for this region
        real mWallTimeSum  = 0.0;    // wall time [s], summed over visits
        real mMemDeltaSum  = 0.0;    // memory change [MB], summed over visits
        real mMemDeltaPeak = 0.0;    // largest single-visit memory growth [MB]
    };

    //------------------------------------------------------------------------------

    /**
     * Consolidated end-of-run performance report. A single global instance
     * (gPerfReporter) is fed by Logger::sign_out for every tracer-bounded region;
     * at the end of the run finalize() aggregates across MPI ranks, prints a console
     * summary table, and writes a JSON report.
     *
     * The reporter is a passive consumer of the existing tracer events - it adds no
     * instrumentation of its own. Everything is a no-op until initialize() is called,
     * so test executables that link IOS but never start a run are unaffected.
     */
    class Performance_Reporter
    {
        //------------------------------------ PRIVATE -----------------------------

      private:
        bool        mEnabled  = false;
        sint        mLevel    = 1;                    // granularity, see initialize()
        std::string mJsonPath = "perf_report.json";
        std::string mMemSource = "proc_rss";          // "tcmalloc" or "proc_rss"

        real mPeakMemMb = 0.0;                         // per-rank high-water memory mark

        std::chrono::system_clock::time_point mRunWallStart;

        // accumulated records keyed by "entity|type|action" (std::map -> sorted,
        // deterministic iteration order across ranks for MPI aggregation)
        std::map< std::string, Module_Performance_Record > mRecords;

        // returns true if a region at the given depth is shown for the current level
        bool include_in_report( uint aDepth ) const;

        //------------------------------------ PUBLIC ------------------------------

      public:
        Performance_Reporter() = default;

        //--------------------------------------------------------------------------

        /**
         * Enable reporting and set the granularity level.
         *
         * Levels:
         *  0 - run total only (total wall time + peak memory)
         *  1 - top-level modules (HMR, XTK, GEN, MDL, ...), one row each
         *  2 - modules + nested sub-regions (e.g. MDL forward/sensitivity, SOL)
         *  3 - full nested tracer tree
         *
         * @param aLevel    granularity level (clamped to [0,3])
         * @param aJsonPath path of the JSON report written on rank 0
         */
        void initialize( sint aLevel, const std::string& aJsonPath );

        //--------------------------------------------------------------------------

        /**
         * Record a completed tracer region. Called from Logger::sign_out.
         *
         * @param aEntity      entity base (e.g. "HMR")
         * @param aType        entity type
         * @param aAction      entity action (e.g. "Perform Forward Analysis")
         * @param aDepth       indentation level of the region
         * @param aWallSeconds wall time spent in the region [s]
         * @param aMemDeltaMb  memory change across the region [MB]
         */
        void record_region(
                const std::string& aEntity,
                const std::string& aType,
                const std::string& aAction,
                uint               aDepth,
                real               aWallSeconds,
                real               aMemDeltaMb );

        //--------------------------------------------------------------------------

        // update the run-level peak (high-water) memory mark [MB]
        void update_peak_rss( real aMemMb );

        //--------------------------------------------------------------------------

        // aggregate across ranks, print the console table and write the JSON report
        void finalize();

        //--------------------------------------------------------------------------

        bool is_enabled() const { return mEnabled; }
        sint level() const { return mLevel; }

        //--------------------------------------------------------------------------

    };    // class Performance_Reporter

}    // namespace moris

// global performance reporter, defined in cl_Performance_Reporter.cpp
extern moris::Performance_Reporter gPerfReporter;

#endif /* MORIS_IOS_CL_PERFORMANCE_REPORTER_HPP_ */
