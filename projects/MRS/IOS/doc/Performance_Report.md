# Run Performance Report

A consolidated end-of-run report of **where a MORIS run spent its time and memory**, broken
down per module. It reuses the existing tracer instrumentation (`Tracer` / `GlobalClock` /
`Logger`), so no new instrumentation is added at call sites — the report is a passive consumer of
the sign-in / sign-out events MORIS already emits at every module boundary (HMR, XTK, GEN, MDL,
SOL, …).

Two artifacts are produced at the end of a run:

1. A **console summary table** (printed via the logger, rank 0).
2. A **JSON report** written by rank 0 to `performance_report_file` (default `perf_report.json`),
   relative to the working directory the run was launched from.

## Quick start

```bash
# bare moris binary
moris input.so --perf-level 2                     # per-module table + JSON
moris input.so --perf-level 3 --perf-report run.json

# via an environment variable (works for any entry point, e.g. EXA example executables
# that do not forward CLI flags)
MORIS_PERF_LEVEL=3 ./Main_Two_Bar_Truss.exe
```

## Configuration

The granularity level and output path are resolved in this order (last wins):

1. OPT parameter list: `performance_report_level` (int, default `1`) and
   `performance_report_file` (string, default `"perf_report.json"`), declared in
   `prm::create_opt_problem_parameter_list()`.
2. Environment variable `MORIS_PERF_LEVEL`.
3. CLI flags `--perf-level N` and `--perf-report <path>`.

A negative level disables the report. Decks that predate these parameters still run unchanged —
the reader guards every lookup with `exists()`.

### Granularity levels

| Level | Contents |
|-------|----------|
| `0` | Run total only — total wall time + peak memory. |
| `1` | Top-level modules (depth ≤ 2): HMR, XTK, GEN, MDL, … one row each. |
| `2` | Modules + nested sub-regions (depth ≤ 4): MDL forward/sensitivity, the solver stack, … |
| `3` | Full nested tracer tree (every region), plus an optional gperftools CPU callgrind profile. |

## Output

### Console table

```
==================== MORIS Performance Report (level 1) ====================
Module (entity - action)               Visits     Wall[s]      Wall max/min     MemD[MB]    MemDpeak
-------------------------------------------------------------------------------------------------------
OPT - Perform                               1       1.341         1.34/1.34          4.0         4.0
MDL - Perform Forward Analysis              2       0.003         0.00/0.00          0.2         0.2
...
-------------------------------------------------------------------------------------------------------
Run total wall (max over ranks): 1.349 s   |   Peak memory (max over ranks): 12.8 MB [tcmalloc]
```

Rows are summed over every visit to a region (so timings accumulate across optimization
iterations) and sorted by wall time. `Wall max/min` are reductions across MPI ranks. The
`MemD[MB]` (growth) column is clamped at zero; the raw signed delta is kept in the JSON.

### JSON schema

```json
{
  "level": 1,
  "run_total_wall_s": 1.349,
  "peak_mem_mb": 12.82,
  "mem_source": "tcmalloc",        // or "proc_rss"
  "ranks_aligned": true,
  "modules": [
    {
      "entity": "OPT", "type": "Manager", "action": "Perform",
      "depth": 1, "visits": 1,
      "wall_s": 1.341, "wall_s_max": 1.341, "wall_s_min": 1.341,
      "mem_delta_mb": 4.03, "mem_delta_peak_mb": 4.03
    }
  ]
}
```

## Memory source

`Logger::current_memory_mb()` provides the numbers:

- **With gperftools** (`MORIS_USE_GPERFTOOLS=ON`): tcmalloc's
  `generic.current_allocated_bytes` via `MallocExtension`. `mem_source` reports `tcmalloc`.
- **Otherwise**: physical resident set size from `/proc/self/statm`. `mem_source` reports
  `proc_rss`.

Per-module values are the **change** in this metric across a region; the run-level peak is the
high-water mark sampled at every boundary. RSS deltas can be negative (allocator page returns) —
that is why the console growth column is clamped while the JSON keeps the raw value.

## gperftools

Building with `-DMORIS_USE_GPERFTOOLS=ON` (which defines `WITHGPERFTOOLS`, links `-ltcmalloc
-lprofiler`, and needs the system `libgoogle-perftools-dev` package) enables two things:

- tcmalloc-sourced memory numbers (above).
- At `--perf-level 3`, the whole run is wrapped in `moris::Profiler`, which emits a callgrind CPU
  profile via `pprof` (falls back to the system `google-pprof` when `$APPS/gperftools/bin/pprof`
  is absent).

Without gperftools, the report still produces the full time + memory summary via the `/proc`
fallback — only the deep CPU profile is unavailable.

## Implementation

| File | Role |
|------|------|
| `cl_Performance_Reporter.{hpp,cpp}` | Accumulator + MPI aggregation + console table + JSON writer; global `gPerfReporter`. |
| `cl_Logger.cpp` | `current_memory_mb()`; `sign_in`/`sign_out` feed `gPerfReporter.record_region` / `update_peak_rss`. |
| `cl_GlobalClock.{hpp,cpp}` | `mMemoryStamps` — per-region sign-in memory snapshot. |
| `fn_WRK_Workflow_Main_Interface.cpp` | Resolves the level/path config and calls `initialize()`; optional gperftools CPU profile. |
| `main.cpp`, `EXA/src/example_main.cpp` | Call `gPerfReporter.finalize()` before MPI finalize. |
| `fn_PRM_OPT_Parameters.hpp` | `performance_report_level` / `performance_report_file` defaults. |

Unit test: `projects/MRS/IOS/test/UT_IOS_Performance_Reporter.cpp`.
