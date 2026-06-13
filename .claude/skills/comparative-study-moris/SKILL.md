---
name: comparative-study-moris
description: Use when running a parametric or comparison study across MORIS runs — sweeping a parameter, mesh resolution, optimizer, or deck variant to quantify its effect on the result. Symptoms — "compare step_size 0.01 vs 0.05", "how does mesh resolution affect the design", "sweep parameter P", "which optimizer is better", "study the effect of X", building a run matrix, comparing converged objectives/designs across runs.
---

# Comparative Study Across MORIS Runs

## Overview

A comparative study isolates the effect of **one** factor on a MORIS result by running a matrix where
everything else is held constant. The entire value is in the controls — a sloppy study yields a
confident but wrong verdict. Two failure modes dominate:

1. **Comparing non-converged runs.** Truncating every arm at a short fixed iteration budget and
   declaring a winner measures *transient descent speed*, not the *optimum*. A method that moves faster
   early can look "better" at iter 50 yet converge to the same or a worse design. Either run to
   convergence, or state explicitly that you are comparing transient behavior.
2. **Confounded or noise-dominated comparison.** More than one thing changed between arms, or the
   measured effect is smaller than the run-to-run noise floor.

Pairs with: `baseline-moris` (capture the reference first) and `triage-moris-runs` (when one arm
differs unexpectedly — the verify/receipt-diff disciplines there apply here too).

## Design the matrix — one factor at a time

- Start every arm from **one common ancestor deck**. Verify it: `sha256sum`/`diff` must show that arms
  differ in **exactly** the intended line(s). Vary one knob; hold mesh, `np`, optimizer, physics,
  seeding, and iteration budget fixed.
- Record the matrix explicitly: a table of {arm → varied value | everything-else fixed}. One dir per
  arm under `runs/studies/<study>/`.
- **Hold `np` constant across arms.** Different rank counts give a different decomposition → ~1 %
  objective spread (the noise floor; see `triage-moris-runs`). If the effect you're measuring is below
  that, it's noise — enlarge the effect, or run replicates.

## Run it

```bash
moris-run arm.cpp --compile-only -o .
mpirun -np N $MORISROOT/build_opt/projects/mains/moris ./arm.so > moris_run.log 2>&1
```
- **Run to convergence** — `max_its` large enough that the design stops changing and KKT `norm_drop`
  is hit — not a short fixed budget, unless transient comparison is the explicit goal. Confirm each log
  ended with `ElapsedWallTime` (clean finish) and the objective plateaued (`grep Objective: log | tail`).
- If you must modify the deck to make it run (library name, solver backend, iteration count), apply the
  **same** change to **every** arm, record it, and verify the modified deck still reproduces a
  known-good baseline before trusting the study — a workaround can interact with the factor under test.

## Verify each arm used the intended settings — from the receipt, not the source

- `Parameter_Receipt.xml`: the varied key holds the intended value in each arm.
- **Receipt-diff the arms section by section** (`awk '/^\t<SEC>/,/^\t<\/SEC>/'` over OPT/HMR/XTK/GEN/
  FEM/SOL): confirm **exactly** the intended parameter differs and nothing else moved. A diff that is
  inert (a `discretization` bound on a `discretization_mesh_index < 0` geometry, a dead property) is
  fine — but you must recognize it as inert (see `triage-moris-runs` D4), not ignore it.
- Cross-check the algorithm's runtime echo in the log if it prints its config.

## Present the comparison

- **Numbers:** one table — final objective, constraint, eval count, wall time per arm — and **state
  whether each converged**.
- **Designs:** render every arm with the **same** script, same color / scale / extent, and compose into
  one figure (side-by-side or grid). Render the final-design VIS exo (`<name>.exo.<np>.<r>`), not the
  `GEN_*.e-s.*` evolution snapshots. Identical rendering is what makes the visual comparison honest —
  different scales/crops make designs look more different than they are.
- **Verdict with caveats:** name the noise floor, say whether the comparison is at convergence or
  transient, and what would strengthen it (longer budget, replicates, a third level).

## Common mistakes

- **Winner from truncated runs.** "0.05 beats 0.01 at iter 50" ≠ "0.05 is better" if neither converged.
  Run to convergence or label the result transient.
- **Confounded matrix.** Two knobs changed (e.g. `step_size` *and* `np`, or the factor *and* a deck
  workaround). Receipt-diff to prove one-factor before reporting.
- **Effect inside the noise floor.** A sub-1 % objective gap at different `np` is decomposition noise,
  not a finding.
- **Inconsistent rendering.** Different color scales/crops exaggerate differences. One script, one
  scale, composed panels.
- **Trusting deck source over the receipt.** Verify the live value per arm; a run/dir name (e.g.
  `step001`) may not reflect what actually ran.
