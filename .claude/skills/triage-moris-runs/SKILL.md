---
name: triage-moris-runs
description: Use when diagnosing a MORIS optimization/simulation run that crashed, aborted, segfaulted, hit a NaN, stalled or diverged, converged to an unexpected design, or differs from another run that "should be the same". Symptoms — "why did this run X", "these two runs give different results", gradient explosion, ADV-consistency assert, "does not match", restart loop, Belos not converged, optimizer froze/lattice instead of truss, comparing two Parameter_Receipt.xml.
---

# Triage a MORIS Run

## Overview

Triaging a MORIS run is **differential debugging**: you explain a run by contrast with a reference
where the behavior differs. Two rules dominate, because violating either produces a confident wrong
answer:

1. **The deck source lies; the receipt and the log are ground truth.** What a `.cpp`/`.py` says is
   not what MORIS ran — pymoris codegen has silently dropped keys (e.g. `step_size`/`penalty`), and a
   `.so` inlines its own parameter lists at *its* compile time. Read `Parameter_Receipt.xml` (what
   MORIS actually used) and `moris_run.log` (what actually happened).
2. **A plausible mechanism is a hypothesis, not a conclusion.** The decisive test is a controlled
   re-run with exactly one knob changed — not a story about why a parameter "should" matter. Run it.

**REQUIRED BACKGROUND for deck structure:** use the `understand-moris-input-deck` skill to read any
deck (criteria ordering, GEN↔FEM string links, mesh-set names).

## The five disciplines

**D1 — Anchor on a falsifier / reference.** If something "should" work, find where it *did* (a
different resolution, an older pinned binary, a reference deck, an earlier commit) — that's your
reference. If none exists, *create* one by toggling to a known-good configuration. With no reference
you are guessing, not triaging.

**D2 — Diff the receipt, section by section, not the deck.** `Parameter_Receipt.xml` has one block per
module: `OPT HMR XTK GEN FEM SOL MSI VIS MORISGENERAL`. Extract and diff each:
```bash
ex(){ awk "/^\t<$2>/,/^\t<\/$2>/" "$1/Parameter_Receipt.xml" | sed 's/^[[:space:]]*//; s/&quot;/"/g'; }
diff <(ex runA HMR) <(ex runB HMR)        # mesh identical? (number_of_elements_per_dimension, orders)
diff <(ex runA OPT) <(ex runB OPT)        # optimizer + algorithm params identical?
```
The receipt also captures what the deck *omitted* (it shows the resolved defaults) — that is how you
catch silent codegen gaps.

**D3 — Change exactly ONE variable.** Hold mesh, MPI rank count, optimizer, and deck constant except
the one knob under test. A comparison that differs in two things (e.g. `np` *and* a parameter) cannot
attribute the result to either. Most "mysterious" run differences are confounded comparisons.

**D4 — Verify a flagged parameter is actually LIVE before blaming it.** A diff is not a cause. Check
the parameter is in force:
- A GEN geometry with `discretization_mesh_index < 0` is **not a design variable** — its
  `discretization_lower/upper_bound` are inert (no ADVs created for it). Only geometries with
  `discretization_mesh_index = 0` (or a real B-spline index) are optimized.
- A FEM `Property`/`IWG`/`CM`/`SP` defined but referenced by nothing is **dead** — `grep` its name; one
  hit (its own definition) = unused.
- User-function *bodies* (geometry fields, load, `compute_objectives`) are **not** in the receipt —
  `diff` them directly between decks before concluding the decks match.

**D5 — Confirm by controlled re-run, never decline it.** Re-run with one knob flipped, or at matched
config. **Bit-for-bit reproduction proves equivalence** (two MMA runs cannot share a 15-digit,
all-iterations-identical objective trajectory unless every live input matches). "No need to re-run, the
mechanism is clear" is the signature of a wrong triage.

## MPI rank count is a real variable

Different `np` → different XTK domain decomposition → different floating-point reduction order →
**different local minimum** in a non-convex TO problem. Expect a run-to-run objective spread of a
percent or so purely from `np`. When comparing runs, hold `np` constant (check `Procs Used` in each
log header) or you will misattribute decomposition noise to the deck.

## Reading the log (grep the diagnostic stream)

```bash
grep -c "Objective:" log            # eval count; tail | grep -v nan for the converged value
grep "GradDebug" log | tail         # obj_norm/obj_max (gradient health), at_lower/at_upper (bound saturation), n_advs (DOF)
grep -c "Number of optimization variables" log   # restart/re-init count; the value = DOF (changes on remesh)
grep -cE "objective scale ->|Gradient clip|non-finite design" log  # which robustness path engaged
grep -nE "Reason:|Segmentation|out_of_range|Aborted" log   # fatal cause
grep "MinADV: -nan" log             # NaN-poisoned design vector (overflow upstream)
grep "ElapsedWallTime =" log        # present == clean finish
```
Healthy XFEM design gradients are O(0.1); `obj_max` to ~1e8 is tolerated raw, ~1e10 overflows the
GCMMA subproblem to NaN. A self-contradictory line (`obj_norm=1e8` with `obj_max=1e-38`) means a clip
fired between the two measurements.

## Band-aid vs root cause

A run kept alive by a gradient clip, every-step reinit, or restart-on-NaN may be **masking an
ill-conditioned setup**, not solving it. Localize the cause:
- Reproduce on a **pinned/older binary** (e.g. `moris_8a328272e`) → if it also fails, the cause is the
  deck/conditioning, not your code change. Run a deck against another binary's headers by compiling with
  that `MORISROOT`.
- Toggle the band-aid off and fix conditioning instead (tight `discretization` bounds ≈ one element
  diagonal, adaptive remeshing) — often the band-aid becomes unnecessary.

## Run mechanics

```bash
spack env activate /home/doble/codes && source ~/.bashrc_moris
# direct binary also needs MKL on the path (moris-run injects it automatically):
export MKL_DIR=/home/doble/codes/spack/opt/spack/linux-x86_64_v4/intel-oneapi-mkl-2024.2.2-4xeqh62vh6x3j5xlqhlyuawrgkwdquib
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_DIR/mkl/2024.2/lib
moris-run deck.cpp --compile-only -o .                       # deck.cpp -> deck.so (receipt-faithful)
mpirun -np 8 $MORISROOT/build_opt/projects/mains/moris ./deck.so > moris_run.log 2>&1
```
- The moris-run launch shell has `set -e`: a leading `pkill`/`pgrep` returning 1 aborts the whole
  compound command before later lines — guard with `|| true` or run as separate commands.
- Run output goes under `$MORIS_RUNS_DIR` (`runs/`), never inside a git repo or `studies/`.
- Render the **final design** from the VIS exo (`<name>.exo.<np>.<r>`); the `GEN_*.e-s.*` files are
  per-iteration *evolution* snapshots (animation only, huge — prune them to reclaim space).

## Common mistakes

- **Confounded comparison** — comparing two runs that differ in `np` *and* a parameter, then blaming
  the parameter. Match everything but one knob (D3).
- **Blaming an inert parameter** — flagging `discretization_*_bound = ±1` on a geometry whose
  `discretization_mesh_index = -2` (not a design var), or a dead `PropFlux`. Verify it's live (D4).
- **Mechanism instead of experiment** — building a regularization story and declining the one-knob
  re-run that would confirm or kill it (D5).
- **Trusting the deck source** — concluding from `.cpp`/`.py` instead of the receipt; the run name
  (e.g. `step001`) may not reflect what ran (it inherited the default 0.01, never set it).
- **Reading the deck linearly** — use `understand-moris-input-deck`; criteria order is GEN's
  `IQI_types`, not FEM IQI order.
- **Treating NaN/`obj_max=1e10` as the disease** — it is a symptom of small-cut-cell ill-conditioning
  upstream; chase the conditioning (bounds, remeshing), not the explosion.
