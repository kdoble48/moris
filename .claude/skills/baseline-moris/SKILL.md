---
name: baseline-moris
description: Use when starting a new MORIS feature, code change, or experiment — BEFORE editing any source. Symptoms — "add a parameter to MORIS", "implement/change the solver/IWG/GEN/optimizer", starting work on a MORIS branch, about to modify projects/, "let's add X to MORIS". Establishes a known-good reference (green build + tests AND a captured reference run) so a later regression is detectable and attributable to your change.
---

# Baseline MORIS Before a Feature

## Overview

Before changing MORIS, capture ground truth so that a later failure is unambiguously attributable to
*your* change and you have a concrete reference to diff against. The trap: **a green build is
necessary but NOT sufficient.** What actually regresses in MORIS is *run behavior* — a converged
objective, a design, a sensitivity — which the unit tests do not cover. And MORIS runs are slow
(minutes to hours), so the reference run must be captured **up front, before you touch code** — after
the fact you can no longer separate your change from the baseline.

Pairs with: `triage-moris-runs` (your baseline IS the reference that makes triage possible) and
`understand-moris-input-deck` (to read the deck you are baselining).

## The baseline checklist — do all of it before editing source

1. **Audit the tree state.** `git status` — what is already uncommitted? Know this so your later diff
   is clean and you do not attribute pre-existing changes to yourself. Do **not** delete or stash
   untracked work that isn't yours (generated `results/`, receipts, or someone's WIP). If not on a
   feature branch, branch first.

2. **Pin the reference point.** Record `git rev-parse --short HEAD` and the build date. For changes
   where "is it the deck or my code?" may come up later, keep a **pre-change binary** (copy
   `build_opt/projects/mains/moris`, or use a pinned clone) so you can reproduce-on-old to localize.

3. **Build green FIRST, on the current tree, before adding anything.** `cd build_opt && ninja -j2
   <target>` must be clean *now*. If it isn't, understand/fix that before your change — otherwise a
   later failure is ambiguous. Use **`-j2`**: the build is memory-hungry; higher `-j` exhausts RAM.

4. **Tests green.** Run the relevant suite and record the count:
   `./projects/<MOD>/test/bin/<MOD>-test.exe` (e.g. `OPT-test.exe`) or `ctest -L <LABEL>`.

5. **Capture a reference RUN — the critical, MORIS-specific step.** Pick the deck/case your change
   will affect and run it NOW, unchanged, saving into a `baseline/` dir:
   - `Parameter_Receipt.xml` — what MORIS actually used,
   - `moris_run.log` — objective/constraint trajectory + eval count,
   - the final-design exo (the VIS `<name>.exo.<np>.<r>` set) or a render.
   This is your falsifier: after the change you diff the new run against it. Without it you cannot tell
   whether your feature changed the answer.

6. **Fix the comparison protocol up front.** Record the **MPI rank count** you will standardize on —
   different `np` gives a different result (decomposition noise; see `triage-moris-runs`) — and run the
   after-change comparison at the *same* `np`, same deck, same mesh.

## Adding a tunable: the three-layer PRM pattern

Exposing a new MORIS parameter is always three edits + a test:
1. `projects/PRM/src/fn_PRM_<MODULE>_Parameters.hpp` — `tParameterList.insert("key", default)`. **The
   default must be current behavior (a no-op)** so the baseline run is reproduced bit-for-bit.
2. The consuming module — read it (`aParameterList.get<T>("key")`, guard with `.exists()` if optional)
   and act on it.
3. `projects/<MODULE>/test/` — a unit test covering the parameter on and off.

Mirror an existing precedent: `sdf_reinit_mode`, `reinit_in_place`, `grad_clip_factor`.

## Common mistakes

- **Green build mistaken for a baseline.** Unit-green ≠ behavior-captured. Run the reference deck
  *before* editing — afterwards you can't separate your change from the baseline. (This is the #1 gap.)
- **Starting on a dirty tree.** Uncommitted unrelated changes fold into your diff and can break your
  run (this has bitten — a pile of stray GEN/MTK edits crashed a run mid-triage). Audit first; build
  green with what's there before adding yours.
- **New-param default changes behavior.** The default must reproduce the baseline run exactly; verify
  the baseline deck output is unchanged after the (default-off) param is added.
- **No `np` recorded.** You later compare at a different rank count and chase decomposition noise.
- **Deleting untracked artifacts to "tidy up".** Generated results/receipts may be a reference. Leave
  them.
