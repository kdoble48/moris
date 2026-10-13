# From a new deck to a trustworthy run

An optimization run that fails late is expensive, and many failures are detectable within the
first few iterations. Work through these stages in order.

## 1. Start from an example

Find an example in `projects/EXA/` with the features you need — the same physics, the same kind
of geometry, the same optimizer — and adapt it. The example decks are run by the test suite, so
they are known to load and run.

## 2. Forward analysis

Run the deck as an analysis (`is_optimization_problem = false`) and inspect the output in
ParaView:

- Do the cuts and the phase assignment look right? Write `xtk_temp.exo`
  (`exodus_output_XTK_ig_mesh`) and color by phase.
- Are the boundary conditions on the intended sets? A Dirichlet condition on a set that does not
  exist does nothing.
- Is the solution qualitatively right — load path, deflection direction, magnitudes?

## 3. A few optimization iterations

Run two or three iterations (`max_its`). This catches crashes in the sensitivity analysis before
any time is spent on a sweep. Read the first iteration's criteria printout: are the objective
terms of the expected magnitude relative to their normalizers, and is the initial design near the
constraint?

## 4. Sensitivity sweep

Verify the analytic gradients against finite differences with the `SWEEP` algorithm:

```cpp
tOpt( OPT::ALGORITHMS ).add_parameter_list( opt::Optimization_Algorithm_Type::SWEEP );
tOpt.set( "hdf5_path", "sweep.hdf5" );
tOpt.set( "num_evaluations_per_adv", "1" );
tOpt.set( "finite_difference_type", "all" );
tOpt.set( "finite_difference_epsilons", "1e-6" );
tOpt.set( "finite_difference_adv_indices", "0,5,17" );   // optional subset of ADVs
```

The sweep writes the analytic and finite-difference objective and constraint gradients for each
ADV to the HDF5 file; compare them per ADV. Pick ADVs near the interface: ADVs far from it have
zero sensitivity and agree trivially. Any HDF5 reader works for the comparison (`h5py`,
HDFView).

## 5. Production run

Launch the full run and watch the first iterations:

- criteria-gradient norms finite (a `nan` never recovers);
- objective changing between iterations (a frozen objective means the design is not moving);
- no sustained explosion of the objective.

If a run must stop, restart it from the last good `ADV_Alg_0_Iter_N.hdf5` via `restart_file`
instead of starting over.

## Tuning GCMMA

| Parameter | Effect |
|---|---|
| `step_size` | the move limit: how far one iteration may move the design |
| `penalty` | how strongly constraint violation is penalized |
| `norm_drop` | the relative-change convergence criterion |
| `max_its` | iteration budget |

A run that oscillates usually needs a smaller `step_size`. A run that stays infeasible needs a
larger `penalty`, or a design start closer to the constraint.

## Level-set initialization

Holes in a pure level-set design can merge and split, but they cannot nucleate. The initial hole
pattern therefore limits the reachable topologies:

- Seed densely (a field array of small holes) rather than sparsely.
- Start near the constraint: an infeasible start forces a feasibility march that tends to merge
  holes early.
- Keep perimeter regularization small (≲ 1% of the objective) while the topology forms.
- If two different starts converge to the same design, that design is the optimum of the
  formulation; change the formulation, not the seed.
