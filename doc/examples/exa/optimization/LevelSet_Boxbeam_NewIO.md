# Level-set box beam (Deck API v2)

A two-phase level-set topology optimization of a 2D cantilever box beam: a fixed outer frame
surrounds a design region whose material layout is optimized for stiffness under a mass
constraint. The example is the {doc}`Deck API v2 <../../../user-guide/deck-api-v2>` twin of
{doc}`LevelSet_Boxbeam <LevelSet_Boxbeam>`: the same problem written with a single entry point,
typed set names, a physics preset, and objective and constraint expressions instead of OPT
callbacks.

```{exa-summary} projects/EXA/optimization/LevelSet_Boxbeam_NewIO
```

## Problem

| | |
|---|---|
| Domain | 3.0 × 1.0 rectangle |
| Phases | 0 void, 1 design interior, 2 frame |
| Material | linear elastic, isotropic (preset defaults: E = 1.0, ν = 0.0), plane problem in UX, UY |
| Supports | left edge of the frame, Nitsche-enforced Dirichlet condition |
| Load | traction (0, 1.0) on the right edge of the frame, applied only below y = `tLoadLimitY` (≈ 0.2) |
| Objective | strain energy of frame and interior, plus 0.2 × the material–void interface perimeter, each normalized by its initial value |
| Constraint | volume of the interior phase ≤ 1.0 |
| Optimizer | GCMMA, step size 0.01, penalty 5.0, 2 iterations by default |

The frame is a superellipse (exponent 24) inset by a wall thickness of 0.05; the design region
is seeded with a 15 × 5 array of small superellipse holes (semi-diameter 0.065, exponent 8).

## Setup in MORIS

**Mesh (HMR).** 45 × 15 elements with two levels of initial refinement (element edge 1/60).
One linear Lagrange mesh and two linear B-spline meshes: mesh 0 carries the level-set ADVs,
mesh 1 the displacement DOFs (`tMsi.set( "UX", 1 )`).

**Geometry and phases (GEN).** Two level-set geometries — the frame superellipse and the
hole field array. The array's B-spline coefficients on mesh 0 are the design variables
(`discretization_mesh_index = 0`), bounded to ±0.06. The phase table `[0, 1, 2, 2]` maps the two
geometry signs to void, interior, and frame. `IQI_types` is not written by hand: the expression
API derives it from the objective and constraint
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of meshes 0 and 1, ghost
stabilization on both material phases, and Exodus output of the integration mesh.

**Physics (FEM).** `fem::presets::linear_elastic()` emits one complete linear-elastic block per
material phase (prefixes `Frame_` and `Interior_`): properties, constitutive model, bulk IWG,
ghost penalty, a bedding stiffness of 10⁻⁶, and strain-energy and volume IQIs. The frame block
adds the Nitsche support and the traction, whose spatial limit is the one custom callback
(`Func_Neumann_U`, registered in-process). Raw parameter lists add the frame–interior Nitsche
interface (`STRUC_LINEAR_INTERFACE_UNSYMMETRIC_NITSCHE`, penalty 100), the perimeter IQI
(`VOLUME` integrated over the material–void interface sets), and displacement IQIs for output.
All set names come from `deck::Phase` handles
([coupling C5](../../../user-guide/interlinks.md#c5-mesh-set-names-xtk-produces-fem-consumes)).

**Solver (SOL).** Direct solve with `Amesos_Klu` (set explicitly; the Amesos default, Pardiso,
is absent from many Trilinos builds), one Newton iteration, one static time step.

**Optimization (OPT).** GCMMA on the expressions passed to `aDeck.objective()` and
`aDeck.constraint()`; their criteria gradients come from reverse-mode differentiation.

**Output (VIS).** `LevelSet_Boxbeam_NewIO.exo` with nodal UX and UY on the frame and interior
sets, written every iteration.

## What the test checks

The ctest runs the full optimization through the MORIS workflow interface and requires a
return code of 0. It checks that the deck loads, cuts, solves, and completes its GCMMA
iterations; it does not compare objective values.

## Run it

```bash
ctest -R LevelSet_Boxbeam_NewIO --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

Two environment variables override the defaults for longer runs:
`BOXBEAM_MAX_ITS` (GCMMA iterations, default 2) and `BOXBEAM_STEP_SIZE` (default 0.01).
