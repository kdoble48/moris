# Cross-module couplings

A MORIS deck is string-coupled: a name defined in one module's parameter list is referenced, as
a string, in another. Nothing checks these references at compile time. Most deck bugs are broken
couplings rather than bad individual keys, and several of them produce wrong results without an
error.

This page lists the eight couplings (C1–C8) and the direction in which each must stay
consistent. The fork checks C1, C4, C5, and C8 at load time (see
[Load-time checks](#load-time-checks)).

## C1 — Criteria order: GEN `IQI_types` → `aCriteria` → OPT callbacks

The **order** of names in GEN's `IQI_types` is the index order of the `aCriteria` vector passed
to `compute_objectives` and `compute_constraints`: `aCriteria(0)` is the first name,
`aCriteria(1)` the second. Reordering `IQI_types` without re-indexing the callbacks changes the
objective silently.

- Every name in `IQI_types` must be an `IQI_name` defined in FEM.
- The callbacks must index `aCriteria` by the position of each name.
- `compute_d*_dcriteria` returns the chain-rule terms of the objective and constraints with
  respect to the criteria; `compute_d*_dadv` is usually zero, because the ADV sensitivities flow
  through the FEM adjoint.

The {doc}`Deck API v2 <deck-api-v2>` removes this coupling: `IQI_types` is generated from the
objective and constraint expressions.

## C2 — Problem declaration: OPT

`is_optimization_problem = true` with `problem = "user_defined"` means the objective and
constraints come from the deck. The algorithm decides what kind of run it is:
`GCMMA`, `MMA`, `ROL`, `LBFGS`, and `SQP` run a design optimization; **`SWEEP` runs a
finite-difference sensitivity check**, not an optimization.

## C3 — How ADVs reach the physics: GEN → FEM

Design variables reach the physics by one of two pathways:

a. **Direct geometry.** A GEN level-set geometry with `discretization_mesh_index` set: its
   B-spline coefficients are the ADVs, and they move the interface.
b. **PDV-coupled property.** A GEN property (often `SCALED_FIELD`) depends on an ADV field and
   publishes a `pdv_type`; FEM properties consume it through `dv_dependencies` and
   `dv_derivative_functions`.

If a deck has design variables but no discretized geometry, the property → `pdv_type` →
`dv_dependencies` chain must be complete at both ends.

## C4 — Inside FEM: bottom-up references

FEM sub-lists are filled in order — properties, constitutive models, stabilization parameters,
IWGs, IQIs — and each later entry references earlier ones by name, paired with a role:

```cpp
// a constitutive model lists its properties:  <your name>,<role the CM expects>
tFem.set( "properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio" );

// an IWG lists its constitutive models and stabilization parameters
tFem.set( "leader_constitutive_models", "CMStrucLinIso,ElastLinIso" );
tFem.set( "stabilization_parameters",  "SPNitscheDirichlet,DirichletNitsche" );
```

The left-hand name is yours; the right-hand role is fixed by the consuming class. An `IQI_name`
is an arbitrary label; its `IQI_type` selects the computation. A dangling name on the left of a
pairing is either a load-time abort or silently missing physics.

## C5 — Mesh-set names: XTK produces, FEM consumes

XTK generates set names from the phase index (`pN`, from GEN's `phase_table`) and the cut
status. FEM's `mesh_set_names` must match them exactly:

| Set name | Meaning |
|---|---|
| `HMR_dummy_c_p2`, `HMR_dummy_n_p2` | bulk cells of phase 2 — `c` cut by the interface, `n` not cut |
| `SideSet_4_n_p2` | external boundary side set 4 (HMR numbering), phase 2 |
| `iside_b0_2_b1_0` | single-sided interface between phase 2 and phase 0 |
| `dbl_iside_p0_1_p1_2` | double-sided interface between phases 1 and 2 (Nitsche coupling) |
| `ghost_p1` | ghost-penalty facets of phase 1 — exist only when `ghost_stab` is on |

Every phase index referenced by FEM must exist in the phase table, and interface sets must name
real phase pairs. The Deck API v2 `deck::Phase` vocabulary generates these strings.

## C6 — Mesh indices: HMR → XTK → GEN

HMR can define several refinement patterns and B-spline meshes. Their indices thread through:
XTK `enrich_mesh_indices` selects the enriched mesh; GEN `discretization_mesh_index` and
`refinement_mesh_index` select the mesh a field is stored on and refined by. An out-of-range
index is a common silent failure.

## C7 — DOF types: FEM ↔ MSI ↔ SOL

The DOF types used by FEM IWGs (`UX,UY`, `TEMP`, `P`, …) must be enabled in MSI and listed in
SOL's `NLA_DofTypes` and `TSA_DofTypes`. Staged multiphysics solves split the DOFs across
NLBGS sub-solvers. A missing DOF type drops equations or stops convergence.

## C8 — Output: VIS → FEM and XTK

VIS `IQI_Names` must be FEM IQIs, and VIS `Set_Names` must be XTK-generated sets. A mismatch
loses output rather than corrupting results.

## Load-time checks

When a deck loads, the fork runs a consistency pass over C1, C4, C5, and C8 and reports
violations with "did you mean" suggestions. Violations are warnings by default; set
`MORIS_STRICT_INPUT=1` in the environment to make them errors.

## Checklist

| # | Check |
|---|---|
| C1 | `IQI_types` ⊆ FEM `IQI_name`s; the callbacks index `aCriteria` in `IQI_types` order |
| C2 | the algorithm matches the intent (`SWEEP` is a gradient check) |
| C3 | discretized geometry, or a complete `pdv_type` ↔ `dv_dependencies` chain |
| C4 | every IWG/IQI/CM name pairing resolves to a defined property, CM, or SP |
| C5 | `mesh_set_names` follow the XTK grammar; phases exist; ghost sets only with ghost on |
| C6 | HMR mesh indices ↔ XTK `enrich_mesh_indices` ↔ GEN mesh indices |
| C7 | FEM DOFs = MSI DOFs = SOL `NLA_DofTypes` |
| C8 | VIS `IQI_Names` ⊆ FEM IQIs; `Set_Names` follow the XTK grammar |
