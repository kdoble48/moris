---
name: understand-moris-input-deck
description: Use when reading, explaining, debugging, modifying, or porting a MORIS C++ input deck — an EXA example or any .cpp that defines OPTParameterList/HMRParameterList/GENParameterList/FEMParameterList etc. Symptoms: "what does this input file do", "where is the objective/geometry/BC defined", unfamiliar aParameterLists.set keys, string names that match across blocks, set names like HMR_dummy_c_p2 / iside_ / ghost_p1.
---

# Understand a MORIS Input Deck

## Overview

A MORIS input deck is **one C++ file** (in `namespace moris`, wrapped in `extern "C"`) that fully
configures a simulation. It has no MORIS logic — only data: top-of-file constants, one
`*ParameterList(Module_Parameter_Lists&)` function per module, and user callbacks. Everything is
**string-coupled**: names defined in one block are referenced by string in others, with no compiler
check. Read by *tracing those links*, not top-to-bottom.

Don't go spelunking through factories/enums to answer structural questions — the structure below is
fixed across every deck. Only open other files to **decode a specific key** (see last section).

## Source-of-truth documentation

Before reverse-engineering, check the authoritative docs (all under `moris/share/doc/`):

- **`Using_Moris.dox`** (`@page UsingMoris`) — the canonical guide to *writing* input decks: input-file
  steps, recommended HMR/SOL settings, convergence tweaks, and the sweep sensitivity-check workflow.
  Start here for "how should this be set up / why."
- **`Wrk_HMR_XTK_GEN_FEM_Workflow.dox`** — the module workflow, and a table mapping each
  `*ParameterList` function to its `fn_PRM_<MODULE>_Parameters.hpp` header (the per-key decode source).
  Companion docs `Wrk_STK_XTK_WorkFlow.dox`, `Wrk_STK_FEM_WorkFlow.dox` cover the non-XFEM workflows.
- **`mesh_generation/main.pdf`** (+ `mesh_generation/examples/`) — the XML-driven mesh-generation
  input format (separate from the optimization decks here).
- **`projects/PRM/src/fn_PRM_<MODULE>_Parameters.hpp`** — the authoritative, code-level list of every key
  and its default (see last section).
- **The EXA decks themselves** — `Using_Moris.dox` explicitly says to start from an existing example with
  the features you want and adapt it; `projects/EXA/*/` is the de-facto reference library.

These are reference, not gospel: defaults and enums live in code, so when a doc and the `fn_PRM_*`/enum
source disagree, trust the source.

## Read order

1. **Top-of-file constants** (`tName`, `tIsOpt`, `tMaxMass`, set-name strings…) — these parameterize
   everything below. The "Derived parameters" section assembles mesh-set name strings.
2. **`OPTParameterList`** — is `is_optimization_problem` true? `problem = "user_defined"` means the
   objective/constraints are the C++ callbacks in this file. Check **which** algorithm `OPT::ALGORITHMS`
   adds: `Optimization_Algorithm_Type::GCMMA`/`SQP`/`LBFGS` is a real design run, but `::SWEEP` is a
   **finite-difference sensitivity check** (keys `num_evaluations_per_adv`, `finite_difference_type`,
   `finite_difference_adv_indices`, `hdf5_path`) — not optimization, even though `is_optimization_problem`
   is true. Don't equate "optimization problem = design run."
3. **`GENParameterList`** — geometry + design variables. `GEN::GEOMETRIES` entries define level-set
   primitives (`gen::Field_Type::SUPERELLIPSE`, `field_array` for tiled holes, or `USER_DEFINED` +
   `field_function_name`). The `discretization_*` keys (set only when `tIsOpt`) make a field's B-spline
   coefficients the **ADVs**, bounded by `discretization_lower/upper_bound`. `phase_table` maps geometry
   sign combinations → material phase index. Two design-variable pathways reach the physics — contract C.
4. **`FEMParameterList`** — the physics (see contract B below).
5. **`SOL`/`MSI`/`VIS`** — solver type, DOF order, Exodus output. Multiphysics decks stage sub-solvers:
   several `NONLINEAR_SOLVERS`/`NONLINEAR_ALGORITHMS` coupled by NLBGS, `TSA_Nonlinear_Solver` selecting
   the driver, `TIME_CONTINUITY_DOF` IWGs for transient fields. Don't assume a single solve.

## The two non-obvious contracts

**A. Criteria ordering links GEN → the OPT callbacks.** The order of names in GEN's
`set("IQI_types", ...)` **is** the index order of the `aCriteria` vector passed to
`compute_objectives`/`compute_constraints`. In `Levelset_Boxbeam.cpp`, `IQI_types` =
`[StrainEnergy_Frame, StrainEnergy_Interior, Volume_Interior, Perimeter_Void]`, so `aCriteria(0)` is
frame strain energy, `aCriteria(2)` is interior volume, etc. The objective/constraint callbacks index
this vector by position — get the order wrong and the math is silently wrong. `compute_d*_dcriteria`
returns the analytic chain-rule terms; `compute_d*_dadv` usually returns **zeros** because ADV
sensitivities flow through the FEM adjoint via the criteria.

**B. The FEM block is built bottom-up and linked by string.** Sub-lists are added in this order via
`aParameterLists( FEM::<SECTION> ).add_parameter_list()`:
`PROPERTIES → CONSTITUTIVE_MODELS → STABILIZATION → IWG → IQI → COMPUTATION`. Each later layer
references earlier ones **by name, paired with a role tag**:
- a CM lists its properties: `"properties", "PropYoungs,YoungsModulus;PropPoisson,PoissonRatio"`
- an IWG lists its CM/SP/properties: `"leader_constitutive_models", "CMStrucLinIso_Frame,ElastLinIso"`,
  `"stabilization_parameters", "SPNitscheDirichletBC,DirichletNitsche"`
- an IQI's `IQI_name` is a user label (referenced by GEN's `IQI_types`); its `IQI_type` enum selects
  the actual computation. The left side of each pair is the name you defined; the right side is the
  fixed role the consumer expects. `hack_for_legacy_fem()` at the top of the block is an ordering quirk,
  not meaningful config.

Properties hold their value in `function_parameters` and a `value_function` (a C++ callback in this
file, e.g. `Func_Const` returns the constant; `Func_Neumann_U` computes a position-dependent traction).

**C. ADVs reach the physics by one of two pathways.** (1) *Direct geometry:* a `GEN::GEOMETRIES`
level-set field with `discretization_*` set — the field's B-spline coefficients are the ADVs and they
move the interface (the Boxbeam case). (2) *PDV-coupled property:* a `GEN::PROPERTIES` entry (often
`gen::Field_Type::SCALED_FIELD`) that `dependencies` on an ADV field and publishes a `pdv_type` (e.g.
`LS1`); FEM properties then consume that PDV via `value_function` **plus** `dv_dependencies` and
`dv_derivative_functions` (which supply the design sensitivities). When a deck's design machinery seems
missing from `GEN::GEOMETRIES`, look for this `PROPERTIES`→`pdv_type`→FEM-`dv_dependencies` chain.

## Mesh-set naming (XTK/GEN convention, not arbitrary)

Set names encode phase and cut status; you don't define them, XTK generates them:
- `HMR_dummy_c_p2` / `HMR_dummy_n_p2` — bulk cells of **phase 2**; `c` = intersected/cut by interface,
  `n` = non-cut. `pN` = phase index from `phase_table`.
- `SideSet_4_n_p2` — external boundary side set 4 (HMR numbers domain sides), phase 2.
- `iside_b0_2_b1_0` — single-sided interface set between phase 2 and phase 0.
- `dbl_iside_p0_1_p1_2` — double-sided interface (phase 1 ↔ phase 2), used for Nitsche coupling IWGs.
- `ghost_p1` — ghost-stabilization set for phase 1 (only when `ghost_stab`/`tUseGhost`).

## Decoding any specific key

A key's **meaning + default** lives in `projects/PRM/src/fn_PRM_<MODULE>_Parameters.hpp` (e.g.
`fn_PRM_XTK_Parameters.hpp` shows `decomposition_type` defaults to `"conformal"`). A key whose value is a
C++ **enum** (`IWG_type`, `constitutive_type`, `Field_Type`, …) — its allowed values are the enum in
`projects/PRM/ENM/src/cl_FEM_Enums.hpp` (and `cl_XTK_Enums.hpp`, etc.), and the enum→class mapping is in
the matching factory, `projects/FEM/INT/src/<SEC>/cl_FEM_<SEC>_Factory.cpp`. Open these only when a
specific key's value or valid set is actually in question.

## Common mistakes

- Reading the file linearly instead of tracing string links — you miss the GEN↔FEM and IWG↔CM↔Property graph.
- Assuming `aCriteria` order matches the FEM IQI definition order. It matches **GEN's `IQI_types`** order.
- Treating `IQI_name`/`constitutive_name` as meaningful — they're arbitrary labels; the `*_type` enum is what computes.
- Re-deriving a key's default by guessing instead of opening its `fn_PRM_*` list.
- Reading a perimeter as a bug because `IQIPerimeter_*` uses `IQI_type::VOLUME` over an *interface* set — integrating 1 over a lower-dimensional set yields its measure (the perimeter). Intentional.
