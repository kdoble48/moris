# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

> Scope: the **core MORIS C++ codebase** (`moris/`). The parent workspace `../CLAUDE.md` covers env
> setup, `pymoris`, `torchtopo`, and `plato`. This file points at the authoritative docs and gives the
> big-picture map that requires reading many files to reconstruct.

MORIS (Multi-physics Optimization Research and Innovation System) is a CU research code for
**PDE-constrained shape/topology optimization** via an isogeometric **XFEM** formulation. Everything
lives under the C++ `moris` namespace.

## Where the docs are

Don't re-derive what's already written — read the relevant doc first:

| Topic | Source of truth |
|-------|-----------------|
| Environment activation | `../CLAUDE.md`, `../.claude/skills/moris-env/SKILL.md` |
| Building | `../.claude/skills/build-moris/SKILL.md`, `share/doc/Building_moris.dox`, `share/install/Build_Instructions.txt` |
| Running / writing tests | `../.claude/skills/run-ctest/SKILL.md`, `../.claude/skills/test-moris/SKILL.md`, `share/doc/Dev_Running_Tests.dox`, `share/doc/Dev_Writing_Tests.dox` |
| Coding style | `share/doc/Dev_Coding_Style.dox` + `.clang-format` (authoritative) |
| Adding a library / executable / project | `share/doc/Dev_Adding_Libraries.dox`, `Dev_Adding_Executables.dox`, `Dev_Adding_Projects.dox` |
| CMake details / known build-system debt | `share/doc/Dev_CMake_*.dox`, `CMAKE_AUDIT.md` |
| Architecture overview & workflows | `share/doc/Overview.dox`, `Modules.dox`, `Wrk_*_Workflow.dox`, `../.claude/rules/agent-knowledge-overview.md` |

## Day-one commands

Env must be active first (see above). Builds are out-of-source (in-source is blocked); `build_opt/`
(release) and `build_dbg/` (debug) already exist.

```bash
cd build_opt && ninja -j4                              # incremental — prefer this; full builds are slow
ctest -j2 --output-on-failure                          # all tests (Catch2)
ctest -L XTK --output-on-failure                       # by label (label == package name)
./projects/XTK/test/bin/XTK-test.exe "[ENRICH_1]" -s   # one test by Catch2 tag — fastest dev loop
```

## Naming conventions (at a glance)

Filename prefix encodes the entity: `cl_` class, `fn_` free function, `op_` operator, `st_` struct,
`enums`/`ts_` enums/typedefs; templates often have a `.tpp` companion. Names are package-prefixed —
`cl_FEM_IWG.hpp`, `cl_XTK_Cut_Integration_Mesh.hpp`, `fn_PRM_FEM_Parameters.hpp` — match the directory's
prefix. Globals `gMorisComm` (`Comm_Manager`) and `gLogger` (`Logger`) are `extern` everywhere.

## The big picture (why so many packages)

Entry point: `fn_WRK_Workflow_Main_Interface` (`projects/mains/main.cpp`). **WRK** orchestrates:
a `Workflow` (default `Workflow_HMR_XTK`, built by `cl_WRK_Workflow_Factory`) drives one run; each module
implements the `Performer` interface (`cl_WRK_Performer.hpp`) and is registered in
`cl_WRK_Performer_Manager.hpp`. Per **OPT** design iteration:

```
new ADVs → GEN.set_advs()  (update geometry, phase table, properties; PDV → nodal fields)
         → HMR refinement   (refine cells the level-set interface cuts)
         → MTK mesh          (interpolation + integration mesh)
         → XTK cut           (cut background cells on the level-set; interface/Nitsche elements)
         → MDL assemble      (per element: IWG residual/Jacobian via Field Interpolators, CMs, SPs)
         → SOL solve         (TSA time loop → NLA Newton → DLA linear solve)
         → OPT criteria       (IQIs evaluate objective/constraints)
         → compute_dcriteria_dadv  (adjoint or finite-difference sensitivities) → OPT update (GCMMA/SQP)
```

A run is configured entirely by **PRM** parameter lists (`parameters.hpp` → `fn_PRM_<MODULE>_Parameters`,
loaded via `Library_IO`). Exposing a new tunable means editing the `fn_PRM_*` list **and** the consuming
module **and** a `UT_*` test — this three-layer change is the recurring pattern (see the branch's
`sdf_reinit` commits).

## Module map

Infrastructure: **MRS** (root services: `COR` types, `CNT` containers, `IOS`/`Logger`, `ASR` asserts) ·
**LINALG** (`Matrix`/`Vector` over Armadillo or Eigen) · **COM** (MPI) · **ALG** (small helpers) ·
**TOL** (tools) · **PRM** (parameter lists).

Mesh / geometry / XFEM: **MTK** (mesh abstraction + STK backend; interpolation + integration mesh) ·
**HMR** (B-spline hierarchical background mesh, level-set h-refinement) · **GEN** (Geometry Engine;
`ADV/` optimizer design variables, `PDV/` ADV→nodal mapping, `SDF` signed-distance fields, phase table) ·
**XTK** (cuts cells against the level-set → `Cut_Integration_Mesh`, enrichment, ghost) · **MIG** (mesh editor, periodic BCs).

Physics / solve / optimize: **FEM/INT** is the physics core — `IWG/` Integrands of the Weak form of
Governing equations, `CM/` Constitutive Models, `SP/` Stabilization Parameters, `IP/` Field
Interpolators, `IQI/` Integrated Quantities of Interest. **FEM/MDL** assembles the global system;
**FEM/MSI** is the Model-Solver Interface (DOF numbering, design-variable sensitivities); **FEM/VIS**
output. **SOL** = `DLA` linear (PETSc/Belos), `NLA` nonlinear (Newton), `TSA` time-stepping, `SOL_CORE`.
**OPT** = algorithms (GCMMA/SQP/LBFGS/Sweep) + sensitivities. Runnable cases live in `projects/EXA/`.
