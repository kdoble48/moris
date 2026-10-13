# Module reference

A curated index of the types and keys each module exposes. Authoritative defaults live in
`projects/PRM/src/fn_PRM_<MODULE>_Parameters.hpp` and valid enum values in
`projects/PRM/ENM/src/`; trust those files when they disagree with this page.

## HMR — background mesh

The hierarchical B-spline background mesh. Key parameters: `number_of_elements_per_dimension`,
`domain_dimensions`, `domain_offset`, `pattern_initial_refinement`, `refinement_buffer`,
`staircase_buffer`, and the sub-lists `HMR::LAGRANGE_MESHES` and `HMR::BSPLINE_MESHES` (one
entry per mesh, each with its `order(s)` and `pattern`).

A deck commonly defines two B-spline meshes: one discretizing the DOFs and one discretizing the
level-set design field. Two practical rules from the legacy documentation:

- The geometry mesh should not have more nodes than the DOF mesh, or the design problem becomes
  ill-conditioned.
- A quadratic B-spline geometry mesh gives a smoother interface than a linear one.

## XTK — immersed cut

| Group | Keys (default) |
|---|---|
| Decomposition | `decompose` (true), `decomposition_type` ("conformal"), `octree_refinement_level`, `triangulate_all` (false), `ig_element_order` (1) |
| Enrichment | `enrich` (true), `use_SPG_based_enrichment` (false), `basis_rank` ("bspline"), `enrich_mesh_indices` ("0") |
| Ghost | `ghost_stab` (false), `visualize_ghost` (false) |
| Agglomeration | `activate_basis_agglomeration`, `activate_cell_agglomeration` (false), `volume_fraction` |
| Output | `output_path`, `output_file` ("xtk_temp.exo"), `exodus_output_XTK_ig_mesh`, `write_enrichment_fields` |

XTK **produces** the mesh-set names that FEM consumes; see coupling C5 in {doc}`interlinks`.

## GEN — geometry and design

- **Geometry types:** `LEVEL_SET`, `SURFACE_MESH`, `VOXEL`.
- **Field types:** `CONSTANT`, `LINE`, `CIRCLE`, `SUPERELLIPSE`, `PLANE`, `SPHERE`,
  `SUPERELLIPSOID`, `SCALED_FIELD`, `COMBINED_FIELDS`, `NODAL_FROM_FILE`,
  `SIGNED_DISTANCE_OBJECT`, `SIGNED_DISTANCE_IMAGE`, `USER_DEFINED`.
- **Field arrays** tile a primitive over a box (`lower/upper_bound_*`, `number_of_fields_*`)
  — the usual way to seed holes in a level-set design.
- **Design variables (ADVs)** enter in one of two ways: a level-set field with
  `discretization_mesh_index` set (its B-spline coefficients are the ADVs, bounded by
  `discretization_lower/upper_bound`), or a GEN property that publishes a PDV type (`LS1`,
  `DENSITY`, `ELASTIC_MODULUS`, …) consumed by FEM properties.
- **`phase_table`** maps each combination of geometry signs to a bulk phase index.
- **`IQI_types`** lists the criteria evaluated for OPT, in index order.

:::{admonition} Design-free geometry in an optimization deck
:class: tip

A geometry that should not carry ADVs (a fixed frame, a static hole) still needs a stored field
in an optimization deck: set `discretization_mesh_index = -1`. The analytic path (`-2`) yields
undefined design gradients for that geometry.
:::

## FEM — physics

Sub-lists are filled bottom-up, and each layer references the earlier ones by name:
`PROPERTIES → CONSTITUTIVE_MODELS → STABILIZATION → IWG → IQI → COMPUTATION`.

**Weak-form integrands (`IWG_Type`, ~63)**, by family:

- *Diffusion:* `SPATIALDIFF_BULK`, Dirichlet/Robin/interface Nitsche variants (symmetric and
  unsymmetric), `_NEUMANN`, `_CONVECTION`, `_RADIATION`, `_GGLS_PC`, `_VW_GHOST`.
- *Linear elasticity:* `STRUC_LINEAR_BULK`, Dirichlet/interface Nitsche variants, `_NEUMANN`,
  `_VW_GHOST`, mixed pressure formulations, Lagrange-multiplier interface variants.
- *Nonlinear elasticity:* `STRUC_NON_LINEAR_BULK_{SE,PF,CAUCHYEPS}`, geometric stiffness,
  Nitsche variants per stress measure.
- *Incompressible and compressible Navier–Stokes, Spalart–Allmaras turbulence.*
- *Contact:* linear Nitsche and penalty variants; nonlinear unbiased variants.
- *Misc:* `TIME_CONTINUITY_DOF`, `GHOST_NORMAL_FIELD`, `L2`, `HJTEST`, `USER_DEFINED`.

**Quantities of interest (`IQI_Type`, ~51):** `VOLUME`, `STRAIN_ENERGY`, `VOLUME_FRACTION`,
`DOF`, `MAX_DOF`, `PROPERTY`, `L2_ERROR_ANALYTIC`, `H1_ERROR_ANALYTIC`, `J_INTEGRAL`,
stress measures (`VON_MISES_STRESS`, `MAX_VON_MISES_STRESS`, `ZIENKIEWICZ_ZHU_VON_MISES_STRESS`,
…), traction and contact quantities, flow quantities (`LIFT_COEFF`, `DRAG_COEFF`,
`MASS_FLOW`, …), `HOMOGENIZED_CONSTITUTIVE`, `MIN_FEATURE_SIZE`, `HEAT_METHOD_PENALTY`.

**Constitutive models (~21):** `DIFF_LIN_ISO`; `STRUC_LIN_ISO`, `STRUC_LIN_MT`,
`STRUC_LIN_ISO_PRESSURE`, `STRUC_LIN_ISO_DAMAGE`; `STRUC_NON_LIN_ISO_SAINT_VENANT_KIRCHHOFF`,
`STRUC_NON_LIN_ISO_COMPRESSIBLE_NEO_HOOKEAN_{BONET,WRIGGERS}`; incompressible and compressible
fluids; turbulence variants.

**Stabilization parameters (~28):** `DIRICHLET_NITSCHE`, `ROBIN_NITSCHE`, `NITSCHE_INTERFACE`,
`GHOST_DISPL`, `GHOST_VW`, `GHOST_NORMAL_FIELD`, `{LEADER,FOLLOWER}_WEIGHT_INTERFACE`,
SUPG/GGLS/crosswind variants, flow ghost parameters, contact penalties.

## SOL — solvers

The SOL list is a chain of indexed sub-lists: time solvers reference nonlinear solvers, which
reference linear solvers, which reference linear algorithms and preconditioners.

| Layer | Types |
|---|---|
| Linear (DLA) | `AMESOS_IMPL` (direct), `BELOS_IMPL` (Krylov), `AZTEC_IMPL`, `PETSC`, `EIGEN_SOLVER`, `SLEPC_SOLVER`, `ML` |
| Preconditioner | `NONE`, `IFPACK`, `ML`, `PETSC` |
| Nonlinear (NLA) | `NEWTON_SOLVER`, `NLBGS_SOLVER` (staged multiphysics), `ARC_LENGTH_SOLVER` |
| Time (TSA) | `MONOLITHIC`, `STAGGERED` |
| Sensitivities | `ADJOINT`, `DIRECT` |

:::{admonition} Set the Amesos backend explicitly
:class: warning

The default Amesos `Solver_Type` is `Amesos_Pardiso`, which many Trilinos installations do not
include; a deck relying on it aborts at the first solve. Set `Solver_Type` to `Amesos_Klu`
(serial) or `Amesos_Superludist` (parallel) explicitly.
:::

Newton robustness options on the nonlinear list include relaxation
(`SolverRelaxationType`), load control (`SolverLoadControlType`), and pseudo-time continuation
(`SolverPseudoTimeControlType`).

## OPT — optimization

Problem keys: `is_optimization_problem`, `workflow` ("HMR_XTK"), `problem` ("user_defined"),
`restart_file`, `finite_difference_type`, `library`.

| Algorithm | Notes |
|---|---|
| `GCMMA` | Svanberg's globally convergent MMA (external library). Keys: `max_its`, `step_size`, `penalty`, `norm_drop`, `asymp_adapt0/b/c` |
| `MMA` | Native in-tree MMA with the same parameter shape as GCMMA |
| `ROL` | Trilinos ROL: trust region / Lin–Moré for bound constraints, augmented Lagrangian with constraints; `xml_file` passes a full ROL parameter list |
| `LBFGS` | Limited-memory BFGS with bounds |
| `SQP` | SNOPT (requires `MORIS_HAVE_SNOPT`) |
| `SWEEP` | **Not a design run**: evaluates criteria and gradients over ADV perturbations and compares analytic and finite-difference sensitivities |

With ROL, `max_its` bounds only the augmented-Lagrangian outer loop; cap the inner
trust-region iterations as well, or one outer iteration can consume hundreds of forward solves.

## MSI and VIS

**MSI** enables DOF types (`UX`, `UY`, `UZ`, `P`, `TEMP`, `VX`, …) and assigns each to a B-spline
mesh index: `tMsi.set( "UX", 1 )`.

**VIS** writes Exodus output: `File_Name`, `Mesh_Type` (`STANDARD`, `STANDARD_WITH_OVERLAP`,
`FULL_DISCONTINUOUS`, `FULL_DISCONTINUOUS_WITH_OVERLAP`), `Set_Names`, `Field_Names` with matching `Field_Type` (`NODAL`, `ELEMENTAL_AVG`,
`GLOBAL`, …) and `IQI_Names`, and `Save_Frequency`.
