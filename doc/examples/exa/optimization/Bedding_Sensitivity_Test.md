# Bedding sensitivity test

A finite-difference sensitivity check for linear elasticity with an elastic bedding term. Two
concentric circles (spheres in 3D) cut a unit square (cube); their radii are the design
variables. A `SWEEP` evaluates the strain energy, the strain energy including the bedding
property, and the volume, and compares their adjoint gradients with central finite
differences. The .dox names the feature under test as "Structural Linear with Bedding": the
bedding contribution of `STRUC_LINEAR_BULK` and of the `STRAIN_ENERGY` IQI.

```{exa-summary} projects/EXA/optimization/Bedding_Sensitivity_Test
```

## Problem

| | |
|---|---|
| Domain | unit square (2D) or unit cube (3D), selected by `gDim` |
| Geometry | two concentric circles (2D) or spheres (3D) centred at (0.5, 0.5[, 0.5]); outer radius 0.3, inner radius 0.1 |
| Phases | phase table `[1, 0, 1, 1]`: three of the four geometry-sign combinations are material (phase 1), one is void (phase 0) |
| Material | linear elastic, isotropic: E = 2.0, ν = 0.0, density 1.0 |
| Bedding | bedding stiffness 2.0 × 10⁻² on the material bulk |
| Supports | side set 4 (left edge), zero displacement, Nitsche-enforced |
| Load | body load (1, 0[, 0]) on the material bulk |
| Analysis | linear, static, one Newton iteration |
| Criteria | strain energy, strain energy with bedding, volume |

## Setup in MORIS

**Mesh (HMR).** 4 elements per direction on a domain of edge 1.0, one level of initial
refinement (8 elements per direction), one linear Lagrange mesh and one linear B-spline mesh.
The displacement DOFs live on B-spline mesh 0 (`UX`, `UY`[, `UZ`] set to 0 in MSI).

**Geometry and phases (GEN).** Two `CIRCLE` (2D) or `SPHERE` (3D) level-set geometries. Each
radius is given as a (lower, initial, upper) triple, which makes it an ADV: the outer radius is
bounded to [0.27, 0.33] and the inner radius to [0.09, 0.11], giving two ADVs.
`IQI_types` is `IQIBulkStrainEnergy`, `IQIBulkStrainEnergyWithBedding`, `IQIBulkVolume`
([coupling C1](../../../user-guide/interlinks.md#c1-criteria-order-gen-iqi-types-acriteria-opt-callbacks)).

**Cut and enrichment (XTK).** Conformal decomposition, enrichment of mesh 0, ghost
stabilization on, Exodus output of the integration mesh.

**Physics (FEM).** One `STRUC_LIN_ISO` constitutive model. The bulk IWG `STRUC_LINEAR_BULK`
takes the `Bedding` and `Load` properties. The support uses
`STRUC_LINEAR_DIRICHLET_UNSYMMETRIC_NITSCHE` with a `DIRICHLET_NITSCHE` penalty of 10.0, and the
ghost term `GHOST_NORMAL_FIELD` uses a `GHOST_DISPL` parameter of 0.005. The three criteria are
two `STRAIN_ENERGY` IQIs, the second with the `PropBedding` property attached, and a `VOLUME`
IQI weighted by density. FEM-internal finite differences use `POINT_3_CENTRAL` with a
perturbation of 10⁻⁵.

**Solver (SOL).** By default a direct `AMESOS_IMPL` solve without preconditioner. Test case 21
switches to PETSc: `fgmres` with `KSPTol` 10⁻¹⁰, a `mumps` preconditioner, and
`SOL_TPL_Type = Petsc`. The nonlinear algorithm runs one iteration (relative residual drop
10⁻⁹).

**Optimization (OPT).** `SWEEP` with `finite_difference_type = "all"` (forward, backward,
central), `finite_difference_epsilons = 1.0e-5`, one evaluation per ADV, and all ADVs perturbed.
The objective is the sum of the three criteria; each criterion is also returned as a
constraint, so the constraint gradients are the per-criterion gradients. Results go to
`Bedding_Sensitivity_Test_<case>_SEN.hdf5`.

**Output (VIS).** `Bedding_Sensitivity_Test_<case>.exo` with nodal displacements and the three
criteria as global fields, written every evaluation.

## What the test checks

The test reads `constraint_gradients eval_1-1 analytical` and
`constraint_gradients eval_1-1 epsilon_1-1 fd_central` from the sweep HDF5 file, requires equal
sizes, and compares every entry. Where |FD| > 10⁻⁶ the relative difference must be below 10⁻⁴;
otherwise the analytical value must match FD within an absolute margin of 10⁻⁶. The Exodus
checks in `check_results` are commented out.

| Proc count | Case index | Dimension | Solver |
|---|---|---|---|
| 1 | 0 | 2D | Amesos |
| 1 | 1 | 3D | Amesos |
| 4 | 10 | 2D | Amesos |
| 4 | 11 | 3D | Amesos |
| 4 | 21 | 3D (`gDim` is still 3) | PETSc; only when `MORIS_HAVE_SLEPC` is defined |

On 4 processors the comparison runs on rank 0 only.

## Run it

```bash
ctest -R Bedding_Sensitivity_Test --output-on-failure     # in a build with MORIS_USE_EXAMPLES=ON
```

The ctest names are `Bedding_Sensitivity_Test-1-procs` and `Bedding_Sensitivity_Test-4-procs`.
Before running, the test removes `*exo*` and `*hdf5*` files from the working directory.

## Notes

- {doc}`Contact_Sensitivity_Test <Contact_Sensitivity_Test>` follows the same sweep-and-compare
  pattern for a contact interface.
- The procedure for writing a sensitivity sweep is described in
  {doc}`../../../user-guide/workflow`.
