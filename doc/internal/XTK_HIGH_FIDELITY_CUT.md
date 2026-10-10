# XTK high-fidelity cut: cell-order intersections (`use_cell_order_intersection`)

2026-10-10, branch `feature/xtk-hifi-cut` (base `fix/xiga-deck-bugs`). Context and measurements of the problem:
`runs/studies/lunch_and_learn_xiga/docs/INTERFACE_JAGGEDNESS.md` (§1 mechanism, §7 sized options). File:line references
are to this branch.

## 1. Problem

At B-spline order p ≥ 2 the level set φ_h is C^(p-1), but with `use_multilinear_interpolation` XTK cuts the bilinear
interpolant of the 4 corner values of each Lagrange cell:

- The multilinear Newton intersection builds a QUAD4/HEX8 interpolator from the first 4/8 background nodes
  (`GEN/src/ADV/cl_GEN_Level_Set_Geometry.cpp:174-235`).
- `dxi/dphi` returns 0 for every node past the corners (`GEN/src/PDV/cl_GEN_Intersection_Node_Bilinear.cpp:119-210`).
- The geometry sets a global linear override on every derived node (`cl_GEN_Design_Field.cpp:61-64`, consumed at
  `cl_GEN_Derived_Node.cpp:35-39`). The template centre node therefore gets the corner average.

On the X9q_quadratic_rcp final design (order 2, 150 iterations), the facet end points are 0.19 h (p99) and 0.92 h (max)
from the φ_h zero set. A closed φ_h island is dropped entirely, at 5.2 h.

## 2. Options evaluated

### (C) "More geometric data at the centre node" (user proposal)

**Q1. What does the centre node get today?** It is created by `Geometry_Engine::create_new_derived_nodes` with all
parent-cell nodes and the cell order (`cl_GEN_Geometry_Engine.cpp:444-500`).

| Case | Centre φ | Spoke intersection |
|---|---|---|
| Multilinear on (any order) | Mean of the 4 corners (global linear override) | Newton on the corner bilinear; the centre value only fixes the sign |
| Multilinear off, order 1 | QUAD4 at (0, 0) = corner mean = exact φ_h | Linear between the end values |
| Multilinear off, order 2/3 | QUAD9/16 basis = exact φ_h (node 8 for QUAD9) | Linear between the end values |

So (C), with multilinear off, is today's default path at order ≥ 2. With multilinear on, giving the centre the exact
value alone changes nothing except the sign test. The Newton still runs on the corner bilinear, so (C) needs the
intersection change of (B) to take effect.

**Q2. Cracks?** None:
- The spokes are interior to the cell.
- An outer-edge intersection depends only on φ along that edge. Its restriction is the same 1D polynomial from both
  neighbouring cells, because QUAD9 restricted to an edge uses only that edge's 3 nodes and φ_h is continuous.
- XTK creates one intersection per parent edge anyway.

**Q3. Sensitivities (C).** The centre is a derived node with Basis_Nodes (N_k(0), k = 0..8). The linear spoke
intersection chains `dxi/dphi_parent · dphi_centre/dphi_k = N_k(0)` through `Design_Field::get_dfield_dadvs`
(`cl_GEN_Design_Field.cpp:155-190`). No new code is needed.

### (B) Intersect the cell's own-order interpolant

Run Newton on φ(ξ) = Σ_k N_k(ξ) φ_k over all cell nodes, with the QUAD9/16 (HEX27/64) Lagrange basis. For a B-spline
field whose Lagrange order equals the B-spline order on the same mesh, this interpolant reproduces φ_h exactly. The
intersections then lie on the true zero set. (C) is included, because the centre gets the exact value.

Sensitivities: the root ξ* of r(ξ) = Σ N_k(ξ(s)) φ_k − c along the edge (edge coordinate s) gives

    ds*/dφ_k = − N_k(ξ*) / (∂r/∂s),   ∂r/∂s = ½ (ξ_2 − ξ_1) · ∇_ξ N · φ

`Intersection_Node_Level_Set::append_dcoordinate_dadv` (`cl_GEN_Intersection_Node_Level_Set.cpp:43-98`) already loops
over all field basis nodes and multiplies by dφ_k/dADV. Generalising the bilinear formula to all k is the whole change.
Parent coordinates do not depend on the ADVs (corners and centre), as before.

### (A) N×N sub-cell subdivision

Subdivide each cut cell into N×N sub-cells before the 4-TRIS split. This adds sample points, so it reduces the chord
turns that (B) leaves.

### Offline prediction

Script: `studies/lunch_and_learn_xiga/scripts/cut_variant_prediction.py`. It rebuilds the XTK 2D cut per variant from a
run's final φ_h. Its `mlin` row reproduces the MORIS measurement (0.195 / 0.918 h, island at 5.2 h), which validates
the model. Data: `runs/studies/lunch_and_learn_xiga/analysis/jaggedness/cut_variant_prediction_X9q_rcp.json`.
Distances are in units of h.

| Variant (X9q_rcp final, order 2) | IG tris | end p99 / max | facet p99 / max | missed φ_h p99 / max | turn at h-edges p99 | turns > 10° | facet-vs-tangent p99 |
|---|---|---|---|---|---|---|---|
| mlin (today) | 2284 | 0.195 / 0.918 | 0.226 / 0.940 | 0.378 / 5.23 | 53.3° | 5.4% | 26.3° |
| C (centre exact, linear edges) | 2284 | 0.146 / 0.339 | 0.153 / 0.863 | 0.145 / 0.20 | 19.2° | 9.9% | 20.6° |
| **B (cell-order Newton)** | 2284 | **0 / 0** | 0.028 / 0.112 | 0.037 / 0.09 | 24.2° | 9.2% | 1.0° |
| A2 (2×2, linear) | 9136 | 0.031 / 0.310 | 0.034 / 0.310 | 0.036 / 0.14 | 9.7° | 4.5% | 14.5° |
| A2 + B | 9136 | 0 / 0 | 0.009 / 0.027 | 0.012 / 0.03 | 14.7° | 4.0% | 0.3° |
| A4 + B | 36544 | 0 / 0 | 0.002 / 0.022 | 0.004 / 0.02 | 5.7° | 1.1% | 0.1° |

The columns, in units of h:
- "facet": the distance of 9 points per facet from φ_h = 0.
- "missed": the distance from the true contour to the cut, which catches dropped islands.

**Choice: (B).** It gives the largest precision gain per effort, at no change in cell count:
- The worst-case facet error falls by 8.4x (0.94 to 0.11 h).
- The dropped island is recovered (5.2 to 0.09 h).
- The facet direction error falls by 26x.
- The cost is about 300 lines in GEN, with no XTK change.

It is robust for optimization:
- The sensitivities stay analytic.
- The Newton iteration is safeguarded by a bracket and cannot abort.
- The default is unchanged.

(C) alone captures part of the gain, but it is already the multilinear-off path and keeps 0.86 h facet errors. (A) adds
4x the integration cells and an XTK template. It is the next step for smoothness, because (B) does not reduce turn
counts (§5).

The "move the centre node" extension, the earlier N2 corner idea, needs a geometry-dependent template node:
- a new node type with coordinate sensitivities;
- XTK's `has_geometric_independent_vertices` = true for regular subdivision (`cl_XTK_Regular_Subdivision_Interface.cpp:38-42`);
- FEM shape sensitivities for a non-interface node.

A2 + B reaches 14.7° h-edge turns without any of this.

## 3. Implementation (B)

| Change | Where |
|---|---|
| PRM key `use_cell_order_intersection` (default false) on design fields | `PRM/src/fn_PRM_GEN_Parameters.hpp:208` |
| `Field_Parameters::mUseCellOrderIntersection`. Read only if present, so decks compiled against older PRM headers run. Implies multilinear (same code paths for derived-node values and ADV ids) | `GEN/src/ADV/cl_GEN_Design_Field.cpp:22-37`, `.hpp:30` |
| Cell-order mode does **not** set the global linear override, and errors if a corner-multilinear geometry already set it | `cl_GEN_Design_Field.cpp:53-64` |
| `Derived_Node::clear_override_linear_interpolation()` / `get_override_linear_interpolation()` | `GEN/src/PDV/cl_GEN_Derived_Node.cpp:168-180` |
| Dispatch | `cl_GEN_Level_Set_Geometry.cpp:170-173` |
| `compute_intersection_local_coordinate_cell_order`: cell interpolant over all nodes; unit-norm scaling with an all-zero guard (MB-020); Newton safeguarded by bisection on the sign-change bracket (100 iterations at most, residual 1e-13, bracket 1e-14; never `MORIS_ERROR`s) | `cl_GEN_Level_Set_Geometry.cpp:463-594` |
| `create_cell_interpolation(dim, nodes)`: QUAD4/9/16, HEX8/27/64 | `cl_GEN_Intersection_Node_Bilinear.cpp:53-77` |
| `get_dxi_dfield_from_ancestor`: −N_k(ξ*) / (∂r/∂s) for every cell node | `cl_GEN_Intersection_Node_Bilinear.cpp:85-117` |

XTK passes all parent-cell nodes to `queue_intersection` (`XTK/src/cl_XTK_Integration_Mesh_Generator.hpp:215-218`), so no
XTK change is needed. The setting is per geometry. Another geometry in the same run that uses multilinear without cell
order is rejected. Geometries with multilinear off are unaffected, except that their derived nodes also keep the cell
order, because the override is global; the MBB pads are one example.

## 4. Verification

**Unit tests** (`GEN/test/ut_GEN_Intersection_Cell_Order.cpp`, tag `[cell_order_intersection]`):
- **Precision.** An analytic circle (r = 0.63) on QUAD9 meshes with h = 0.5 and 0.25. All 4 outer edges and 4 spokes of
  every cell are queued, with a real centre derived node. Maximum distance of the intersections from the circle:

  | Interpolation | h = 0.5 | h = 0.25 | Observed rate |
  |---|---|---|---|
  | multilinear | 5.39e-2 | 1.28e-2 | 2.08 |
  | linear (C) | 4.71e-2 | 1.12e-2 | — |
  | cell order | 7.07e-3 | 6.79e-4 | 3.38 |

- **Order 1.** Cell order equals multilinear to 1e-10 (QUAD4).
- **Sensitivities.** Quadratic B-spline circle on 3×3 QUAD9 cells, 28 intersection nodes on outer edges and spokes.
  Analytic d(x, y)/dADV, summed per ADV, against central FD (ε = 1e-6): max error / max |analytic| = 4.0e-10.

**Default unchanged:**
- GEN-test: 81109 assertions in 26 cases before, 81759 in 28 after (the 650 new assertions in the 2 new cases), all passing.
- XTK-test: 3285 in 29 before and after, all passing.
- End to end: the worktree binary without the key reproduces X9q_quadratic_rcp, built with the main LTO binary,
  bit for bit over the first 11 objective values (15 digits printed).

**End to end** (`runs/studies/lunch_and_learn_xiga/X80_hifi/`):
- Setup: order 2, X9r SIMP seed, native GCMMA + in-place SDF every step + multilinear, 10 iterations, nice 10,
  OMP_NUM_THREADS = 2.
- Peak memory: 0.4 GB (ref) and 0.6 GB (hifi). Wall time: under 1 min per arm.

| Arm | end p99 from φ_h | facet-vs-φ_h tangent rms | turns > 10° (edge / interior vertices) | Objective at 11 evals |
|---|---|---|---|---|
| X9q_quadratic_rcp (150 its, reference) | 0.192 h | 6.0° | 5.4% (11.0 / 1.2) | — |
| X80 ref (10 its, key off) | 0.142 h | 5.5° | 5.4% (11.7 / 0.9) | 1.00046 |
| X80 hifi (10 its, key on) | **0.000 h** | **0.3°** | 15.4% (12.4 / 17.6) | 0.99845 |

The predictor's B variant on the hifi field gives 15.6% turns above 10° and an end distance of 0, matching MORIS.

SWEEP: central FD on 10 interface ADVs at the hifi iteration-10 design, with the perimeter weight pinned.

| Arm | Objective, ε = 1e-4 | Objective, ε = 1e-5 | Constraint, ε = 1e-4 | Constraint, ε = 1e-5 |
|---|---|---|---|---|
| hifi | 3.8e-5 | 9.8e-5 | 8.2e-6 | 8.2e-8 |
| ref | 9.1e-6 | 4.8e-5 | 7.8e-7 | 1.2e-8 |

Each entry is the max per-ADV |analytic − FD| / |analytic|. Summary: `X80_hifi/sweep_summary.txt`.

## 5. What (B) does and does not fix

(B) makes the cut faithful to φ_h:
- The vertices lie on the zero set.
- Facet midpoints lie within 0.03 h (p99) of it.
- Islands are no longer dropped.

It does **not** make the polygon smoother. The vertices now sample the true C1 curve, whose sub-cell curvature
(design corners of radius about h, plus sub-cell wiggles of the SDF-maintained field) the corner-bilinear cut used to
flatten. That curve is drawn with 2–3 chords per cell. On the hifi field, the default cut would also turn more than 10°
at 8.4% of vertices, against 5.4% on the ref field. So the optimizer, which now sees the sub-cell features, keeps more
of them in 10 iterations. Watch this in longer runs: the perimeter regularization may need to act at that scale.

Smoothness at constant fidelity comes from more samples per cell: A2 + B gives 5.7% turns above 10° and h-edge turn
p99 13.5° on the same field.

## 6. Remaining work

1. **(A) 2×2 sub-cell template in 2D**, to use with (B). This is about 3–5 days.
   - Add a `Regular_Subdivision_16_TRIS` template with 4 mid-edge nodes and 5 interior nodes (centre and 4 sub-centres).
     The XTK template API already supports edge, face and cell nodes (`cl_XTK_Regular_Subdivision_Interface.hpp:60-98`).
   - `make_new_vertex_requests_trivial` handles only face and single cell-interior nodes. Add edge (2D facet) requests,
     and multiple interior nodes with secondary identifiers (the octree path's `register_new_request(..., hash, ...)`
     pattern, `cl_XTK_Regular_Subdivision_Interface.cpp:405-440`).
   - Register it behind an XTK key `cut_subdivision_level` (0 = today's 4-TRIS), selected in
     `Regular_Subdivision_Interface` ctor/`perform_impl_vertex_requests` (`:24-36`, `:63-74`).
   - Mid-edge nodes on edges shared with uncut cells leave hanging nodes in the IG mesh. This is non-conformal, as
     NC_REGULAR_SUBDIVISION_HEX8 already is in 3D; check side-set and ghost facet construction.
2. **3D:** (B) works as is for HEX27/64, but is not yet tested in 3D.
3. **Thesis-grade check:** a longer (≥ 40-iteration) A/B against X70_o2_lag2 (Lagrange-only refinement), on objective,
   perimeter and kink statistics.
4. `get_dfield_dcoordinates` for intersection-on-intersection stays unimplemented, as for multilinear
   (`cl_GEN_Level_Set_Geometry.cpp:600-606`).
