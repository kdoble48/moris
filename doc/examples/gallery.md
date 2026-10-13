# Results gallery

Results computed with MORIS. Every image is rendered from the output of a MORIS run —
Exodus files from the VIS module or the XTK integration mesh — and the numbers in the captions
are read from the run logs. The level-set designs use the immersed formulation throughout: a
B-spline level-set field on an HMR background mesh, cut by XTK, with ghost stabilization and
Nitsche-enforced boundary conditions.

Click an image to open it at full size. To add a result, see
[Adding a result](#adding-a-result).

## Level-set topology optimization in 2D

:::{container} gallery-grid

```{figure} gallery/mbb_final.png
:target: ../_images/mbb_final.png
:alt: MBB half-beam level-set design

**MBB half-beam.** Compliance minimization under a mass constraint, GCMMA, linear B-splines,
perimeter penalty ramped from 0.5% to 5% of the objective. Strain energy 77.91 → 76.63 over
151 evaluations; one connected body. The design was started from a MORIS SIMP solution.
```

```{figure} gallery/mbb_evolution.webp
:target: ../_images/mbb_evolution.webp
:alt: Animation of the MBB design evolution

**Design evolution** of the same run, from the SIMP-seeded start to evaluation 151 (every
second evaluation shown).
```

```{figure} gallery/mbb_frames.png
:target: ../_images/mbb_frames.png
:alt: MBB design at evaluations 1, 20, 40, 60, 100, and 151

**Snapshots** of the same run at evaluations 1, 20, 40, 60, 100, and 151.
```

```{figure} gallery/mbb_cheese_final.png
:target: ../_images/mbb_cheese_final.png
:alt: Full MBB beam started from a regular hole pattern

**Full MBB beam from a hole pattern.** A pure level-set run started from a regular 9 × 3 array
of holes, with the level set redistanced to a signed-distance field every 10 iterations.
Strain energy 210.2 → 76.14 over 167 evaluations. Holes in a level-set design can merge but
not nucleate, so the layout is reached by merging the seeded holes.
```

```{figure} gallery/cantilever_bridge.png
:target: ../_images/cantilever_bridge.png
:alt: Cantilever and bridge designs

**Cantilever and bridge,** each with and without signed-distance upkeep and multilinear
level-set interpolation. Cantilever strain energy 35.66 and 36.17; bridge 2.177 and 2.178
(151 evaluations each, started from SIMP density solutions).
```

```{figure} gallery/mbb_refinement.png
:target: ../_images/mbb_refinement.png
:alt: MBB designs on 60x20, 120x40, and 240x80 meshes

**Mesh refinement.** The same problem on 60 × 20, 120 × 40, and 240 × 80 meshes. Top: the
length scale fixed in elements, so finer meshes resolve finer members; bottom: the length scale
fixed physically, so the layout persists under refinement (strain energy 68.52, 68.55, 67.95).
```

```{figure} gallery/mbb_adaptive.png
:target: ../_images/mbb_adaptive.png
:alt: MBB design on an adaptively refined HMR mesh

**Adaptive refinement.** HMR refines the background mesh along the moving interface and
remeshes during the run. Top: one refinement level, strain energy 76.73 on 12,999 cells.
Bottom: two levels, 78.19 on 33,036 cells after eight remeshes; its interface is visibly less
smooth.
```

```{figure} gallery/lbracket_order.png
:target: ../_images/lbracket_order.png
:alt: L-bracket designs with linear and quadratic B-splines

**L-bracket, linear vs quadratic B-splines,** at 32% of the L area. Normalized compliance
c\* = 273.5 (linear) and 254.7 (quadratic). A conformal re-analysis of the final geometries
gives 273.5 and 253.3.
```

:::

## Level-set topology optimization in 3D

:::{container} gallery-grid gallery-wide

```{figure} gallery/boxbeam_3d.jpg
:target: ../_images/boxbeam_3d.jpg
:alt: 3D cantilever box-beam design, exterior and quarter cutaway

**3D cantilever box beam,** 48 × 16 × 16 background mesh, mass budget 20% of the box. Exterior
(left) and a quarter cutaway (right). Strain energy 7.294 after 155 evaluations, about three
hours on 8 MPI ranks. The level set was seeded from a SIMP density solution.
```

:::

## Verification

:::{container} gallery-grid

```{figure} gallery/kirsch_order.png
:target: ../_images/kirsch_order.png
:alt: Kirsch plate convergence by B-spline order

**Kirsch plate with a hole, immersed.** Convergence of the hoop stress and the displacement
against the analytic Kirsch solution for B-spline orders 1–3. The displacement converges at
O(h^2.2), O(h^2.8), and O(h^3.4); the stress rates are limited by the piecewise-linear
representation of the cut circle.
```

```{figure} gallery/hertz_profile.png
:target: ../_images/hertz_profile.png
:alt: Hertz contact pressure profile against the analytic solution

**Hertz contact,** immersed Nitsche contact
([`Struc_Hertz_Contact`](https://github.com/kdoble48/moris/tree/develop/projects/EXA/structure/contact/Struc_Hertz_Contact)).
On a 160 × 160 mesh the peak pressure matches the Hertz value (p₀/p₀,H = 1.00), the
contact half-width is 0.98 of the analytic one, and the profile's relative L₂ error is 0.031.
Mesh convergence and load scaling: {doc}`exa/structure/Struc_Hertz_Contact`.
```

:::

## Adding a result

1. Render the image from the run's output (Exodus or XTK mesh) — never from a mock-up — and
   downsize it to at most 1400 px wide (`.jpg` for shaded 3D renders, `.png` for flat 2D).
2. Copy it to `doc/examples/gallery/` and record its source path in
   `doc/examples/gallery/SOURCES.json`.
3. Add a `figure` to this page with a caption that states the problem, the method, and the
   numbers from the run log. Report a run that stopped at its iteration limit as "after *N*
   evaluations", not as converged.
