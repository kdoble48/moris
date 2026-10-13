# MORIS

**Multi-physics Optimization Research and Innovation System** — a C++ research platform for
PDE-constrained shape and topology optimization, built on an isogeometric formulation of the
extended finite element method (XIGA): level-set geometry on hierarchical B-spline meshes, cut by
an immersed interface, with weakly enforced boundary and interface conditions.

MORIS is developed at the University of Colorado's Aerospace Mechanics Research Center (AMReC) in
Kurt Maute's group. This site documents the [`kdoble48/moris`](https://github.com/kdoble48/moris)
fork and tracks its `develop` branch.

## Where to start

| If you want to… | Read |
|---|---|
| see what MORIS computes | {doc}`examples/gallery` · {doc}`examples/catalog` · {doc}`examples/graph` |
| build MORIS and run an example | {doc}`getting-started/index` |
| write or debug an input deck | {doc}`user-guide/index` |
| understand why a deck fails silently | {doc}`user-guide/interlinks` |
| change the C++ core, run tests, report a bug | {doc}`developer/index` |
| look up a class | {doc}`api/index` · <a href="doxygen/index.html">full Doxygen index</a> |

## How a MORIS run works

Every design iteration of an optimization passes through the same module chain:

```text
new ADVs → GEN   geometry, phase table, design properties
         → HMR   refine the B-spline background mesh where the interface passes
         → XTK   cut the mesh, enrich the basis, build interface / Nitsche / ghost sets
         → MDL   assemble residuals and Jacobians (FEM: IWGs, CMs, SPs)
         → SOL   time → nonlinear → linear solve
         → OPT   criteria (IQIs) → objective and constraints → sensitivities → update
```

The input deck configures each module with a parameter list. The modules are coupled by **name**:
a string defined in one module is referenced in another, without a compiler check. The
{doc}`user-guide/interlinks` page lists every coupling.

```{toctree}
:hidden:
:caption: Getting started

getting-started/index
getting-started/install
getting-started/build
getting-started/running
```

```{toctree}
:hidden:
:caption: User guide

user-guide/index
user-guide/deck-anatomy
user-guide/modules
user-guide/interlinks
user-guide/deck-api-v2
user-guide/workflow
```

```{toctree}
:hidden:
:caption: Gallery

examples/gallery
examples/example-gallery
examples/catalog
examples/graph
```

```{toctree}
:hidden:
:caption: Developer guide

developer/index
developer/testing
developer/coding-style
developer/fork-workflow
developer/documentation
```

```{toctree}
:hidden:
:caption: API reference

api/index
api/input
api/gen
api/hmr
api/xtk
api/fem
api/sol
api/opt
```

```{toctree}
:hidden:
:caption: Design records

design/index
```
