# Example gallery

Every example in `projects/EXA` that leaves output behind, rendered from that output: the solution
field of the VIS output where there is one, otherwise the XTK integration mesh colored by phase.
Each card opens the example's page, where the scene can be rotated and zoomed in the browser.

## Example map

Every example, linked to the ones it most resembles. Two examples are connected when their decks
share detected features (the physics, the optimizer, the geometry representation, ghost
stabilization, Nitsche enforcement, the heat method, adaptive refinement, …) and deck anatomy
(DOFs, solver stack, mesh size, block structure); an edge is wider the more they share. Drag to
pan, scroll to zoom, hover a node or an edge for what it stands for, click a node to open the
example's page. The menu highlights every example with one feature; the legend entries hide and
show whole families. Each example page carries the same graph reduced to its own neighbours.

```{exa-graph}
```

## Scenes

```{include} _generated/example_gallery.md
```
