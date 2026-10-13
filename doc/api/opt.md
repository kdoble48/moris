# OPT — optimization

`opt::Manager` runs the optimization loop: it owns the `opt::Problem` (objective, constraints,
and their gradients, evaluated through the criteria interface) and the list of algorithms.

```{doxygenclass} moris::opt::Manager
```

```{doxygenclass} moris::opt::Problem
```

## Algorithms

```{doxygenclass} moris::opt::Algorithm
```

```{doxygenclass} moris::opt::Algorithm_MMA
```

```{doxygenclass} moris::opt::Algorithm_ROL
```
