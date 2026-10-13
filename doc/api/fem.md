# FEM — physics

The physics layer. An **IWG** (integrand of the weak form of the governing equations)
contributes a residual and Jacobian on a set; an **IQI** (integrated quantity of interest)
evaluates a criterion; both draw on **constitutive models**, **stabilization parameters**, and
**properties**. `mdl::Model` assembles them over the cut mesh.

## Weak-form integrands

```{doxygenclass} moris::fem::IWG
```

## Quantities of interest

```{doxygenclass} moris::fem::IQI
```

## Constitutive models

```{doxygenclass} moris::fem::Constitutive_Model
```

## Stabilization parameters

```{doxygenclass} moris::fem::Stabilization_Parameter
```

## Properties

```{doxygenclass} moris::fem::Property
```

## Model assembly

```{doxygenclass} moris::mdl::Model
```
