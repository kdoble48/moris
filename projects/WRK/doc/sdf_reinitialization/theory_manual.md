# Theory Manual: Eikonal Equation and SDF Reinitialization

## 1. Introduction

### 1.1 The Problem

In level-set topology optimization, the material domain $\Omega$ is implicitly defined by a level-set function $\phi(\mathbf{x})$:

$$\Omega = \{\mathbf{x} : \phi(\mathbf{x}) < 0\}$$

The **interface** $\Gamma$ (material boundary) is the zero level-set:

$$\Gamma = \{\mathbf{x} : \phi(\mathbf{x}) = 0\}$$

During optimization, $\phi$ evolves according to shape sensitivities. Over many iterations, $\phi$ loses its **signed distance property**, causing:

- Numerical instabilities (steep/flat gradients)
- Inaccurate curvature computation
- Spurious zero-crossings (floating artifacts)

---

## 2. Signed Distance Functions

### 2.1 Definition

A **signed distance function** (SDF) $\phi$ satisfies:

$$\phi(\mathbf{x}) = \begin{cases}
-d(\mathbf{x}, \Gamma) & \text{if } \mathbf{x} \in \Omega \\
+d(\mathbf{x}, \Gamma) & \text{if } \mathbf{x} \notin \Omega \\
0 & \text{if } \mathbf{x} \in \Gamma
\end{cases}$$

where $d(\mathbf{x}, \Gamma) = \min_{\mathbf{y} \in \Gamma} \|\mathbf{x} - \mathbf{y}\|$ is the Euclidean distance.

### 2.2 The Eikonal Property

A signed distance function satisfies the **Eikonal equation**:

$$|\nabla \phi| = 1 \quad \text{almost everywhere}$$

This means:
- The gradient magnitude is unity everywhere
- $\phi$ increases/decreases at constant rate (1 unit per unit distance)
- Isolines are equally spaced

### 2.3 Geometric Properties of SDFs

| Property           | Formula                                         | Benefit                        |
| ------------------ | ----------------------------------------------- | ------------------------------ |
| Unit gradient      | $                                               | \nabla\phi                     | = 1$          | Numerical stability |
| Normal vector      | $\mathbf{n} = \nabla\phi /                      | \nabla\phi                     | = \nabla\phi$ | Direct computation  |
| Curvature          | $\kappa = \nabla \cdot \mathbf{n} = \Delta\phi$ | Simplified formula             |
| Interface velocity | $v_n = -\partial\phi/\partial t$                | Hamilton-Jacobi interpretation |

---

## 3. The Eikonal Equation

### 3.1 Standard Form

The Eikonal equation in its general form:

$$|\nabla T(\mathbf{x})| = F(\mathbf{x})$$

where:
- $T(\mathbf{x})$ is the travel time (or distance) field
- $F(\mathbf{x})$ is the slowness (reciprocal of speed)

For SDFs with uniform speed: $F = 1$.

### 3.2 Boundary Value Problem

To compute $\phi$ as an SDF given the interface $\Gamma$:

$$\begin{aligned}
|\nabla \phi| &= 1 \quad \text{in } \mathbb{R}^n \setminus \Gamma \\
\phi &= 0 \quad \text{on } \Gamma
\end{aligned}$$

This is a **static Hamilton-Jacobi equation** (no time derivative).

### 3.3 Characteristic Curves

The Eikonal equation has characteristics (rays) emanating from $\Gamma$:

$$\frac{d\mathbf{x}}{ds} = \nabla\phi, \quad \mathbf{x}(0) \in \Gamma$$

Along each ray, $\phi$ increases linearly with arc length $s$.

---

## 4. Sussman-Sethian Reinitialization

### 4.1 The Reinitialization Equation

Sussman, Smereka, and Osher (1994) introduced a **pseudo-time evolution** to reinitialize $\phi$ to an SDF:

$$\boxed{\frac{\partial \phi}{\partial \tau} + S(\phi_0)(|\nabla\phi| - 1) = 0}$$

where:
- $\tau$ is pseudo-time (not physical time)
- $\phi_0 = \phi(\mathbf{x}, \tau=0)$ is the **initial** (distorted) level-set
- $S(\phi_0)$ is the **sign function** of the original field

### 4.2 Sign Function

The sign function determines the direction of propagation:

$$S(\phi_0) = \begin{cases}
+1 & \text{if } \phi_0 > 0 \quad \text{(outside)} \\
-1 & \text{if } \phi_0 < 0 \quad \text{(inside)} \\
0 & \text{if } \phi_0 = 0 \quad \text{(interface)}
\end{cases}$$

**Important**: $S$ depends on $\phi_0$, not $\phi(\tau)$. This anchors the interface.

### 4.3 Smoothed Sign Function

For numerical stability, use a regularized version:

$$S_\epsilon(\phi_0) = \frac{\phi_0}{\sqrt{\phi_0^2 + \epsilon^2}}$$

where $\epsilon \sim h$ (mesh size). Properties:
- Smooth transition through zero
- $S_\epsilon \to \text{sgn}(\phi_0)$ as $\epsilon \to 0$
- Prevents division by zero at interface

Alternative (Peng et al., 1999):

$$S_\epsilon(\phi_0) = \frac{\phi_0}{\sqrt{\phi_0^2 + |\nabla\phi_0|^2 \epsilon^2}}$$

### 4.4 Steady-State Solution

At steady state ($\partial\phi/\partial\tau = 0$):

$$S(\phi_0)(|\nabla\phi| - 1) = 0$$

This implies:
- $|\nabla\phi| = 1$ where $S \neq 0$ (away from interface)
- Interface preserved: $\phi = 0$ where $\phi_0 = 0$

### 4.5 Physical Interpretation

The reinitialization equation is a **Hamilton-Jacobi equation** with Hamiltonian:

$$H(\nabla\phi) = S(\phi_0)(|\nabla\phi| - 1)$$

It propagates information **outward** from the interface ($S > 0$) and **inward** toward the interface ($S < 0$), resetting $\phi$ to exact distances.

---

## 5. Weak Form for FEM Implementation

### 5.1 Strong Form

$$\frac{\partial \phi}{\partial \tau} + S(\phi_0)(|\nabla\phi| - 1) = 0 \quad \text{in } \Omega$$

### 5.2 Weighted Residual Statement

Multiply by test function $w$ and integrate:

$$\int_\Omega w \frac{\partial \phi}{\partial \tau} \, d\Omega + \int_\Omega w \cdot S(\phi_0)(|\nabla\phi| - 1) \, d\Omega = 0$$

### 5.3 Semi-Discrete Form

Using finite element approximation $\phi^h = \sum_j N_j \hat{\phi}_j$:

$$\int_\Omega N_i \frac{\partial \phi^h}{\partial \tau} \, d\Omega + \int_\Omega N_i \cdot S(\phi_0)(|\nabla\phi^h| - 1) \, d\Omega = 0$$

In matrix form:

$$\mathbf{M} \frac{d\hat{\boldsymbol{\phi}}}{d\tau} = -\mathbf{R}(\hat{\boldsymbol{\phi}})$$

where:
- $\mathbf{M}_{ij} = \int_\Omega N_i N_j \, d\Omega$ (mass matrix)
- $R_i = \int_\Omega N_i \cdot S(\phi_0)(|\nabla\phi^h| - 1) \, d\Omega$ (residual)

### 5.4 Time Discretization

**Forward Euler** (explicit):

$$\hat{\boldsymbol{\phi}}^{n+1} = \hat{\boldsymbol{\phi}}^n - \Delta\tau \, \mathbf{M}^{-1} \mathbf{R}(\hat{\boldsymbol{\phi}}^n)$$

**Backward Euler** (implicit):

$$\mathbf{M} \frac{\hat{\boldsymbol{\phi}}^{n+1} - \hat{\boldsymbol{\phi}}^n}{\Delta\tau} + \mathbf{R}(\hat{\boldsymbol{\phi}}^{n+1}) = 0$$

Requires solving nonlinear system at each step (Newton iteration).

---

## 6. Numerical Discretization

### 6.1 Gradient Magnitude

The term $|\nabla\phi|$ requires special care. For structured grids:

$$|\nabla\phi| = \sqrt{\phi_x^2 + \phi_y^2 + \phi_z^2}$$

### 6.2 Upwind Differencing (Godunov Scheme)

For the hyperbolic nature of the equation, use **upwind differencing**:

$$|\nabla\phi|^2 \approx \begin{cases}
\max(D^{-x}\phi, 0)^2 + \min(D^{+x}\phi, 0)^2 + \ldots & \text{if } S > 0 \\
\min(D^{-x}\phi, 0)^2 + \max(D^{+x}\phi, 0)^2 + \ldots & \text{if } S < 0
\end{cases}$$

where:
- $D^{+x}\phi = (\phi_{i+1} - \phi_i)/h$ (forward difference)
- $D^{-x}\phi = (\phi_i - \phi_{i-1})/h$ (backward difference)

### 6.3 Godunov Hamiltonian

The Godunov numerical flux for $H = S(|\nabla\phi| - 1)$:

$$\hat{H}^{\text{God}} = S \left( \sqrt{\sum_d [H_d^+]^2 + [H_d^-]^2} - 1 \right)$$

where:
$$H_d^+ = \begin{cases} \max(D^{-d}\phi, -D^{+d}\phi, 0) & \text{if } S > 0 \\ \max(D^{+d}\phi, -D^{-d}\phi, 0) & \text{if } S < 0 \end{cases}$$

### 6.4 CFL Condition

For explicit time stepping, stability requires:

$$\Delta\tau \leq \frac{h}{|S| \cdot \max|\nabla\phi|}$$

Since $|\nabla\phi| \to 1$ at steady state and $|S| \leq 1$:

$$\Delta\tau \leq h$$

---

## 7. Convergence and Termination

### 7.1 Convergence Criterion

The iteration converges when the Eikonal property is satisfied:

$$\| |\nabla\phi| - 1 \|_{L^\infty} < \epsilon_{\text{tol}}$$

Or in discrete form:

$$\max_{i} \left| |\nabla\phi|_i - 1 \right| < \epsilon_{\text{tol}}$$

Typical values: $\epsilon_{\text{tol}} = 10^{-3}$ to $10^{-4}$.

### 7.2 Number of Iterations

The pseudo-time required to reach steady state:

$$\tau_{\text{final}} \approx \text{bandwidth}$$

If reinitializing in a narrow band of width $\beta h$:

$$N_{\text{iter}} \approx \frac{\beta h}{\Delta\tau} = \beta$$

Typically 5-20 iterations suffice for a narrow band.

### 7.3 Narrow Band Approach

For efficiency, only reinitialize within distance $\beta$ of the interface:

$$\text{Reinitialize only where } |\phi_0| < \beta \cdot h$$

Common choice: $\beta = 3$ to $6$ mesh cells.

---

## 8. Interface Preservation

### 8.1 The Challenge

The original Sussman-Sethian method can cause **interface drift** (the zero level-set moves).

### 8.2 Sussman-Fatemi Fix

Sussman and Fatemi (1999) introduced a **mass-conserving correction**:

$$\phi^{n+1} = \phi^{n+1}_{\text{reinit}} + \lambda \cdot \delta_\epsilon(\phi_0)$$

where:
- $\delta_\epsilon$ is a regularized delta function
- $\lambda$ is chosen to preserve $\int \delta_\epsilon(\phi) \, d\Omega$

### 8.3 Constrained Reinitialization

Add a Lagrange multiplier to preserve interface:

$$\frac{\partial \phi}{\partial \tau} + S(\phi_0)(|\nabla\phi| - 1) = \lambda \delta(\phi_0)$$

---

## 9. Fast Sweeping Method (Alternative)

### 9.1 Overview

The Fast Sweeping Method (Zhao, 2005) solves the static Eikonal equation directly without pseudo-time iteration.

### 9.2 Algorithm

1. **Initialize**: 
   - $\phi = 0$ on interface nodes
   - $\phi = \infty$ elsewhere

2. **Gauss-Seidel Sweeps**: For each sweep direction $(±1, ±1, ±1)$:
   ```
   for i in sweep_order_x:
       for j in sweep_order_y:
           for k in sweep_order_z:
               φ_new = solve_local_eikonal(φ_neighbors)
               φ[i,j,k] = min(φ[i,j,k], φ_new)
   ```

3. **Repeat** sweeps until convergence (typically 2-4 full sweep cycles)

### 9.3 Local Solver

At each node, solve the quadratic:

$$\max(\phi - \phi_x^{\min}, 0)^2 + \max(\phi - \phi_y^{\min}, 0)^2 + \max(\phi - \phi_z^{\min}, 0)^2 = h^2$$

where $\phi_x^{\min} = \min(\phi_{i-1,j,k}, \phi_{i+1,j,k})$, etc.

### 9.4 Comparison

| Property               | Sussman-Sethian | Fast Sweeping       |
| ---------------------- | --------------- | ------------------- |
| Complexity             | $O(N \cdot M)$  | $O(N)$              |
| Parallelism            | ✅ Natural       | ⚠️ Sequential sweeps |
| Unstructured mesh      | ✅ Yes           | ⚠️ Harder            |
| FEM integration        | ✅ Natural       | ⚠️ External to FEM   |
| Interface preservation | ⚠️ Needs care    | ⚠️ Needs resampling  |

---

## 10. Implementation Summary

### 10.1 Algorithm: FEM-Based Sussman-Sethian

```
Input: φ_current (distorted level-set)
Output: φ_reinit (signed distance function)

1. Store φ₀ = φ_current
2. Compute S_ε(φ₀) using smoothed sign function
3. for n = 1 to max_iterations:
   4. Assemble residual R_i = ∫ N_i · S_ε(φ₀)(|∇φⁿ| - 1) dΩ
   5. Solve: φⁿ⁺¹ = φⁿ - Δτ · M⁻¹ R
   6. if ||∇φⁿ⁺¹| - 1|_∞ < tol:
        break
7. return φ_reinit = φⁿ⁺¹
```

### 10.2 Key Parameters

| Parameter       | Symbol                  | Typical Value | Notes                        |
| --------------- | ----------------------- | ------------- | ---------------------------- |
| Pseudo-timestep | $\Delta\tau$            | $0.5h$        | CFL stability                |
| Smoothing       | $\epsilon$              | $h$           | Sign function regularization |
| Convergence tol | $\epsilon_{\text{tol}}$ | $10^{-4}$     | $\|                          | ∇φ | - 1\|$ |
| Max iterations  | $M$                     | 20-50         | Narrow band: fewer           |
| Bandwidth       | $\beta$                 | $3h$ - $6h$   | Narrow band width            |

---

## 11. References

1. **Osher, S., & Sethian, J. A.** (1988). Fronts propagating with curvature-dependent speed. *Journal of Computational Physics*, 79(1), 12-49.

2. **Sussman, M., Smereka, P., & Osher, S.** (1994). A level set approach for computing solutions to incompressible two-phase flow. *Journal of Computational Physics*, 114(1), 146-159.

3. **Sussman, M., & Fatemi, E.** (1999). An efficient, interface-preserving level set redistancing algorithm. *SIAM Journal on Scientific Computing*, 20(4), 1165-1191.

4. **Peng, D., Merriman, B., Osher, S., Zhao, H., & Kang, M.** (1999). A PDE-based fast local level set method. *Journal of Computational Physics*, 155(2), 410-438.

5. **Zhao, H.** (2005). A fast sweeping method for eikonal equations. *Mathematics of Computation*, 74(250), 603-627.

6. **Sethian, J. A.** (1999). *Level Set Methods and Fast Marching Methods*. Cambridge University Press.

7. **Allaire, G., Jouve, F., & Toader, A. M.** (2004). Structural optimization using sensitivity analysis and a level-set method. *Journal of Computational Physics*, 194(1), 363-393.

---

## Appendix A: Derivation of Weak Form Jacobian

For Newton iteration on the implicit scheme, we need:

$$\frac{\partial R_i}{\partial \hat{\phi}_j} = \int_\Omega N_i \cdot S(\phi_0) \frac{\partial |\nabla\phi|}{\partial \hat{\phi}_j} \, d\Omega$$

where:

$$\frac{\partial |\nabla\phi|}{\partial \hat{\phi}_j} = \frac{\nabla\phi \cdot \nabla N_j}{|\nabla\phi|}$$

Thus:

$$\boxed{J_{ij} = \int_\Omega N_i \cdot S(\phi_0) \frac{\nabla\phi \cdot \nabla N_j}{|\nabla\phi|} \, d\Omega}$$

**Regularization**: When $|\nabla\phi| < \epsilon$, replace with:

$$\frac{1}{|\nabla\phi|} \to \frac{1}{\max(|\nabla\phi|, \epsilon)}$$

---

## Appendix B: 2D Finite Difference Stencil

For structured grid with spacing $h$, explicit Godunov update:

$$\phi_{i,j}^{n+1} = \phi_{i,j}^n - \Delta\tau \cdot S_{i,j} \cdot (G_{i,j} - 1)$$

where:

$$G_{i,j} = \sqrt{[\max(a, -b, 0)]^2 + [\max(c, -d, 0)]^2} \quad \text{if } S > 0$$

with:
- $a = (\phi_{i,j} - \phi_{i-1,j})/h$
- $b = (\phi_{i+1,j} - \phi_{i,j})/h$
- $c = (\phi_{i,j} - \phi_{i,j-1})/h$
- $d = (\phi_{i,j+1} - \phi_{i,j})/h$
