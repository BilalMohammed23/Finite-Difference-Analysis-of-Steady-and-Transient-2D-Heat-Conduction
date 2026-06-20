# 2D Heat Conduction Solver — Steady and Transient State

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Domain and Boundary Conditions](#2-domain-and-boundary-conditions)
3. [Governing Equations](#3-governing-equations)
4. [Steady-State Discretisation](#4-steady-state-discretisation)
5. [Transient Discretisation — Explicit Scheme](#5-transient-discretisation--explicit-scheme)
6. [Transient Discretisation — Implicit Scheme](#6-transient-discretisation--implicit-scheme)
7. [Iterative Solvers](#7-iterative-solvers)
8. [Successive Over-Relaxation (SOR)](#8-successive-over-relaxation-sor)
9. [Convergence Criterion](#9-convergence-criterion)
10. [File Structure](#10-file-structure)
11. [Parameters](#11-parameters)

---

## 1. Problem Statement

Solve the **2D heat conduction equation** on a unit square domain under both **steady-state** and **transient** (unsteady) conditions. The problem is solved using finite difference discretisation with five iterative solvers:

- Jacobi
- Gauss–Seidel
- SOR over Jacobi
- SOR over Gauss–Seidel
- Explicit (transient only)

Each solver is compared in terms of number of iterations and wall-clock time to convergence.

---

## 2. Domain and Boundary Conditions

$$x \in [0, 1], \quad y \in [0, 1]$$

A uniform Cartesian grid of $n \times n$ nodes is used (collocated, cell-vertex arrangement).

**Dirichlet boundary conditions** (fixed temperatures on all four walls):

| Boundary | Temperature (K) |
|---|---|
| Bottom ($y = 0$) | $T = 900$ |
| Top ($y = 1$) | $T = 600$ |
| Left ($x = 0$) | $T = 400$ |
| Right ($x = 1$) | $T = 800$ |

**Corner nodes** are set to the average of the two adjacent wall temperatures:

$$T_{\text{corner}} = \frac{T_{\text{wall 1}} + T_{\text{wall 2}}}{2}$$

e.g., bottom-right corner: $T = (900 + 800)/2 = 850$ K.

**Initial condition (transient):** Interior nodes initialised to $T = 300$ K; steady-state interior initialised to $T = 1$ K (uniform, arbitrary — iterates to solution).

---

## 3. Governing Equations

### Steady-State (Laplace Equation)

For steady-state heat conduction with no internal heat generation:

$$\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2} = 0$$

### Transient (Diffusion Equation)

For unsteady heat conduction:

$$\frac{\partial T}{\partial t} = \alpha \left(\frac{\partial^2 T}{\partial x^2} + \frac{\partial^2 T}{\partial y^2}\right)$$

where $\alpha$ is the thermal diffusivity. The solution is marched in time until it reaches the same steady-state temperature distribution.

---

## 4. Steady-State Discretisation

Applying **2nd-order central finite differences** to the Laplace equation at interior node $(i, j)$:

$$\frac{T_{i-1,j} - 2T_{i,j} + T_{i+1,j}}{\Delta x^2} + \frac{T_{i,j-1} - 2T_{i,j} + T_{i,j+1}}{\Delta y^2} = 0$$

Solving for $T_{i,j}$:

$$\boxed{T_{i,j} = \frac{1}{k \Delta x^2}\left(T_{i-1,j} + T_{i+1,j}\right) + \frac{1}{k \Delta y^2}\left(T_{i,j-1} + T_{i,j+1}\right)}$$

where the coefficient $k$ is:

$$k = \frac{2(\Delta x^2 + \Delta y^2)}{\Delta x^2 \, \Delta y^2}$$

> This is the discrete harmonic mean stencil — each neighbouring node contributes inversely proportional to the square of the grid spacing in that direction. For a uniform grid ($\Delta x = \Delta y$), this reduces to the standard 5-point Laplacian average: $T_{i,j} = (T_{i-1,j} + T_{i+1,j} + T_{i,j-1} + T_{i,j+1})/4$.

---

## 5. Transient Discretisation — Explicit Scheme

Discretising the diffusion equation in time using a **forward Euler** (explicit) scheme:

$$\frac{T_{i,j}^{n+1} - T_{i,j}^n}{\Delta t} = \alpha \left[\frac{T_{i-1,j}^n - 2T_{i,j}^n + T_{i+1,j}^n}{\Delta x^2} + \frac{T_{i,j-1}^n - 2T_{i,j}^n + T_{i,j+1}^n}{\Delta y^2}\right]$$

Solving for $T_{i,j}^{n+1}$:

$$\boxed{T_{i,j}^{n+1} = T_{i,j}^n + k_1 \underbrace{\left(T_{i-1,j}^n - 2T_{i,j}^n + T_{i+1,j}^n\right)}_{H} + k_2 \underbrace{\left(T_{i,j-1}^n - 2T_{i,j}^n + T_{i,j+1}^n\right)}_{V}}$$

where:

$$k_1 = \frac{\alpha \, \Delta t}{\Delta x^2}, \qquad k_2 = \frac{\alpha \, \Delta t}{\Delta y^2}$$

$H$ and $V$ are the discrete second derivatives in $x$ and $y$ respectively — directly the second-difference stencils evaluated at the old time level.

**Stability constraint (von Neumann):**

$$\Delta t \leq \frac{1}{2\alpha\left(1/\Delta x^2 + 1/\Delta y^2\right)}$$

Violating this causes the explicit scheme to diverge. The implicit schemes below are unconditionally stable.

**Implemented in:** `transient__explicit_1.m`

---

## 6. Transient Discretisation — Implicit Scheme

Discretising with a **backward Euler** (fully implicit) scheme at time level $n+1$:

$$\frac{T_{i,j}^{n+1} - T_{i,j}^n}{\Delta t} = \alpha \left[\frac{T_{i-1,j}^{n+1} - 2T_{i,j}^{n+1} + T_{i+1,j}^{n+1}}{\Delta x^2} + \frac{T_{i,j-1}^{n+1} - 2T_{i,j}^{n+1} + T_{i,j+1}^{n+1}}{\Delta y^2}\right]$$

Collecting all $n+1$ terms to the left and rearranging:

$$T_{i,j}^{n+1}\underbrace{\left(1 + 2k_1 + 2k_2\right)}_{1/\text{term}_1} = T_{i,j}^n + k_1\left(T_{i-1,j}^{n+1} + T_{i+1,j}^{n+1}\right) + k_2\left(T_{i,j-1}^{n+1} + T_{i,j+1}^{n+1}\right)$$

Solving for $T_{i,j}^{n+1}$:

$$\boxed{T_{i,j}^{n+1} = \underbrace{\frac{1}{1 + 2k_1 + 2k_2}}_{\text{term}_1} \cdot T_{i,j}^n + \underbrace{\frac{k_1}{1 + 2k_1 + 2k_2}}_{\text{term}_2} \cdot H + \underbrace{\frac{k_2}{1 + 2k_1 + 2k_2}}_{\text{term}_3} \cdot V}$$

where:

$$\text{term}_1 = \frac{1}{1 + 2k_1 + 2k_2}, \quad \text{term}_2 = k_1 \cdot \text{term}_1, \quad \text{term}_3 = k_2 \cdot \text{term}_1$$

and $H$, $V$ are the sum of neighbouring values at the new time level $n+1$ (or old level $n$ depending on the solver variant). Since $T^{n+1}$ appears on both sides, this requires an iterative solve at each time step — which is where the Jacobi/Gauss–Seidel/SOR solvers come in.

**Implemented in:** `transient_jac.m`, `transient_gauss.m`, `transient_sor_jac.m`, `transient_sor_gauss.m`

---

## 7. Iterative Solvers

All solvers iterate until the convergence criterion is met (see Section 9). The key difference between them is **which values of $T$ are used for the neighbouring nodes** during the update sweep.

### 7.1 Jacobi Method
 
Uses exclusively **old iteration** values for all neighbours. The entire field is updated simultaneously — the new value at $(i,j)$ uses only $T^{(k)}$ (the previous iterate) for all neighbours:
 
**Steady state:**
$$T_{i,j}^{(k+1)} = \frac{1}{k \Delta x^2}\left(T_{i-1,j}^{(k)} + T_{i+1,j}^{(k)}\right) + \frac{1}{k \Delta y^2}\left(T_{i,j+1}^{(k)} + T_{i,j-1}^{(k)}\right)$$
 
**Transient (implicit):**
$$T_{i,j}^{(k+1)} = \text{term}_1 \cdot T_{\text{initial}} + \text{term}_2 \cdot \left(T_{i-1,j}^{(k)} + T_{i+1,j}^{(k)}\right) + \text{term}_3 \cdot \left(T_{i,j-1}^{(k)} + T_{i,j+1}^{(k)}\right)$$
 
> Jacobi converges more slowly because it never uses updated neighbours within the same sweep. It requires more iterations than Gauss–Seidel but is straightforwardly parallelisable.
 
**Implemented in:** `steady_jac.m`, `transient_jac.m`
 
---
 
### 7.2 Gauss–Seidel Method
 
Uses **immediately updated** values as they become available within the same sweep. When computing $T_{i,j}$, the neighbours $T_{i-1,j}$ and $T_{i,j-1}$ (already updated in the current sweep) are used directly:
 
**Steady state:**
$$T_{i,j}^{(k+1)} = \frac{1}{k \Delta x^2}\left(T_{i-1,j}^{(k+1)} + T_{i+1,j}^{(k)}\right) + \frac{1}{k \Delta y^2}\left(T_{i,j+1}^{(k)} + T_{i,j-1}^{(k+1)}\right)$$
 
**Transient (implicit):**
$$T_{i,j}^{(k+1)} = \text{term}_1 \cdot T_{\text{initial}} + \text{term}_2 \cdot \left(T_{i-1,j}^{(k+1)} + T_{i+1,j}^{(k)}\right) + \text{term}_3 \cdot \left(T_{i,j-1}^{(k+1)} + T_{i,j+1}^{(k)}\right)$$
 
> Gauss–Seidel typically converges in roughly half the iterations of Jacobi because in-sweep updates propagate information faster. The in-place update is visible in the code: `T(i-1,j)` and `T(i,j-1)` already hold the new-iterate values when the $(i,j)$ update is computed.
 
**Implemented in:** `steady_gauss.m`, `transient_gauss.m`
 
---

## 8. Successive Over-Relaxation (SOR)

SOR accelerates convergence by **blending the iterative update with the previous value** using a relaxation factor $\alpha$:

$$T_{i,j}^{\text{SOR}} = (1 - \alpha)\,T_{i,j}^{(k)} + \alpha\,T_{i,j}^{*}$$

where $T_{i,j}^{*}$ is the standard Jacobi or Gauss–Seidel update, and $\alpha$ is the over-relaxation factor.

| $\alpha$ value | Behaviour |
|---|---|
| $\alpha = 1$ | Equivalent to standard Jacobi / Gauss–Seidel |
| $1 < \alpha < 2$ | Over-relaxation — accelerates convergence |
| $\alpha > 2$ | Unstable |
| $0 < \alpha < 1$ | Under-relaxation — slower but more stable |

> **Note from code:** `SOR_Jacobi` is noted to be unstable for $\alpha > 1$ in this implementation. `SOR_Gauss–Seidel` with $\alpha = 1.2$ (the value used here) typically converges faster than pure Gauss–Seidel.

**Implemented in:** `steady_sor_jac.m`, `steady_sor_gauss.m`, `transient_sor_jac.m`, `transient_sor_gauss.m`

---

## 9. Convergence Criterion

All iterative solvers use an **absolute maximum error** criterion:

$$\varepsilon = \max_{i,j} \left|T_{i,j}^{(k+1)} - T_{i,j}^{(k)}\right| < \text{tol}$$

where $\text{tol} = 10^{-4}$ K. The `max(max(...))` operation in MATLAB first takes the column-wise maximum of the absolute difference matrix, then the overall maximum — giving the single largest change across the entire field. Iteration continues until this falls below the tolerance.

For the **transient solvers**, this convergence check is applied at every time step — the inner iteration loop runs until the implicit equation converges at that time level, then the time is advanced and the process repeats.

---

## 10. File Structure

| File | Solver type | Problem | Called by |
|---|---|---|---|
| `main.m` | Driver | Steady-state | — |
| `main_unsteady.m` | Driver | Transient | — |
| `steady_jac.m` | Jacobi | Steady | `main.m` |
| `steady_gauss.m` | Gauss–Seidel | Steady | `main.m` |
| `steady_sor_jac.m` | SOR–Jacobi | Steady | `main.m` |
| `steady_sor_gauss.m` | SOR–Gauss–Seidel | Steady | `main.m` |
| `transient_jac.m` | Implicit Jacobi | Transient | `main_unsteady.m` |
| `transient_gauss.m` | Implicit Gauss–Seidel | Transient | `main_unsteady.m` |
| `transient_sor_jac.m` | Implicit SOR–Jacobi | Transient | `main_unsteady.m` |
| `transient_sor_gauss.m` | Implicit SOR–Gauss–Seidel | Transient | `main_unsteady.m` |
| `transient__explicit_1.m` | Explicit forward Euler | Transient | `main_unsteady.m` |

### Call Graph

```
main.m  (steady)
│
├── steady_jac.m          ← option 1: Jacobi
├── steady_gauss.m        ← option 2: Gauss–Seidel
├── steady_sor_gauss.m    ← option 3: SOR–Gauss–Seidel
└── steady_sor_jac.m      ← option 4: SOR–Jacobi

main_unsteady.m  (transient)
│
├── transient_jac.m            ← option 1: Implicit Jacobi
├── transient_gauss.m          ← option 2: Implicit Gauss–Seidel
├── transient_sor_gauss.m      ← option 3: Implicit SOR–Gauss–Seidel
├── transient_sor_jac.m        ← option 4: Implicit SOR–Jacobi
└── transient__explicit_1.m    ← option 5: Explicit forward Euler
```

---

## 11. Parameters

### Steady State (`main.m`)

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Grid nodes | $n$ | User input (≤ 30) | Square grid $n \times n$ |
| Domain | — | $[0,1]^2$ | Unit square |
| Tolerance | tol | `1e-4` | Convergence criterion |
| Relaxation factor | $\alpha$ | `1.2` | SOR acceleration (GS); unstable > 1 for Jacobi |
| Stiffness coefficient | $k$ | $2(\Delta x^2 + \Delta y^2)/(\Delta x^2 \Delta y^2)$ | Steady-state discrete Laplacian coefficient |

### Transient (`main_unsteady.m`)

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Grid nodes | $n$ | User input (≤ 30) | Square grid $n \times n$ |
| Thermal diffusivity | $\alpha_d$ | `1.4` | Diffusion coefficient |
| Time step | $\Delta t$ | `1e-3` | Time increment |
| Tolerance | tol | `1e-4` | Per-time-step convergence |
| Relaxation factor | $\alpha$ | `1.2` | SOR factor |
| Time steps | — | `1400` | Total time steps |
| $k_1$ | $k_1$ | $\alpha_d \Delta t / \Delta x^2$ | x-direction diffusion coefficient |
| $k_2$ | $k_2$ | $\alpha_d \Delta t / \Delta y^2$ | y-direction diffusion coefficient |
| term₁ | — | $1/(1 + 2k_1 + 2k_2)$ | Implicit scheme central coefficient |
| term₂ | — | $k_1 \cdot \text{term}_1$ | Implicit scheme x-neighbour weight |
| term₃ | — | $k_2 \cdot \text{term}_1$ | Implicit scheme y-neighbour weight |

---

### Reference

> Incropera, F. P., DeWitt, D. P., Bergman, T. L., & Lavine, A. S. (2007). *Fundamentals of Heat and Mass Transfer* (7th ed.). Wiley.
>
> Anderson, J. D. (1995). *Computational Fluid Dynamics: The Basics with Applications.* McGraw-Hill. (Chapter 4 — Elliptic PDEs and iterative methods)

---

*Solver: 2D Heat Conduction — Steady (Laplace) and Transient (Diffusion), Jacobi / Gauss–Seidel / SOR iterative methods, Explicit and Implicit time discretisation. Implemented in MATLAB.*
