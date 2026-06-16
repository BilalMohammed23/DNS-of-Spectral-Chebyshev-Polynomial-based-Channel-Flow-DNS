# Pressure Poisson Solver — Fourier–Chebyshev Spectral Method

---

## Table of Contents

1. [Problem Statement](#1-problem-statement)
2. [Domain and Boundary Conditions](#2-domain-and-boundary-conditions)
3. [Fourier Transform in x and z](#3-fourier-transform-in-x-and-z)
4. [Chebyshev Expansion in y](#4-chebyshev-expansion-in-y)
5. [Reducing to a 1D ODE per Fourier Mode](#5-reducing-to-a-1d-ode-per-fourier-mode)
6. [Expressing Second Derivatives — Coefficient Recursions](#6-expressing-second-derivatives--coefficient-recursions)
7. [Eliminating $b_n$ and $d_n$ — The Tridiagonal System for $a_n$](#7-eliminating-bn-and-dn--the-tridiagonal-system-for-an)
8. [Even–Odd Decoupling](#8-evenodd-decoupling)
9. [Boundary Condition Enforcement — Tau Method](#9-boundary-condition-enforcement--tau-method)
10. [Algorithm — Full Solution Procedure](#10-algorithm--full-solution-procedure)
11. [File Structure](#11-file-structure)
12. [Parameters](#12-parameters)
13. [Results](#13-results)

---

## 1. Problem Statement

Solve the **3D Poisson equation** for pressure in a channel flow geometry:

$$\frac{\partial^2 P}{\partial x^2} + \frac{\partial^2 P}{\partial y^2} + \frac{\partial^2 P}{\partial z^2} = Q(x, y, z)$$

where $Q$ is a prescribed source term (the RHS arising from the velocity divergence or momentum residual).

The manufactured solution used for verification is:

$$Q = \left(-2 - 8\pi^2(1 - y^2)\right)\cos(x)$$

which corresponds to an exact pressure solution $P = (1 - y^2)\cos(x)$, enabling direct error assessment.

---

## 2. Domain and Boundary Conditions

The domain is a 3D channel:

$$x \in [0, 2\pi], \quad y \in [-1, 1], \quad z \in [0, 2\pi]$$

| Direction | Treatment | Reason |
|---|---|---|
| $x$ | Periodic — Fourier | Streamwise, fully developed flow |
| $z$ | Periodic — Fourier | Spanwise, homogeneous |
| $y$ | Non-periodic — Chebyshev | Wall-normal; bounded by walls at $y = \pm 1$ |

**Wall boundary conditions** (Neumann — normal pressure gradient prescribed):

$$\left.\frac{\partial^2 \hat{P}}{\partial y^2}\right|_{y=-1} = \hat{f}(k_1, k_3), \qquad \left.\frac{\partial^2 \hat{P}}{\partial y^2}\right|_{y=+1} = \hat{g}(k_1, k_3)$$

> At the walls, only the **normal pressure gradient** can be prescribed — this is the physically correct Neumann condition. Here $\hat{f} = \hat{g} = 0$ (homogeneous).
>
> The BCs are imposed on the **second derivative** (not $\hat{P}$ itself) so that Chebyshev polynomials of the same length $N$ can be used throughout.

---

## 3. Fourier Transform in x and z

Since $x$ and $z$ are periodic, expand $P$ as a double Fourier series:

$$P(x, y, z) = \sum_{k_1} \sum_{k_3} \hat{P}(k_1, k_3, y) \, e^{i(k_1 x + k_3 z)}$$

Substituting into the Poisson equation and invoking orthogonality in $x$ and $z$:

$$-k_1^2 \hat{P} + \frac{\partial^2 \hat{P}}{\partial y^2} - k_3^2 \hat{P} = \hat{Q}$$

$$\Rightarrow \quad \frac{\partial^2 \hat{P}}{\partial y^2} - \underbrace{(k_1^2 + k_3^2)}_{\alpha} \hat{P} = \hat{Q}$$

This reduces the 3D PDE to a **family of decoupled 1D ODEs in $y$**, one for each Fourier mode $(k_1, k_3)$, parameterised by:

$$\alpha = k_1^2 + k_3^2$$

> After the Fourier transform we have a purely **algebraic** structure in $x$ and $z$ — what remains is an ODE in $y$ only.

**Implemented in:** `main_PPE.m` — FFT over dimensions 1 and 3, wavenumber array $\alpha(k_1, k_3)$ assembled in a double loop.

---

## 4. Chebyshev Expansion in y

For each Fourier mode, $\hat{P}(y)$ is expanded in Chebyshev polynomials on $y \in [-1, 1]$:

$$\hat{P}(y) = \sum_{n=0}^{N} a_n T_n(y)$$

The Gauss–Lobatto nodes and the $\theta$ substitution $y = \cos\theta$ are used (identical to the advection solver):

$$y_j = \cos\left(\frac{j\pi}{N}\right), \quad j = 0, 1, \ldots, N$$

The derivatives expand as:

$$\frac{\partial \hat{P}}{\partial y} = \sum_{n=0}^{N} b_n T_n(y), \quad b_N = 0$$

$$\frac{\partial^2 \hat{P}}{\partial y^2} = \sum_{n=0}^{N} d_n T_n(y), \quad d_N = d_{N-1} = 0$$

> Each differentiation reduces the polynomial degree by 1, hence $b_N = 0$ and $d_N = d_{N-1} = 0$.

**Forward Chebyshev transform of $\hat{Q}(y)$ for each mode:** `qn_cheb_coeff_RK4.m`

---

## 5. Reducing to a 1D ODE per Fourier Mode

Substituting the Chebyshev expansions into the Fourier-transformed Poisson equation:

$$\sum_{n=0}^{N} d_n T_n(y) - \alpha \sum_{n=0}^{N} a_n T_n(y) = \sum_{n=0}^{N} q_n T_n(y)$$

Invoking orthogonality of $T_n(y)$, coefficient by coefficient:

$$\boxed{d_n - \alpha \, a_n = q_n}, \quad n = 0, 1, 2, \ldots, N$$

This gives $N+1$ governing equations. We now have **3N unknowns** in total:

| Coefficient set | Count | Equations available |
|---|---|---|
| $a_n$ | $N+1$ | $N+1$ governing equations |
| $b_n$ | $N$ | $N$ recursion equations |
| $d_n$ | $N-1$ | $N-1$ recursion equations |

Total: $3N$ unknowns, $3N$ equations. However, the two wall BCs are not yet included — so **2 equations (for $n=0$ and $n=N$) must be dropped** and replaced by the BCs.

---

## 6. Expressing Second Derivatives — Coefficient Recursions

From the derivative recursion (same as in the advection solver):

$$2n \, b_n = c_{n-1} \, d_{n-1} - d_{n+1}$$

$$2n \, a_n = c_{n-1} \, b_{n-1} - b_{n+1}$$

where $c_n = 2$ if $n = 0$, else $c_n = 1$.

From the PDE $d_n = \alpha a_n + q_n$, substituting into the $b_n$ recursion:

$$b_n = \frac{1}{2n}\left[c_{n-1} \, q_{n-1} - q_{n+1} + \alpha\left(c_{n-1} a_{n-1} - a_{n+1}\right)\right]$$

Substituting $b_n$ into the $a_n$ recursion and collecting terms yields the final tridiagonal equation for $a_n$:

$$\frac{\alpha \, c_{n-1} c_{n-2}}{2(n-1)} a_{n-2} - \left(\frac{\alpha c_{n-1}}{2(n-1)} + \frac{\alpha c_n}{2(n+1)} + 2n\right) a_n + \frac{\alpha}{2(n+1)} a_{n+2} = \text{RHS}(q_n)$$

where the right-hand side is:

$$\text{RHS}(q_n) = \frac{1}{2(n+1)}\left(c_n q_n - q_{n+2}\right) - \frac{1}{2(n-1)}\left(c_{n-1} c_{n-2} q_{n-2} - c_{n-1} q_n\right)$$

> **Key structural observation:** The equation for $a_n$ couples only to $a_{n-2}$ and $a_{n+2}$ — **even and odd modes are completely decoupled**. This means a single $(N+1)\times(N+1)$ system splits into two independent $(N/2)\times(N/2)$ tridiagonal systems.

---

## 7. Eliminating $b_n$ and $d_n$ — The Tridiagonal System for $a_n$

After eliminating $b_n$ and $d_n$ analytically, the system for $a_n$ takes the block-tridiagonal form:

$$\mathbf{M}_1 \, \mathbf{a} = \mathbf{M}_2$$

where $\mathbf{M}_1$ is a **tridiagonal matrix** (with entries coupling $a_{n-2}$, $a_n$, $a_{n+2}$ on sub-, main, and super-diagonal respectively) and $\mathbf{M}_2$ is the RHS vector built from $q_n$ and the boundary values.

A representative $16\times16$ example of the even-mode matrix has the structure:

```
Row 0  (BC):  [ 0    8    32   72  ...  512  ]
Row 1  (n=2): [ 2α/8  -(1+α/16)  α/24  0   0  ... ]
...
```

> The off-diagonal elements can be large, requiring careful pivoting or reordering. The even–odd split reduces matrix size and improves conditioning.

**Implemented in:** `tridiag_struc_even.m`, `tridiag_struc_odd.m`, `solve_tridiag_an.m`

---

## 8. Even–Odd Decoupling

Because the coupling in the $a_n$ equation is between $a_n$ and $a_{n \pm 2}$, the system naturally splits:

**Even system** — involves $a_0, a_2, a_4, \ldots, a_{N-2}$:

$$\mathbf{M}_1^{\text{even}} \, [a_0, a_2, \ldots, a_{N-2}]^T = \mathbf{M}_2^{\text{even}}$$

**Odd system** — involves $a_1, a_3, a_5, \ldots, a_{N-1}$:

$$\mathbf{M}_1^{\text{odd}} \, [a_1, a_3, \ldots, a_{N-1}]^T = \mathbf{M}_2^{\text{odd}}$$

Each system is of size $N/2 \times N/2$, solved independently and then interleaved to recover the full coefficient vector $a_n$.

**Implemented in:** `solve_tridiag_an.m` — solves even and odd systems separately via `m1_even \ m2_even` and `m1_odd \ m2_odd`, then assembles $a_n$.

---

## 9. Boundary Condition Enforcement — Tau Method

The BCs on the second derivative at $y = \pm 1$:

$$\left.\frac{\partial^2 \hat{P}}{\partial y^2}\right|_{y=-1} = \hat{f} \quad \Rightarrow \quad \sum_{n=0}^{N} b_n \, T_n(-1) = \hat{f} \quad \Rightarrow \quad \sum_{n=0}^{N} b_n (-1)^n = \hat{f}; \quad b_N = 0$$

$$\left.\frac{\partial^2 \hat{P}}{\partial y^2}\right|_{y=+1} = \hat{g} \quad \Rightarrow \quad \sum_{n=0}^{N} b_n \, T_n(+1) = \hat{g} \quad \Rightarrow \quad \sum_{n=0}^{N} b_n (1)^n = \hat{g}; \quad b_N = 0$$

Since $T_n(\pm 1) = (\pm 1)^n$, these sum to:

- **Odd system BC:** $\hat{f} + \hat{g}$ (symmetric combination, $n=1$ row replaced)
- **Even system BC:** $\hat{f} - \hat{g}$ (anti-symmetric combination, $n=0$ row replaced)

The $n = 0$ equation (even) and $n = 1$ equation (odd) in each tridiagonal system are **replaced** by these BC conditions — this is the Tau method applied to the second-derivative BCs.

In `main_PPE.m`: `ghat = 0; fhat = 0;` (homogeneous Neumann BCs).

---

## 10. Algorithm — Full Solution Procedure

```
1. Set up 3D grid
      x, z : uniform Fourier grid,  Nx points on [0, 2π]
      y     : Chebyshev–Gauss–Lobatto nodes on [-1, 1]

2. Compute source term Q(x, y, z) on the 3D grid

3. FFT Q in x and z  →  Q̂(k₁, y, k₃)
      Q_fft = fft(fft(Q, [], 1), [], 3) / (Nx*Nx)

4. Compute wavenumber array  α(k₁, k₃) = k₁² + k₃²

5. For each Fourier mode (k₁, k₃):
      a. Chebyshev transform Q̂(:,y,:) in y  →  q_n        [qn_cheb_coeff_RK4]
      b. Solve even tridiagonal system  →  a_n (even)       [tridiag_struc_even]
      c. Solve odd  tridiagonal system  →  a_n (odd)        [tridiag_struc_odd]
      d. Interleave even + odd  →  full a_n                 [solve_tridiag_an]
      e. Reconstruct P̂(y) from a_n via Clenshaw            [cheb_eval_series_RK4]

6. Inverse FFT P̂ in x and z  →  P(x, y, z)
      P = ifft(ifft(P̂, [], 1), [], 3, 'symmetric')

7. Visualise pressure slice at z = z_mid
```

---

## 11. File Structure

| File | Role | Called by |
|---|---|---|
| `main_PPE.m` | Driver: grid setup, FFT, mode loop, IFFT, plot | — |
| `qn_cheb_coeff_RK4.m` | Forward DCT in $y$: $\hat{Q}(y) \to q_n$ for each mode | `main_PPE.m` |
| `solve_tridiag_an.m` | Assembles and solves even + odd tridiagonal systems for $a_n$ | `main_PPE.m` |
| `tridiag_struc_even.m` | Builds $\mathbf{M}_1^{\text{even}}$ and $\mathbf{M}_2^{\text{even}}$ with BC in row 0 | `solve_tridiag_an.m` |
| `tridiag_struc_odd.m` | Builds $\mathbf{M}_1^{\text{odd}}$ and $\mathbf{M}_2^{\text{odd}}$ with BC in row 1 | `solve_tridiag_an.m` |
| `cheb_eval_series_RK4.m` | Inverse Chebyshev: $a_n \to \hat{P}(y)$ via Clenshaw recurrence | `main_PPE.m` |
| `misc.m` | Utility script — tests even/odd index extraction from arrays | — |

### Call Graph

```
main_PPE.m
│
├── FFT in x, z  →  Q̂(k₁, y, k₃)
│
└── Loop over (k₁, k₃) modes:
    ├── qn_cheb_coeff_RK4.m      ← Q̂(y) → q_n  (DCT in y)
    ├── solve_tridiag_an.m       ← q_n + α + BCs → a_n
    │   ├── tridiag_struc_even.m    ← builds M1_even, M2_even
    │   └── tridiag_struc_odd.m     ← builds M1_odd,  M2_odd
    └── cheb_eval_series_RK4.m   ← a_n → P̂(y)  (Clenshaw)
│
└── IFFT in x, z  →  P(x, y, z)
```

---

## 12. Parameters

| Parameter | Symbol | Value | Description |
|---|---|---|---|
| Grid points | $N_x$ | `128` | Points in $x$ and $z$ (Fourier) |
| Domain length | $L$ | $2\pi$ | Periodic box length in $x$, $z$ |
| Chebyshev points | $N$ | $N_x - 1 = 127$ | Chebyshev modes in $y$ |
| Wavenumber range | $k_1, k_3$ | $[0, N_x/2-1, -N_x/2, \ldots, -1]$ | Standard FFT ordering |
| Wall BC (bottom) | $\hat{f}$ | `0` | Neumann BC at $y = -1$ |
| Wall BC (top) | $\hat{g}$ | `0` | Neumann BC at $y = +1$ |

---

## 13. Results

### Pressure Field — 3D Fourier–Chebyshev Poisson Solver ($N_x = 128$)

Pressure slice at the $z$-midplane, showing the reconstructed pressure field from the manufactured solution test case $P = (1-y^2)\cos(x)$.

![3D Pressure Poisson Solver result](Results/3D_Pressure_Poisson_solver.png)

> The spectral method recovers the smooth pressure field with exponential accuracy. The Fourier–Chebyshev approach is exact for smooth periodic-in-$x$, $z$ and non-periodic-in-$y$ problems, consistent with the channel flow geometry.

---

*Solver: 3D Pressure Poisson Equation — Fourier (x, z) × Chebyshev (y) spectral decomposition with even–odd tridiagonal solve and Tau BC enforcement.*
