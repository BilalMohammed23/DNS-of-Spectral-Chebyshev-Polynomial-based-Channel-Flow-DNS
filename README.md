# Fully Spectral Chebyshev Polynomial-Based Channel Flow DNS

This repository documents the development of a **fully spectral incompressible Navier–Stokes solver** for plane channel flow, built from first principles using Fourier discretisation in the periodic streamwise direction and Chebyshev polynomial expansion in the wall-normal direction.

The project is structured as a **progressive solver development pathway** — starting from a simple scalar advection problem to build intuition for the Chebyshev spectral framework, then extending to a 3D pressure Poisson solver, and finally assembling a fully coupled DNS solver for channel flow. Each stage is self-contained and documented with the underlying numerical theory, derivations, and MATLAB implementation.

---

## Solver Progression

| Stage | Folder | What it solves |
|---|---|---|
| 1 | `Advection Solver` | 1D linear advection — Chebyshev pseudospectral + RK4 |
| 2 | `Pressure Poisson Solver` | 3D Poisson equation — Fourier (x,z) × Chebyshev (y) |
| 3 | `Channel Flow Solver` | 2D incompressible N–S — coupled velocity–pressure DNS |

---

## How to Navigate

Each folder contains its own `README.md` with the full numerical theory, derivations, algorithm, file structure, and results. It is recommended to go through the folders **in order** — the spectral differentiation framework, Chebyshev coefficient recursions, and even–odd tridiagonal solve developed in the earlier stages carry forward directly into the channel flow solver.

If you are already familiar with Chebyshev spectral methods and only want the channel flow formulation, go directly to `Channel Flow Solver`.

---

## Implementation

All solvers are implemented in **MATLAB**. No external toolboxes are required beyond the base MATLAB installation.

---

## Key Numerical Features

- Chebyshev–Gauss–Lobatto collocation in the wall-normal direction
- Fourier pseudospectral discretisation in periodic directions
- 3/2-rule de-aliasing for nonlinear terms
- Adams–Bashforth (explicit) for nonlinear terms, Crank–Nicolson (implicit) for viscous terms
- Global coupled velocity–pressure matrix solve — no pressure boundary condition required
- Even–odd decoupled tridiagonal systems for spectral efficiency

---

## Reference

> Moin, P., & Kim, J. (1980). On the numerical solution of time-dependent viscous incompressible fluid flows involving solid boundaries. Journal of computational physics, 35(3), 381-392.
