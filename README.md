# Steady-State Heat Equation Solver (2D FEM)

This project implements a **Finite Element Method (FEM)** solver for the **steady-state heat equation** on a 2D rectangular plate using **triangular elements** and **linear interpolation polynomials**.

The governing PDE is:
\[
-50 \nabla^2 u = 500
\]
with **Dirichlet boundary conditions**
[
u = 100 \text{ on all four sides.}
]

---

## Overview

The code discretizes the domain into triangular elements and assembles the global stiffness matrix and load vector. It then solves for the nodal temperatures using standard FEM formulation.

To validate the results, the computed temperature distribution is compared with MATLAB’s built-in **PDE Toolbox** solution. The comparison gives a close match:

[
| \text{Max}(T) - \text{Max}(U) | = 0.004
]

where

* **T** = maximum temperature from MATLAB’s PDE Toolbox
* **U** = maximum temperature from this implementation

---

## Features

* 2D FEM formulation using **linear triangular elements**
* **Dirichlet boundary conditions** applied on all sides
* Comparison with MATLAB’s **PDE Toolbox**
* Modular code structure for further development

---

## Results

The solver reproduces the temperature field with high accuracy when compared to MATLAB’s reference solution.
The small deviation (~0.004) indicates strong numerical consistency between the implementations.

---

## Future Work

This is the **first installment** of an ongoing project. Upcoming goals include:

* Extending support for **different boundary conditions** (Neumann, Robin, or mixed)
* Allowing **custom source/sink terms** instead of constants
* Enhancing **code generality** and **user configurability**

---

## Notes

This project is primarily educational — a step toward understanding and building general FEM solvers for heat transfer and related PDEs.
