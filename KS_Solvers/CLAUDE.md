# KS_Solvers - Kohn-Sham Eigensolvers

## Purpose
Parallel iterative diagonalization algorithms for the Kohn-Sham Hamiltonian. The Hamiltonian is represented as an operator (H|psi>) rather than a matrix, and these solvers find eigenvalues/eigenvectors iteratively.

## Language
Fortran 90.

## Solver Algorithms

Each solver lives in its own subdirectory:

### Davidson (`Davidson/`)
Block Davidson iterative diagonalization. The most commonly used solver.
- `regterg.f90` - Real (gamma-point) generalized eigenvalue solver
- `cegterg.f90` - Complex (k-point) generalized eigenvalue solver
- `pregterg.f90` / `pcegterg.f90` - Parallel variants

### CG (`CG/`)
Band-by-band Conjugate Gradient minimization.
- `rcgdiagg.f90` - Real CG diagonalization
- `ccgdiagg.f90` - Complex CG diagonalization

### RMM (`RMM/`)
Residual Minimization Method (RMM-DIIS).
- `rrmmdiagg.f90` - Real RMM
- `crmmdiagg.f90` - Complex RMM

### ParO (`ParO/`)
Parallel Orbital minimization.
- `paro_gamma.f90` / `paro_gamma_new.f90` - Gamma-point variants
- `paro_k.f90` / `paro_k_new.f90` - k-point variants

### PPCG (`PPCG_legacy/`)
Preconditioned Conjugate Gradient with GPU support.

### DENSE (`DENSE/`)
Dense matrix diagonalization (full matrix construction and LAPACK solve).
- `rotate_HSpsi_gamma.f90`, `rotate_HSpsi_k.f90` - Subspace rotation
- `rotate_wfc_gamma.f90`, `rotate_wfc_k.f90` - Wavefunction rotation

### JacobiDavidson (`JacobiDavidson/`)
Jacobi-Davidson eigensolver. Processes one eigenvalue at a time with explicit deflation and a TPA (diagonal) preconditioner.
- `rjdsym.f90` - Real (gamma-point) Jacobi-Davidson
- `cjdsym.f90` - Complex (k-point) Jacobi-Davidson

### Davidson_RCI (`Davidson_RCI/`)
Reverse Communication Interface for Davidson.
- `david_rci.f90` - RCI-based Davidson solver

## Interface
- `ks_solver_interfaces.h` - C interface header

## Solver Selection
The solver is selected at runtime via `pw.x` input (`diagonalization` keyword): `'david'`, `'cg'`, `'rmm-davidson'`, `'rmm-paro'`, `'ppcg'`, `'jd'`.

## Dependencies
- LAXlib (matrix operations)
- UtilXlib (MPI, timing)
- ELPA, ScaLAPACK (optional, for subspace diagonalization)

## Build
- Produces `libqeks.a` or CMake target `qe_ks_solvers`
- Each subdirectory has its own Makefile
