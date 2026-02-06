# LR_Modules - Linear Response Modules

## Purpose
Shared library of core algorithms for linear-response DFPT calculations. Used by PHonon, HP, TDDFPT, KCW, and other packages that need perturbation-theory response functions.

## Language
Fortran 90.

## Key Components

### Response Kernels
- `response_kernels.f90` - Exchange-correlation and Hartree response kernels
- `dfpt_kernels.f90` - DFPT-specific kernel computations

### DFPT Core
- `dfpt_type.f90` - DFPT data types
- `dfpt_tetra_mod.f90` - Tetrahedron method for DFPT
- `apply_dpot.f90` - Apply perturbation potential to wavefunctions
- `dv_of_drho.f90` - Potential response from density response

### Linear Solvers
- `cgsolve_all.f90` - Conjugate gradient solver (real)
- `ccgsolve_all.f90` - Conjugate gradient solver (complex)

### Symmetry in Linear Response
- `symdvscf.f90` - Symmetrize SCF potential response
- `sym_def.f90` - Symmetry definitions for response
- `sym_dns.f90` - Symmetrize density matrix response

### Perturbation Handling
- `apply_dpot.f90` - Apply perturbation potential
- `addusdbec.f90` - Ultrasoft augmentation contribution to response
- `orthogonalize.f90` - Orthogonalization of response wavefunctions

### Special Functions
- Lanczos methods for iterative response
- Orthogonalization and projection routines

## Dependencies
- Modules, FFTXlib, LAXlib, PW (pwlibs)

## Build
- Produces `libqelrmod.a` or CMake target `qe_lr_modules`
- Built via `make lrmods` from root
