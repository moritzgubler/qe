# PHonon - Vibrational and Dielectric Properties

## Purpose
Calculates phonon frequencies, dielectric tensors, Born effective charges, Raman tensors, and electron-phonon coupling using Density-Functional Perturbation Theory (DFPT). One of the most widely used QE packages.

## Language
Fortran 90.

## Main Executable
`ph.x` - Phonon calculation program.

## Key Features
- Phonon frequencies and eigenvectors at arbitrary q-points
- Dynamical matrices and interatomic force constants
- Dielectric constant (electronic and ionic contributions)
- Born effective charges
- Raman tensors (non-resonant)
- Electron-phonon coupling coefficients
- Infrared and Raman spectra
- Third-order derivatives (anharmonic)
- Support for metals (with smearing) and insulators

## Source Organization

### PH (`PH/`) - Main phonon calculation
- `phq_readin.f90` - Input parsing
- `phq_setup.f90` - Setup perturbations and symmetry
- `solve_linter.f90` - Solve linear system for DFPT response
- `dynmat0.f90` / `dynmat_us.f90` - Dynamical matrix computation
- `elphon.f90` - Electron-phonon coupling
- `drho.f90` - Density response calculation
- `dvpsi_e.f90` - Electric field perturbation
- `q_points_wannier.f90` - q-point grid generation

### Gamma (`Gamma/`) - Gamma-point specialized routines
Optimized routines for zone-center (Gamma-point) phonons only.

### FD (`FD/`) - Finite displacement method
Alternative frozen-phonon approach using finite atomic displacements.

## Additional Executables
- `dynmat.x` - Dynamical matrix analysis and IR/Raman spectra
- `q2r.x` - q-space to real-space force constants
- `matdyn.x` - Phonon dispersion from force constants
- `fqha.x` - Quasi-harmonic approximation
- `lambda.x` - Electron-phonon coupling lambda
- `alpha2f.x` - Eliashberg function

## Workflow
1. Run `pw.x` SCF calculation
2. Run `ph.x` for DFPT at desired q-points
3. Run `q2r.x` to get real-space force constants
4. Run `matdyn.x` for phonon dispersion and DOS

## Input Format
Namelist-based (see `INPUT_PH.txt`):
- `&INPUTPH` - Main phonon input (q-point, convergence, etc.)

## Documentation
- `Doc/developer_man.pdf` - Developer manual
- `INPUT_PH.txt` / `INPUT_PH.def` - Input documentation

## Examples
Multiple examples: basic phonons, electron-phonon, image parallelism, grid recovery, partial calculations, tetrahedron method, etc.

## Dependencies
- PW (pwlibs), LR_Modules, Modules, FFTXlib, LAXlib

## Build
- `make ph` from root
- CMake target: `ph`
