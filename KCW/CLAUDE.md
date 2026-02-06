# KCW - Koopmans-Compliant Wannier

## Purpose
Implements Koopmans-compliant functionals in a Wannier function representation. These functionals correct the systematic errors of standard DFT in predicting ionization energies and electron affinities, providing orbital-dependent potentials that satisfy Koopmans' theorem.

## Language
Fortran 90.

## Main Executable
`kcw.x`

## Key Features
- Koopmans screening parameter calculation
- Wannier-based effective Hamiltonian
- DFT reference calculation interface
- Band structure from Koopmans functionals
- Post-processing tools

## Source Organization (`src/`)

### Screening
- `kcw_screened.f90` - Screened Koopmans corrections
- `kcw_screening.f90` - Screening parameter computation

### Hamiltonian
- `kcw_hamiltonian.f90` - Koopmans Hamiltonian construction

### Post-Processing (`PP/`)
- `kcw_bands.f90` - Band structure calculation
- `merge_Umat.f90` - Wannier rotation matrix utilities

## Workflow
1. `pw.x` SCF calculation
2. Wannier90 to generate Wannier functions
3. `kcw.x` for screening parameters and Koopmans Hamiltonian
4. Post-processing for band structures

## Input Format
Namelist-based (see `INPUT_kcw.txt`).

## Examples
5 examples including band structure calculations, screening computations, and symmetry tests.

## Dependencies
- PW (pwlibs), LR_Modules, Wannier90, Modules, FFTXlib, LAXlib

## Build
- `make kcw` from root
- CMake target: `kcw`
