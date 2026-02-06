# TDDFPT - Time-Dependent DFPT

## Purpose
Linear-response time-dependent DFPT calculations for spectroscopy: optical absorption, electron energy loss spectroscopy (EELS), inelastic X-ray scattering, and magnetic excitations (magnons).

## Language
Fortran 90.

## Main Executables
- `turbo_lanczos.x` - Lanczos-based linear response
- `turbo_davidson.x` - Davidson-based linear response
- `turbo_spectrum.x` - Spectrum post-processing
- `turbo_magnon.x` - Magnetic excitation spectra (spin waves)
- `turbo_eels.x` - Electron energy loss spectra

## Key Features
- Optical absorption spectra (UV-Vis)
- EELS (bulk and surface)
- Inelastic X-ray scattering
- Resonant inelastic X-ray scattering
- Molecular NEXAFS (near-edge X-ray absorption)
- Spin-wave (magnon) spectra
- Lorentzian broadening
- Continuum solvation model (PCM) support

## Source Organization (`src/`)

### Lanczos Solvers
- `lanczos_nonhermitian.f90` - Non-Hermitian Lanczos (general case)
- `lanczos_pseudohermitian.f90` - Pseudo-Hermitian Lanczos (optimized)

### Spectrum Processing
- `spectrum.f90` - Spectrum computation from Lanczos coefficients
- `absorption.f90` - Absorption spectrum specific routines
- `print_spectrum.f90` - Output formatting

### Tools
- `ColorCalculator/` - Java tool for color visualization from absorption spectra

## Workflow
1. `pw.x` SCF calculation
2. `turbo_lanczos.x` or `turbo_davidson.x` for linear response
3. `turbo_spectrum.x` to extract spectra

## Input Format
Namelist-based. Multiple INPUT documentation files for each executable.

## Examples
CH4 absorption, benzene, Cu L2,3 edges, NiO magnetic, SiO2, diamond.

## Dependencies
- PW, LR_Modules, Modules, FFTXlib, LAXlib

## Build
- `make tddfpt` from root
- CMake target: `tddfpt`
