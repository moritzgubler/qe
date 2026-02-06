# dft-d3 - Grimme DFT-D3 Dispersion Correction

## Purpose
Implementation of Grimme's DFT-D3 dispersion correction for van der Waals interactions. Adds semi-empirical London dispersion energy, forces, and stress to DFT calculations.

## Language
Fortran 90.

## Key Features
- DFT-D3 dispersion energy
- Analytical forces and stress tensor
- Multiple damping functions:
  - Zero-damping (original D3)
  - Becke-Johnson (BJ) damping
  - Modified zero-damping (-zerom)
  - Modified BJ damping (-bjm)

## Source Files
- `api.f90` - Public API module
- `core.f90` - Core D3 algorithm
- `common.f90` - Shared data and utilities
- `pars.f90` - D3 parameters (C6 coefficients, etc.)
- `sizes.f90` - Array dimensions
- `dftd3_qe.f90` - QE integration interface
- `test_code.f90` - Standalone test

## Integration
The QE interface (`dftd3_qe.f90`) connects D3 to `pw.x` and `cp.x`. Activated by setting `vdw_corr = 'dft-d3'` in the PW input.

## Origin
Based on the original DFT-D3 code by S. Grimme, adapted for QE. See README for references.

## Dependencies
- None (standalone, but integrated with Modules for QE interface)

## Build
- Produces `libdftd3qe.a` or CMake target `qe_dftd3`
