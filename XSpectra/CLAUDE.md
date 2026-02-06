# XSpectra - X-ray Absorption Spectra

## Purpose
Calculates X-ray absorption near-edge structure (XANES/NEXAFS) spectra using the Lanczos approach within the PAW formalism. Enables simulation of K, L2,3, and other absorption edges.

## Language
Fortran 90.

## Main Executable
`xspectra.x`

## Key Features
- X-ray absorption spectra (XANES/NEXAFS)
- K-edge and L2,3-edge calculations
- Dipole and quadrupole transition matrix elements
- Core-hole effects (via supercell approach)
- PAW augmentation for transition operators
- DFT+U support
- Cross-section calculations

## Source Organization (`src/`)
- `xspectra.f90` - Main driver
- `init_gipaw_1.f90`, `init_gipaw_2.f90` - GIPAW-style initialization
- `gaunt_mod.f90` - Gaunt coefficients for angular momentum coupling
- `stdout_routines.f90` - Output routines
- `spectra_correction.f90` - Spectral corrections and post-processing
- `read_k_points.f90` - k-point handling

## Workflow
1. `pw.x` SCF with core-hole pseudopotential on absorbing atom
2. `xspectra.x` to compute absorption spectrum

## Input Format
Namelist-based (see `INPUT_XSPECTRA.txt`).

## Examples
Diamond (K-edge), Cu (L2,3-edge), NiO (dipole/quadrupole), SiO2.

## Dependencies
- PW (pwlibs), Modules, upflib

## Build
- `make xspectra` from root
- CMake target: `xspectra`
