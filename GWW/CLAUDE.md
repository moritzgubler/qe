# GWW - GW with Wannier Functions and Lanczos

## Purpose
Many-body perturbation theory in the GW approximation using ultra-localized Wannier functions and Lanczos chains. Also includes Bethe-Salpeter equation (BSE) for excitonic effects.

## Language
Fortran 90.

## Key Components

### GWW (`gww/`) - Core GW
- `basic_structures.f90` - Fundamental data types
- `input_gw.f90` - GW input parsing
- `green_function.f90`, `go_green.f90` - Green's function evaluation
- `polarization.f90`, `lanczos_polarization.f90` - Polarization via Lanczos
- `self_energy.f90`, `do_self_lanczos*.f90` - Self-energy calculation
- `start_end.f90` - Initialization/finalization

### BSE (`bse/`) - Bethe-Salpeter Equation
- `bse_main.f90` - BSE main driver
- `spectrum.f90` - Optical spectrum from BSE
- `diago_exc.f90` - Excitonic eigenvalue problem

### PW4GWW (`pw4gww/`) - PWscf Interface
Interface routines to extract Wannier functions and matrix elements from PWscf calculations for use in GW.

### MINPACK (`minpack/`)
External minimization library algorithms.

## Workflow
1. `pw.x` SCF calculation
2. `pw4gww.x` to generate Wannier functions and matrix elements
3. `gww.x` for GW quasiparticle energies
4. Optionally `bse.x` for excitonic spectra

## Examples
Methane, Silicon with reference outputs.

## Dependencies
- PW, Modules, LAXlib, FFTXlib

## Build
- `make gwl` from root
- CMake target: `gwl`
