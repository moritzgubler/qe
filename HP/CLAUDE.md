# HP - Hubbard Parameters from DFPT

## Purpose
Calculates Hubbard U and V parameters self-consistently using density-functional perturbation theory (DFPT). Enables first-principles determination of DFT+U correction parameters without empirical fitting.

## Language
Fortran 90.

## Main Executable
`hp.x`

## Key Features
- Self-consistent Hubbard U parameters
- Extended Hubbard V parameters (inter-site)
- DFPT-based linear response approach
- Support for ultrasoft and PAW pseudopotentials
- Multiple inequivalent atomic sites
- Automatic symmetry handling

## Source Organization (`src/`)

### Setup and Input
- `hp_readin.f90` - Input parsing
- `hp_init.f90` - Initialization
- `hp_setup_q.f90` - q-point setup

### Perturbation and Response
- `hp_dvpsi_pert.f90` - Perturbation potential application
- `hp_calc_chi.f90` - Response function (chi) calculation
- `hp_solve_linear_system.f90` - DFPT linear system solver

### Symmetry
- `hp_psymdvscf.f90` - Symmetrize potential response
- `hp_symdnsq.f90` - Symmetrize occupation matrix response

### Output
- `hp_postproc.f90` - Post-processing of results
- `hp_write_chi.f90` - Write chi matrices
- `hp_write_dnsq.f90` - Write occupation response

### Parallelism
- `hp_allocate_q.f90` - Allocate q-point data
- `hp_clean_q.f90` - Clean up q-point data

## Workflow
1. `pw.x` SCF with `lda_plus_u = .true.`
2. `hp.x` to compute Hubbard parameters
3. Iterate: update U values in `pw.x`, rerun `hp.x` until self-consistent

## Input Format
Namelist-based (see `INPUT_HP.txt`):
- `&INPUTHP` - HP calculation parameters

## Examples
10 examples: LiCoO2, NiO, CrI3, Ni, Ni2MnGa, and others demonstrating various Hubbard parameter calculations.

## Dependencies
- PW (pwlibs), LR_Modules, PHonon (phlibs), Modules, FFTXlib, LAXlib

## Build
- `make hp` from root
- CMake target: `hp`
