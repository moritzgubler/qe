# PWCOND - Ballistic Transport

## Purpose
Calculates ballistic (coherent) transport properties using the Landauer formalism and scattering-state approach. Computes transmission coefficients and conductance for nanostructures (molecular junctions, nanowires, tunnel junctions).

## Language
Fortran 90.

## Main Executable
`pwcond.x`

## Key Features
- Transmission coefficients as function of energy
- Ballistic conductance (Landauer formula)
- Complex band structure
- Scattering state wavefunctions
- Local density of states in transport direction
- Lead-scatterer-lead geometry

## Source Organization (`src/`)
- `condmain.f90` - Main driver program
- `do_cond.f90` - Core conductance calculation
- `scatter_forw.f90` - Forward scattering calculation
- `scatt_states_plot.f90` - Scattering state visualization
- `init_cond.f90` - Initialization
- `init_orbitals.f90` - Orbital setup for transport
- `local.f90` - Local potential handling
- `poten.f90` - Potential setup
- `form_zk.f90` - Complex k-point formation

## Input Format
Namelist-based (see `INPUT_PWCOND.txt`).

## Examples
Al nanowire, Pt chain, Au-CO molecular junction.

## Dependencies
- PW (pwlibs), Modules

## Build
- `make pwcond` from root
- CMake target: `pwcond`
