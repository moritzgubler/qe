# PIOUD - Path Integral Molecular Dynamics

## Purpose
Path integral molecular dynamics (PIMD) using the PIOUD algorithm. Enables quantum nuclear effects in ab initio simulations through ring-polymer representation of atomic nuclei.

## Language
Fortran 90.

## Key Features
- Ring polymer MD (RPMD) for quantum nuclear effects
- Thermostat integration (Nose-Hoover, Langevin)
- Centroid molecular dynamics
- Constant pressure PIMD
- FCP (Fictitious Cell Parameter) optimization
- Trajectory analysis tools

## Source Organization (`src/`)

### Core PIMD
- `pimd_variables.f90` - PIMD data structures (bead positions, forces, etc.)
- `pimd_subrout.f90` - Core PIMD algorithms
- `pimd_utils.f90` - Utility routines

### Engine
- `trpmd_base.f90` - Thermostatted RPMD base routines
- `trpmd_io_routines.f90` - RPMD I/O

### Input
- `ring_input_parameters_module.f90` - Input parameter definitions
- `ring_formats.f90` - Input format handling

### FCP
- `fcp_variables.f90` - Fixed cell parameter variables
- `fcp_opt_routines.f90` - FCP optimization routines

### SCF Interface
- `compute_scf_pioud.f90` - SCF calculation for each bead

## Dependencies
- PW, Modules

## Build
- `make pioud` from root
- CMake target: `pioud`
