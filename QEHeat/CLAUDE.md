# QEHeat - Energy Current for Thermal Transport

## Purpose
Computes energy current (heat flux) in insulators for thermal transport calculations within DFT. Implements the approach for computing thermal conductivity via Green-Kubo relations from ab initio molecular dynamics trajectories.

## Language
Fortran 90.

## Main Executable
`all_currents.x`

## Key Features
- Energy current from wavefunctions and forces
- Thermal transport properties via Green-Kubo
- Support for Car-Parrinello (CP) trajectories
- Heat flux decomposition (kinetic, potential, convective)
- Single-point and trajectory analysis modes

## Source Organization (`src/`)
- `all_currents.f90` - Main driver
- `kohn_sham_mod.f90` - KS energy current contribution
- `averages.f90` - Trajectory averaging
- `scf_result.f90` - SCF result handling
- `cpv_traj.f90` - CP trajectory reader
- `traj_object.f90` - Trajectory data structure
- `ec_functionals.f90` - Energy current functional evaluation
- `hartree_xc_mod.f90` - Hartree and XC contributions to current

## Workflow
1. Run CP or PW molecular dynamics to generate trajectory
2. Run `all_currents.x` to compute energy current at each snapshot
3. Post-process to get thermal conductivity via Green-Kubo

## Examples
H2O trajectory, SiO2 single point, small H2O with CP.

## Dependencies
- PW, CPV, Modules

## Build
- `make all_currents` from root
- CMake target: `all_currents`
