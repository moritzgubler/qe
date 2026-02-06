# CPV - Car-Parrinello Molecular Dynamics

## Purpose
Car-Parrinello molecular dynamics (CPMD) code. Performs MD simulations using fictitious electron dynamics instead of explicit SCF at each step, enabling efficient ab initio MD for large systems and long timescales.

## Language
Fortran 90.

## Main Executables
- `cp.x` - Main CP molecular dynamics program
- `manycp.x` - Multiple simultaneous CP runs
- `cppp.x` - CP post-processing
- `wfdd.x` - Wannier function dynamics/analysis

## Key Features
- Fictitious electron dynamics (Car-Parrinello scheme)
- Wavefunction optimization (conjugate gradient, damped dynamics)
- Multiple thermostat options (Nose-Hoover, etc.)
- Variable-cell dynamics
- EXX (exact exchange) support
- PAW and ultrasoft pseudopotentials
- Hubbard U corrections
- Wannier function calculations
- Autopilot module for complex simulation protocols

## Source Organization (`src/`)

### Core MD Engine
- `cpr_mod.f90` / `cpr.f90` - Main CP module and driver
- `cpr_loop.f90` - CP time-step loop
- `cp_emass.f90` - Fictitious electron mass management
- `cp_wavefunctions.f90` - Wavefunction data structures

### Dynamics and Forces
- `forces.f90` - Force computation
- `stress.f90` - Stress tensor
- `move_electrons.f90` - Electron dynamics integration
- `move_ions.f90` - Ion dynamics integration

### Post-Processing
- `cppp.f90` - Post-processing main program
- `wfdd.f90` - Wannier function computation

### Autopilot
- Uses `autopilot.f90` from Modules to allow runtime parameter changes during long MD runs

## Input Format
Fortran namelist-based (see `INPUT_CP.txt`):
- `&CONTROL` - Job type (cp, vc-cp, etc.)
- `&SYSTEM` - System parameters
- `&ELECTRONS` - Electron dynamics parameters (emass, dt, etc.)
- `&IONS` - Ion dynamics parameters
- `&CELL` - Cell dynamics parameters

## Documentation
- `Doc/user_guide.md` - CP user guide
- `Doc/autopilot_guide.md` - Autopilot documentation

## Examples
9 example directories covering basic MD, Wannier functions, variable cell, EXX, PAW, restart, external fields, etc.

## Dependencies
- Modules, FFTXlib, LAXlib, upflib, XClib

## Build
- `make cp` from root
- CMake target: `cp`
