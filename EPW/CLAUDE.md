# EPW - Electron-Phonon Wannier

## Purpose
Calculates electron-phonon coupling properties using Wannier function interpolation. Enables dense Brillouin zone sampling for transport, superconductivity, and optical properties from coarse DFPT grids.

## Language
Fortran 90. Version: EPW v5.9.

## Main Executable
`epw.x`

## Key Features
- Electron-phonon matrix elements via Wannier interpolation
- Phonon-limited carrier transport (resistivity, mobility)
- Superconductivity (Migdal-Eliashberg, anisotropic)
- Phonon-assisted optical absorption
- Electron self-energy and spectral functions
- Phonon self-energy and linewidths
- Wannier function generation (interface to Wannier90)

## Source Organization (`src/`)

### Core Modules
- `global_var.f90` - Global variables and parameters
- `ep_constants.f90` - Physical constants
- `input.f90` - Input parsing

### Wannier Interpolation
- `wannier.f90` - Wannier function management
- `bloch2wannier.f90` - Bloch to Wannier transformation
- `wannier2bloch.f90` - Wannier to Bloch interpolation
- `pw2wan.f90` - PWscf to Wannier90 interface

### Transport and Properties
- `transport.f90` - Carrier transport calculations
- `transport_legacy.f90` - Legacy transport routines
- `supercond.f90` - Superconductivity calculations
- `selfen.f90` - Self-energy calculations

### I/O (`io/`)
- `io.f90` - General I/O
- `io_transport.f90` - Transport data I/O
- `io_selfen.f90` - Self-energy I/O

### Post-Processing
- `ZG/` directory with additional tools:
  - Joint density of states
  - k-point generation utilities
  - Special displacement method

## Workflow
1. `pw.x` SCF on coarse k-grid
2. `ph.x` DFPT on coarse q-grid
3. `epw.x` interpolates to dense grids using Wannier functions

## Input Format
Namelist-based input (see `Doc/INPUT_EPW.txt`).

## Examples
Diamond, GaN, MgB2, Pb, SiC, LiF - covering transport, superconductivity, and phonon properties.

## Dependencies
- PW, PHonon, Modules, Wannier90 (external), ELPA (optional)

## Build
- `make epw` from root (also builds Wannier90)
- CMake target: `epw`
- Executable linked to `bin/epw.x`
