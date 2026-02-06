# PP - Post-Processing

## Purpose
Collection of post-processing tools for analyzing and visualizing results from PWscf calculations. Covers electronic structure analysis, charge density manipulation, spectral properties, and data format conversion.

## Language
Fortran 90.

## Key Executables (in `src/`)

### Electronic Structure
- `dos.x` - Density of states
- `projwfc.x` - Projected density of states (PDOS) and projected band structure
- `bands.x` - Band structure extraction and plotting
- `fermi_velocity.x` - Fermi velocity calculation
- `fermi_proj.x` - Fermi surface projections

### Charge Density and Potentials
- `pp.x` - General post-processing (charge density, potentials, wavefunctions to various formats)
- `average.x` - Planar/macroscopic average of potentials
- `chdens.x` - Charge density analysis

### Visualization and Export
- `plotrho.x` - 2D charge density plots
- `plotband.x` - Band structure plotting utility
- `open_grid.x` - Unfold k-points from reduced to full grid

### Wannier Functions
- `pw2wannier90.x` - Interface to Wannier90 for MLWF generation

### Other Tools
- `pw2bgw.x` - Interface to BerkeleyGW
- `pw2critic.x` - Interface to CRITIC2 (Bader analysis)
- `pw2gw.x` - GW interface
- `epsilon.x` - Dielectric function (RPA)
- `molecularpdos.x` - Molecular projected DOS
- `sumpdos.x` - Sum partial DOS files
- `initial_state.x` / `final_state.x` - Core-level shift calculations
- `wfck2r.x` - Wavefunction k-space to real-space conversion

### Simple Transport
- `simple_transport/` - Boltzmann transport (Fermi integrals)

## Input Format
Most tools use namelist-based input. See `INPUT_PP.txt` for `pp.x`, etc.

## Examples
Si band structure, Ni DOS, NiO, water molecule dipole, Wannier functions, work function, molecular DOS, core-level shifts, Fermi surface.

## Dependencies
- PW (pwlibs), Modules, upflib

## Build
- `make pp` from root
- CMake target: `pp`
