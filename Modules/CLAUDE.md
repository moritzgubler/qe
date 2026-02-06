# Modules - Core Data Structures and Utilities

## Purpose
Central module providing the fundamental data structures, constants, I/O infrastructure, and parallelization framework used by all QE application packages. This is the main "glue" library.

## Language
Fortran 90.

## Key Data Structures

### Crystal and Cell
- `cell_base.f90` - Lattice vectors, volume, reciprocal lattice, cell dynamics
- `ions_base.f90` - Atomic positions, velocities, masses, species, constraints
- `symm_base.f90` - Crystal symmetry operations, space group detection

### Electronic Structure
- `electrons_base.f90` - Electron count, occupations, spin
- `gvecw.f90` - Wavefunction G-vector set
- `start_k.f90` / `kpoint_grid.f90` - k-point generation and management
- `control_flags.f90` - Calculation control parameters (scf, relax, md, etc.)

### FFT and Grids
- `fft_base.f90` - FFT descriptor initialization and management for QE
- `fft_rho.f90` - Charge density FFT operations
- `fft_wave.f90` - Wavefunction FFT operations
- `recvec.f90` - Reciprocal space vector management

### I/O
- `io_files.f90` - File path and unit management
- `io_base.f90` - Base I/O routines
- `io_global.f90` - Global I/O (stdout, ionode management)
- `read_input.f90` - Master input reader
- `input_parameters.f90` - All input namelist parameters (central input definition)

### Parallelism
- `mp_global.f90` - Global MPI communicator setup
- `mp_pools.f90` - k-point pool communicators
- `mp_bands.f90` - Band group communicators
- `mp_images.f90` - Image parallelism communicators
- `mp_exx.f90` - Exact exchange parallelism

### Constants and Kinds
- `constants.f90` - Physical constants in atomic units
- `kind.f90` - Fortran kind parameters (DP, etc.)

### Dispersion Corrections
- `mm_dispersion.f90` - Grimme DFT-D2 dispersion
- `lj_correction.f90`, `london_module.f90` - Lennard-Jones / London corrections

### RISM (Solvent Model)
- `1drism.f90`, `3drism.f90`, `lauerism.f90` - Reference Interaction Site Model implementation

### Other Key Modules
- `noncol.f90` - Non-collinear magnetism data
- `extfield.f90` - External fields (electric, gate, etc.)
- `constraints_module.f90` - Geometric constraints for MD/relaxation
- `autopilot.f90` - On-the-fly parameter modification during MD
- `bz_form.f90` - Brillouin zone utilities
- `latgen.f90` - Lattice generation from ibrav parameters
- `compute_dipole.f90` - Dipole moment calculations
- `pw_dot.f90` - Plane-wave dot products

## Dependencies
- UtilXlib, FFTXlib, LAXlib, upflib, XClib

## Build
- Produces `libqemod.a` or CMake target `qe_modules`
