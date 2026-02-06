# upflib - Pseudopotential Library

## Purpose
Handles all pseudopotential-related operations: reading/writing UPF (Unified Pseudopotential Format) files, format conversion, projector/beta function setup, PAW support, and matrix element generation.

## Language
Fortran 90.

## Key Components

### UPF Format I/O
- `read_upf_new.f90` - Modern UPF reader (v2 XML format)
- `read_upf_v1.f90` - Legacy UPF v1 reader
- `write_upf_new.f90` - UPF writer
- `upf_to_internal.f90` - Convert UPF data to internal representation
- `pseudo_types.f90` - Core data types for pseudopotential storage

### Format Conversion (other formats -> UPF)
- `read_ps.f90` - Generic pseudopotential reader (dispatches to specific formats)
- `read_cpmd.f90`, `cpmd2upf.f90` - CPMD format
- `read_fhi.f90`, `fhi2upf.f90` - FHI format
- `read_ncpp.f90`, `ncpp2upf.f90` - NCPP format
- `read_uspp.f90`, `uspp2upf.f90` - Vanderbilt USPP format
- `casino2upf.f90` - CASINO format
- `hgh2qe.f90` - HGH/GTH Goedecker format
- `read_psml.f90` - PSML format

### Beta Functions and Projectors
- `beta_mod.f90` - Beta function data structures
- `gen_us_dj.f90`, `gen_us_dy.f90` - Generate US pseudopotential projectors
- `init_us_0.f90`, `init_us_b0.f90` - Initialize US PP setup
- `qvan2.f90`, `dqvan2.f90` - Q (augmentation) function computation

### Radial Grid and Spherical Harmonics
- `radial_grids.f90` - Radial grid management
- `sph_bes.f90` - Spherical Bessel functions
- `ylmr2.f90` - Real spherical harmonics
- `qrad_mod.f90` - Radial integration for augmentation

### PAW Support
- PAW augmentation and one-center terms handled through various `init_*` and `gen_*` routines

### Executables
- `upfconv.x` - Pseudopotential format converter
- `virtual_v2.x` - Virtual crystal approximation PP generator
- `casino2upf.x` - CASINO to UPF converter

## Dependencies
- UtilXlib (MPI, error handling)
- LAPACK (for matrix operations)
- DeviceXlib (optional, GPU)

## Build
- Produces `libqeupf.a` or CMake target `qe_upflib`
