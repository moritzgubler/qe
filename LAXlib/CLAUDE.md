# LAXlib - Linear Algebra Library

## Purpose
Parallel distributed dense-matrix diagonalization and linear algebra operations. Supports ELPA, ScaLAPACK, and a custom parallel algorithm as backends.

## Language
Fortran 90.

## Key Components

### Core Module
- `la_module.f90` - Main public interface module
- `la_types.f90` - Data types for distributed matrices
- `la_param.f90` - Parameters and constants

### Diagonalization
- `cdiaghg.f90` - Complex Hermitian generalized eigenvalue problem
- `rdiaghg.f90` - Real symmetric generalized eigenvalue problem
- These routines dispatch to ELPA, ScaLAPACK, or custom parallel algorithm based on build configuration

### Parallel Tools
- `ptoolkit.f90` - Parallel linear algebra toolkit (custom algorithm)
- `distools.f90` - Distributed matrix tools
- `la_helper.f90` - Helper routines

### Headers
- `laxlib.h` - C header for external interface
- `laxlib_kinds.h` - Kind parameter definitions
- `laxlib_param.h` - Parameter header
- `laxlib_mid.h`, `laxlib_low.h`, `laxlib_hi.h` - Interface levels

## Testing
```bash
# After building:
cd LAXlib && make TEST
mpirun -np 4 ./la_test.x
```
- `test.f90` - Test program with configurable matrix sizes
- `tests/` - Unit test directory

## Dependencies
- UtilXlib (MPI wrappers)
- LAPACK (always required)
- ScaLAPACK (optional, `__SCALAPACK`)
- ELPA (optional, `__ELPA`)
- DeviceXlib (optional, for GPU)

## Build
- Produces `libqela.a` or CMake target `qe_laxlib`
