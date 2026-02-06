# external - External Third-Party Libraries

## Purpose
Manages external dependencies as git submodules or vendored code. These libraries are optionally built as part of QE or can be provided externally.

## Components

### Wannier90
Maximally localized Wannier function code. Used by EPW, KCW, and PP/pw2wannier90.
- Submodule in `wannier90/`
- Can use internal (submodule) or external installation (`WANNIER90_ROOT`)

### LibMBD (Many-Body Dispersion)
Many-body dispersion correction library.
- Configured via `mbd.cmake`
- Can use internal or external (`MBD_ROOT`)

### DeviceXlib
Low-level GPU device management utilities.
- Configured via `devxlib.cmake`
- Can use internal or external (`DEVICEXLIB_ROOT`)

### FoX (Fortran XML)
Fortran XML I/O library for reading/writing XML data files.
- Configured via `fox.cmake`
- Optional, controlled by `QE_ENABLE_FOX`

### Environ
Continuum solvation library for embedded calculations.
- Configured via `environ.cmake`
- Controlled by `QE_ENABLE_ENVIRON` (NO/INTERNAL/EXTERNAL)

### LAPACK (reference)
Reference LAPACK implementation as fallback.
- Used when `QE_LAPACK_INTERNAL=ON`

## Submodule Management
- `submodule_commit_hash_records` - Tracks pinned commit hashes
- Initialize submodules: `git submodule update --init external/<name>`

## Build
- CMakeLists.txt handles conditional inclusion based on build options
- Each library has its own cmake configuration file
