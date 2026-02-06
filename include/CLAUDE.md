# include - Global Header Files

## Purpose
Global include files shared across the entire QE codebase.

## Files

### `qe_version.h`
Defines the QE version number. Included by various packages to report version info.

### `cpv_device_macros.h`
GPU/device-related preprocessor macros for the CPV (Car-Parrinello) code. Defines macros for GPU memory management and kernel launching.

### `defs.h.README`
Documentation file explaining the preprocessor definitions used throughout QE (the `__FLAG` macros). Reference for understanding what each `#ifdef` controls.

## Usage
Included automatically via `-I$(TOPDIR)/include` or `include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")` in CMake.
