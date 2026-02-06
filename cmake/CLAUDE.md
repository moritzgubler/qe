# cmake - CMake Build System Modules

## Purpose
CMake configuration modules, compiler-specific flag files, find-package scripts, and helper macros for the QE build system.

## Key Files

### Helper Macros
- `qeHelpers.cmake` - Core QE CMake helper functions:
  - `qe_add_global_compile_definitions()` - Add preprocessor flags globally
  - `qe_install_targets()` - Standard target installation
  - `qe_preprocess_source()` - Fortran source preprocessing
  - `qe_enable_cuda_fortran()` - Enable CUDA Fortran for targets

### Compiler-Specific Flags
- `GNUFortranCompiler.cmake` - GCC/gfortran flags
- `IntelFortranCompiler.cmake` - Intel ifort/ifx flags
- `NVFortranCompiler.cmake` - NVHPC nvfortran flags
- `CrayFortranCompiler.cmake` - Cray ftn flags
- `IBMFortranCompiler.cmake` - IBM XL flags

### Find-Package Scripts
- `FindSCALAPACK.cmake` - Find ScaLAPACK
- `FindELPA.cmake` - Find ELPA library
- `FindLibxc.cmake` - Find Libxc

### Other
- `unit_test.cmake` - Unit test framework setup
- `GitInfo.cmake` - Git version/branch info extraction
- `AddRPATH.cmake` - RPATH management for installed binaries
- `CMakeGraphVizOptions.cmake` - Dependency graph visualization
- `quantum_espresso.pc.in` - pkg-config template
- `qeConfig.cmake.in` - CMake package config template

## Usage
Included automatically by the root `CMakeLists.txt` via `set(CMAKE_MODULE_PATH ...)`.
