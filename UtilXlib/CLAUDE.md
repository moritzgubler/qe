# UtilXlib - Utility Library

## Purpose
Lowest-level library providing MPI abstraction, error handling, timing, and memory management. All other QE libraries and packages depend on this.

## Language
Fortran 90, C (for system-level operations like timing and memory stats).

## Key Components

### MPI Wrappers (`mp_*.f90`)
Abstraction layer over MPI operations. All QE MPI communication goes through these wrappers:
- `mp.f90` - Core MP interface module
- `mp_base.f90` - Base MPI operations (send, recv, bcast, reduce, etc.)
- `mp_bands_util.f90` - Band-level parallelism utilities

### Error Handling
- `error_handler.f90` - Global error handler with MPI-safe abort

### Timing
- `clocks_handler.f90` - Fortran timing infrastructure
- `cptimer.c` - C-level high-resolution timer

### Memory
- `mem_counter.f90` - Memory allocation tracking
- `memstat.c` - System memory statistics (C interface)

### Device/GPU Support
- `device_helper.f90` - GPU device management utilities
- `nvtx_wrapper.f90` - NVIDIA NVTX profiler integration

### Other Utilities
- `divide.f90` - Work distribution utilities
- `export_gstart_2_solvers.f90` - Interface for KS_Solvers

## Preprocessor Flags
- `__MPI` - Enable MPI support
- `__CUDA` - Enable CUDA device helpers
- `__GPU_MPI` - GPU-aware MPI operations
- `__TRACE` - Execution tracing output

## Dependencies
- MPI library (optional, controlled by `__MPI`)
- CUDA toolkit (optional, for GPU features)

## Build
- `CMakeLists.txt` and `Makefile` present
- Produces library `libutil.a` (Make) or CMake target `qe_utilxlib`
- Unit tests in `tests/` subdirectory
