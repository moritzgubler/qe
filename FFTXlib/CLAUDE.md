# FFTXlib - Parallel FFT Library

## Purpose
Parallel distributed 3D FFT library with load-balanced data distribution. Handles plane-wave/G-vector distribution across MPI ranks and supports multiple FFT backends including GPU.

## Language
Fortran 90, C (for low-level FFT wrappers and FFTW interfaces).

## Architecture

### Descriptor Types
- `fft_types.f90` - Main FFT descriptor type (`fft_type_descriptor`), defines grid dimensions, data distribution, stick decomposition
- `stick_base.f90` - Stick (1D column) decomposition for 3D FFT parallelization
- `fft_smallbox_type.f90` - Small-box FFT descriptors (for augmentation charges)
- `fft_param.f90` - Compile-time parameters (max dimensions, etc.)

### Parallel Execution
- `fft_interfaces.f90` / `fft_fwinv.f90` - High-level forward/inverse FFT interfaces
- `fft_parallel.f90` - MPI-parallel 3D FFT implementation
- `scatter_mod.f90` - MPI scatter/gather for FFT data redistribution
- `tg_gather.f90` - Task group gather operations
- `fft_interpolate.f90` - Grid interpolation between different FFT meshes
- `fft_smallbox.f90` - Small-box FFT for augmentation charges

### Backend Wrappers (`fft_scalar.*.f90`)
Each file wraps a specific FFT library:
- `fft_scalar.FFTW3.f90` - FFTW3 (most common)
- `fft_scalar.DFTI.f90` - Intel MKL DFTI
- `fft_scalar.ESSL.f90` - IBM ESSL
- `fft_scalar.SX6.f90` - NEC SX
- `fft_scalar.cuFFT.f90` - NVIDIA cuFFT (GPU)
- C files: `fft_stick.c`, `fftw.c`, `fftw_dp.c`, `fftw_sp.c` - Internal scalar FFTW

### Helper Routines
- `fft_ggen.f90` - G-vector generation and Miller index mapping
- `fft_helper_subroutines.f90` - Common helper operations
- `fft_support.f90` - Utility functions
- `fft_error.f90` - FFT-specific error handling

## Testing
```bash
# After configuring QE:
cd FFTXlib && make TEST
mpirun -np 4 ./fft_test.x -ecutwfc 80 -alat 20 -nbnd 128 -ntg 4
```
- `fft_test.f90` / `test.f90` - Test and benchmark program
- `gen_test_params.py` - Extract FFT parameters from pw.x output

## Key Concepts
- **Stick decomposition:** 3D FFT grid split into 1D "sticks" along z-axis, distributed across MPI ranks
- **Task groups:** Additional parallelism level for FFT operations (`-ntg` flag)
- **Dual-grid:** Support for separate wavefunction and charge density grids (ecutwfc vs ecutrho)

## Dependencies
- UtilXlib (MPI wrappers)
- External FFT library (FFTW3, MKL, etc.) or internal implementation
- CUDA toolkit (optional, for cuFFT)

## Build
- Produces `libqefft.a` or CMake target `qe_fftxlib`
