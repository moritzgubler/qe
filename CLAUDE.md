# Quantum ESPRESSO - Project Guide

## Overview

Quantum ESPRESSO (QE) is an open-source suite for electronic-structure calculations and materials modeling at the nanoscale, based on density-functional theory (DFT), plane waves, and pseudopotentials. Version 7.5, licensed under GPLv2.

- **Languages:** Primarily Fortran 90 with some C. Build tooling in CMake, Make, shell, Python, Tcl/Perl.
- **Repository:** https://gitlab.com/QEF/q-e
- **Documentation:** https://www.quantum-espresso.org/Doc/user_guide/
- **Wiki:** https://gitlab.com/QEF/q-e/-/wikis/home

## Build System

Two build systems are supported:

### CMake (recommended for new work)
```bash
mkdir build && cd build
cmake -DCMAKE_Fortran_COMPILER=mpif90 -DCMAKE_C_COMPILER=mpicc ..
make -jN
```

### Make (legacy)
```bash
./configure [options]
make [-jN] target    # targets: pw, ph, cp, pp, hp, epw, all, etc.
```

### Key CMake Options
| Option | Default | Description |
|--------|---------|-------------|
| `QE_ENABLE_MPI` | ON | MPI parallelization |
| `QE_ENABLE_OPENMP` | OFF | OpenMP threading |
| `QE_ENABLE_CUDA` | OFF | GPU acceleration (requires NVHPC compiler) |
| `QE_ENABLE_OPENACC` | OFF | OpenACC (auto-ON with CUDA) |
| `QE_ENABLE_SCALAPACK` | OFF | Distributed linear algebra |
| `QE_ENABLE_ELPA` | OFF | ELPA eigensolver (requires SCALAPACK) |
| `QE_ENABLE_LIBXC` | OFF | External XC functionals (>=5.1.2) |
| `QE_ENABLE_HDF5` | OFF | HDF5 I/O support |
| `QE_ENABLE_TEST` | ON | Unit and system tests |
| `QE_LAPACK_INTERNAL` | OFF | Use bundled reference LAPACK |
| `QE_FFTW_VENDOR` | AUTO | FFT backend: Intel_DFTI, FFTW3, ArmPL, etc. |

### Preprocessor Macros
Feature flags used throughout the codebase (set by build system):
- `__MPI` - MPI enabled
- `__CUDA` - CUDA GPU acceleration
- `__OPENMP_GPU` - OpenMP offload
- `__SCALAPACK` - ScaLAPACK enabled
- `__ELPA` / `__ELPA_2016` / `__ELPA_2015` - ELPA version
- `__HDF5` / `__HDF5_SERIAL` - HDF5 support
- `__LIBXC` - Libxc external XC library
- `__TRACE` - Execution tracing
- `__GPU_MPI` - GPU-aware MPI
- `__ENVIRON` - Environ solvation library
- `__OSCDFT` - OS-CDFT support
- `__MPI_MODULE` - Use Fortran `mpi` module vs `mpif.h`

## Architecture and Directory Structure

### Core Libraries (built first, used by all packages)
| Directory | Purpose |
|-----------|---------|
| `UtilXlib/` | MPI wrappers, error handling, timing, memory management |
| `FFTXlib/` | Parallel 3D FFTs with MPI/OpenMP, multiple backends |
| `LAXlib/` | Parallel dense-matrix diagonalization (ELPA, ScaLAPACK, custom) |
| `XClib/` | Exchange-correlation functionals (LDA, GGA, meta-GGA, hybrids) |
| `upflib/` | Pseudopotential I/O, UPF format handling, PAW support |
| `Modules/` | Central data structures: cell, ions, wavefunctions, k-points, symmetry, I/O, parallelization |
| `KS_Solvers/` | Kohn-Sham eigensolvers: Davidson, CG, RMM, ParO, PPCG |

### Main Application Packages
| Directory | Package | Executable(s) | Purpose |
|-----------|---------|---------------|---------|
| `PW/` | PWscf | `pw.x` | Core DFT SCF, structural optimization, MD |
| `CPV/` | CP | `cp.x`, `cppp.x` | Car-Parrinello molecular dynamics |
| `PHonon/` | PHonon | `ph.x` | Phonons and dielectric properties via DFPT |
| `PP/` | PostProc | `pp.x`, `projwfc.x`, `dos.x`, `bands.x` | Post-processing and analysis |
| `EPW/` | EPW | `epw.x` | Electron-phonon coupling with Wannier functions |
| `TDDFPT/` | turboTDDFT | `turbo*.x` | Time-dependent DFPT spectroscopy |
| `HP/` | HP | `hp.x` | Hubbard parameters from DFPT |
| `NEB/` | PWneb | `neb.x` | Nudged Elastic Band transition states |
| `PWCOND/` | PWCOND | `pwcond.x` | Ballistic transport |
| `XSpectra/` | XSpectra | `xspectra.x` | X-ray absorption spectra |
| `GWW/` | GWL | `gww.x`, `bse.x` | GW many-body perturbation theory |
| `QEHeat/` | QEHeat | `all_currents.x` | Energy current for thermal transport |
| `KCW/` | KCW | `kcw.x` | Koopmans-compliant functionals |
| `PIOUD/` | PIOUD | various | Path integral MD |
| `atomic/` | LD1 | `ld1.x` | Pseudopotential generation |
| `COUPLE/` | COUPLE | library | Interface for coupling to external codes |

### Support Directories
| Directory | Purpose |
|-----------|---------|
| `LR_Modules/` | Shared linear-response DFPT algorithms |
| `dft-d3/` | Grimme DFT-D3 dispersion corrections |
| `external/` | Git submodules: Wannier90, LibMBD, FoX, DeviceXlib |
| `include/` | Global headers: `qe_version.h`, GPU macros |
| `cmake/` | CMake modules, compiler flags, find-package scripts |
| `install/` | Configure support, external library build scripts |
| `pseudo/` | Sample pseudopotential files (UPF format) |
| `test-suite/` | Comprehensive integration tests for all packages |
| `dev-tools/` | Developer utilities: helpdoc, memory analysis, code normalization |
| `Doc/` | General documentation: user guide, Brillouin zones, changelogs |
| `GUI/` | PWgui (Tcl/Tk) graphical interface and Emacs modes |
| `build/` | Out-of-source CMake build directory |
| `archive/` | Archived/legacy files |

## Build Dependency Order

```
UtilXlib -> FFTXlib
UtilXlib -> LAXlib -> KS_Solvers
UtilXlib -> upflib
           XClib
UtilXlib + FFTXlib + LAXlib + upflib + XClib -> Modules
Modules + KS_Solvers + dft-d3 -> PW
PW -> {PP, NEB, PWCOND, XSpectra, PIOUD, COUPLE}
Modules -> LR_Modules
PW + LR_Modules -> PHonon -> {EPW, GWW}
PW + LR_Modules -> {TDDFPT, HP, KCW}
PW + CPV -> QEHeat
Modules -> CPV
Modules -> atomic
```

## Coding Conventions

- **Fortran standard:** Fortran 90/95 free-format (`.f90` extension). Some legacy fixed-format files exist.
- **Preprocessing:** C preprocessor directives (`#if`, `#ifdef`) used extensively for feature flags. Files needing preprocessing typically have `.f90` extension and are preprocessed by the build system.
- **Module naming:** Fortran modules generally match their filename (e.g., `cell_base.f90` contains `MODULE cell_base`).
- **Parallelism model:** MPI for distributed memory (primary), OpenMP for shared memory (secondary), OpenACC/CUDA Fortran for GPU.
- **Data distribution:** Plane waves, G-vectors, and real-space grids distributed across MPI ranks. Managed by FFTXlib and Modules.
- **Units:** Atomic (Rydberg) units internally: energies in Ry, lengths in Bohr, etc.
- **I/O format:** XML for restart/data files, Namelist-based Fortran input files, UPF for pseudopotentials.

## Parallelism Levels

QE supports multiple levels of MPI parallelism (flags for `mpirun ... executable -flags`):
- `-ni` / `-nimage`: Image parallelism (NEB, PHonon)
- `-nk` / `-npools`: k-point pools
- `-nb` / `-nbgrp`: Band groups
- `-nt` / `-ntg`: Task groups (FFT parallelism)
- `-nd` / `-ndiag`: Linear algebra processors

## Testing

```bash
# CMake
cd build && ctest

# Make
cd test-suite && make run-tests

# Specific package
cd test-suite && make run-tests-pw
```

Test naming convention: `pw_*`, `ph_*`, `cp_*`, `epw_*`, `hp_*`, `tddfpt_*`, `xspectra_*`, etc.

## Key Environment Variables

- `BIN_DIR` - Directory with QE executables
- `PSEUDO_DIR` - Pseudopotential file directory
- `TMP_DIR` - Scratch/temporary directory for calculations
- `PARA_PREFIX` - MPI launcher prefix (e.g., `mpirun -np 4`)
- `PARA_POSTFIX` - Parallelization flags (e.g., `-nk 1 -nd 1`)
- `OMP_NUM_THREADS` - OpenMP thread count

## GPU Support

Requires NVHPC (nvfortran) compiler. Uses CUDA Fortran + OpenACC.
```bash
./configure --with-cuda=$CUDA_HOME --with-cuda-cc=70 --with-cuda-runtime=11.0 --enable-openmp
# or with CMake:
cmake -DQE_ENABLE_CUDA=ON -DCMAKE_Fortran_COMPILER=nvfortran ...
```
