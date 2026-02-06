# XClib - Exchange-Correlation Library

## Purpose
Implements exchange-correlation (XC) functionals for DFT calculations. Supports LDA, GGA, meta-GGA, and hybrid functionals, with optional Libxc integration.

## Language
Fortran 90, C (for BEEF functional interface).

## Architecture

### Drivers (top-level dispatch)
- `qe_drivers_lda_lsda.f90` - LDA/LSDA functional drivers
- `qe_drivers_gga.f90` - GGA functional drivers
- `qe_drivers_mgga.f90` - Meta-GGA functional drivers
- `qe_drivers_d_lda_lsda.f90`, `qe_drivers_d_gga.f90`, `qe_drivers_d_mgga.f90` - Derivative drivers

### Functional Implementations
- `qe_funct_corr_lda_lsda.f90` - LDA/LSDA correlation (PZ, VWN, PW, etc.)
- `qe_funct_exch_lda_lsda.f90` - LDA/LSDA exchange (Slater, etc.)
- `qe_funct_corr_gga.f90` - GGA correlation (PBE, LYP, etc.)
- `qe_funct_exch_gga.f90` - GGA exchange (PBE, B88, etc.)
- `qe_funct_mgga.f90` - Meta-GGA (SCAN, TPSS, rSCAN, etc.)

### Libxc Wrappers
- `xc_wrapper_lda_lsda.f90`, `xc_wrapper_gga.f90`, `xc_wrapper_mgga.f90` - Interface to external Libxc library when `__LIBXC` is defined

### Configuration
- `dft_setting_routines.f90` - DFT functional selection and configuration
- `dft_setting_params.f90` - Functional parameter definitions
- `xclib_utils_and_para.f90` - Utility and parallelization helpers

### BEEF-vdW
- `xc_beef_interface.f90` - Fortran interface
- `beef_interface.c`, `beefun.c` - C implementation of BEEF functional

### Testing
- `xclib_test.f90` - Standalone XC library test program
- `xclib_test_defs.py` - Test definition generator
- `xclib_test_references/` - Reference data for validation

## Dependencies
- Libxc >= 5.1.2 (optional, `__LIBXC` flag)
- No internal QE dependencies (standalone library)

## Build
- Produces `libqexclib.a` or CMake target `qe_xclib`
