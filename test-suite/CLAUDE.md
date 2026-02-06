# test-suite - Integration and Regression Tests

## Purpose
Comprehensive test suite for validating all QE packages. Contains input files, reference outputs, and comparison scripts for automated testing.

## Running Tests

### With CMake
```bash
cd build && ctest
```

### With Make
```bash
cd test-suite
make run-tests          # Run all tests
make run-tests-pw       # PW tests only
make run-tests-ph       # PHonon tests only
make run-tests-cp       # CP tests only
make run-tests-epw      # EPW tests only
# etc.
```

## Test Naming Convention
Tests are organized by package prefix:
- `pw_*` - PWscf tests (SCF, NSCF, MD, relax, DFT+U, EXX, etc.)
- `ph_*` - PHonon DFPT tests
- `cp_*` - Car-Parrinello MD tests
- `epw_*` - Electron-phonon tests
- `hp_*` - Hubbard parameter tests
- `neb_*` - NEB pathway tests
- `tddfpt_*` - TD-DFPT spectroscopy tests
- `xspectra_*` - X-ray spectra tests
- `gw_*` - GW calculation tests
- `kcw_*` - Koopmans functional tests
- `pp_*` - Post-processing tests

## Test Structure
Each test directory contains:
- Input files (`.in`, `.pwi`, `.phi`)
- Reference output files for comparison
- Test-specific configuration

## Comparison Framework
- `testcode/` - External test comparison framework (testcode)
- Numerical tolerance-based comparison (not exact text diff)
- Handles platform-dependent variations in FFT grids, convergence, etc.

## Notes
- Tests are designed for serial or small parallel (2-4 cores) execution
- Do not use for benchmarking parallel performance
- Some tests require specific optional features (HDF5, Libxc, etc.)
