# PW - Plane-Wave Self-Consistent Field (PWscf)

## Purpose
The core DFT engine of Quantum ESPRESSO. Performs self-consistent electronic structure calculations using plane waves and pseudopotentials. Supports structural optimization, molecular dynamics, band structure calculations, and many advanced features.

## Language
Fortran 90 (~250+ source files in `src/`).

## Main Executable
`pw.x` - The primary QE executable for ground-state DFT calculations.

## Calculation Types (`calculation` input keyword)
- `'scf'` - Self-consistent field
- `'nscf'` - Non-self-consistent (band structure, DOS)
- `'bands'` - Band structure along high-symmetry paths
- `'relax'` - Ionic relaxation (BFGS, damped dynamics)
- `'vc-relax'` - Variable-cell relaxation
- `'md'` - Born-Oppenheimer molecular dynamics
- `'vc-md'` - Variable-cell MD

## Source Organization (`src/`)

### SCF Engine
- `electrons.f90` - Main SCF loop
- `c_bands.f90` - Band structure calculation (calls KS_Solvers)
- `scf_mod.f90` - SCF data structures (charge density, potentials)
- `mix_rho.f90` - Charge density mixing (Broyden, etc.)
- `v_of_rho.f90` - Potential from density (Hartree + XC + external)

### Hamiltonian Application
- `h_psi.f90` - H|psi> (Hamiltonian acting on wavefunction)
- `s_psi.f90` - S|psi> (overlap operator for ultrasoft/PAW)
- `vloc_psi.f90` - Local potential contribution
- `add_vuspsi.f90` - Ultrasoft/PAW nonlocal contribution

### Forces and Stress
- `forces.f90` - Total force calculation
- `stress.f90` - Total stress tensor
- `force_us.f90` - Ultrasoft PP force contributions
- `force_cc.f90` - Core-correction forces

### Wavefunctions
- `wfcinit.f90` - Wavefunction initialization
- `rotate_wfc.f90` - Subspace rotation
- `orthoatwfc.f90` - Atomic wavefunction orthogonalization

### Special Methods
- `exx.f90` - Exact exchange (EXX / hybrid functionals)
- `ldaU.f90` - DFT+U (Hubbard correction)
- `esm.f90` - Effective Screening Medium (slab/wire boundary conditions)
- `add_efield.f90` - External electric field (sawtooth potential)
- `add_bfield.f90` - External magnetic field
- `bp_c_phase.f90` - Berry phase calculations

### I/O
- `punch.f90` - Write data files
- `pw_restart_new.f90` - Restart file I/O (XML format)
- `read_file_new.f90` - Read saved data

### Setup and Initialization
- `input.f90` - Input parsing
- `setup.f90` - Calculation setup
- `init_run.f90` - Run initialization

## Additional Tools (`tools/`)
- `ibrav2cell.x` - Convert ibrav to explicit cell vectors
- `cell2ibrav.x` - Convert cell vectors to ibrav
- `ev.x` - Equation of state fitting (energy vs. volume)
- `kpoints.x` - k-point generation utility
- `pwi2xsf.x` - Convert PW input to XSF visualization format
- `scan_ibrav.x` - Scan ibrav parameters

## Input Format
Fortran namelist-based input (see `INPUT_PW.txt` or `Doc/INPUT_PW.html`):
- `&CONTROL` - Job type, file paths, convergence
- `&SYSTEM` - Crystal structure, cutoffs, XC functional, smearing
- `&ELECTRONS` - SCF convergence parameters, mixing
- `&IONS` - Relaxation/MD parameters
- `&CELL` - Variable-cell parameters
- Card sections: `ATOMIC_SPECIES`, `ATOMIC_POSITIONS`, `K_POINTS`, `CELL_PARAMETERS`

## Examples
Located in `examples/`: SCF (Si, Al), band structure, MD, EXX, magnetism, VCS relaxation, clusters, electric fields, DFT+U, etc.

## Dependencies
- Modules, FFTXlib, LAXlib, KS_Solvers, upflib, XClib, dft-d3
- External: Wannier90 (optional), LibMBD (optional)

## Build
- `make pw` from root, or `make all` in `PW/`
- Produces `pw.x` and utility executables in `PW/src/` and `PW/tools/`
- CMake target: `pw`
