# atomic - Pseudopotential Generation (LD1)

## Purpose
Atomic code for all-electron calculations and pseudopotential generation. Generates norm-conserving, ultrasoft, and PAW pseudopotentials in UPF format.

## Language
Fortran 90.

## Main Executable
`ld1.x`

## Key Features
- All-electron atomic calculations (any atom in periodic table)
- Norm-conserving pseudopotential generation
- Ultrasoft pseudopotential generation
- PAW dataset generation
- Pseudopotential testing (logarithmic derivatives, transferability)
- Multiple exchange-correlation functionals
- Scalar-relativistic and fully-relativistic calculations
- Semi-core state handling

## Source Organization (`src/`)
- All-electron solver and radial equation integration
- Pseudopotential construction routines
- Transferability testing
- UPF output generation

## Workflow
1. Configure atomic parameters (element, XC, cutoff radii)
2. Run `ld1.x` to generate pseudopotential
3. Test transferability with logarithmic derivative comparison
4. Use generated `.UPF` file in PW/CP calculations

## Examples
`examples/paw_examples/` - PAW pseudopotential generation examples.

## Dependencies
- Modules (basic), UtilXlib

## Build
- `make ld1` from root
- CMake target: `ld1`
