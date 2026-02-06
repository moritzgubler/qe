# NEB - Nudged Elastic Band

## Purpose
Calculates minimum energy pathways and transition states between two configurations using the Nudged Elastic Band (NEB) method. Used for studying reaction barriers, diffusion, and phase transitions.

## Language
Fortran 90.

## Main Executable
`neb.x`

## Key Features
- Nudged Elastic Band (NEB) method
- Climbing Image NEB (CI-NEB) for accurate transition states
- String method variant
- Path reparametrization
- Automatic interpolation between endpoints
- Image parallelism (`-ni` flag)

## Source Organization (`src/`)

### Path Engine
- `path_variables.f90` - Path data structures and parameters
- `path_io_tools.f90` - Path I/O utilities
- `path_reparametrisation.f90` - Path reparametrization algorithms
- `path_io_units_module.f90` - I/O unit management

### Engine Interface
- `engine_to_path_pos.f90` - Position interface between DFT engine and path
- `engine_to_path_alat.f90` - Cell parameter interface
- `engine_to_path_nat.f90` - Atom count interface
- `engine_to_path_tot_charge.f90` - Charge interface
- `engine_to_path_fix_atom_pos.f90` - Fixed atom handling

## Input Format
Namelist-based (see `INPUT_NEB.txt`):
- `&PATH` - NEB parameters (num_of_images, CI_scheme, etc.)
- `BEGIN ... END` blocks for initial and final configurations

## Workflow
1. Prepare initial and final atomic configurations
2. Run `neb.x` to find minimum energy pathway
3. Analyze energy profile and transition state geometry

## Examples
H2+H reaction barrier with various image configurations.

## Dependencies
- PW (pwlibs), Modules, UtilXlib

## Build
- `make neb` from root
- CMake target: `neb`
