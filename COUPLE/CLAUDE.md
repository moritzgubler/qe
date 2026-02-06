# COUPLE - External Code Coupling Interface

## Purpose
Library interface for coupling Quantum ESPRESSO with external codes. Provides a C-compatible API to drive PW and CP calculations from other programs.

## Language
Fortran 90 with C-compatible interfaces.

## Key Features
- Library mode interface for PW and CP
- Position/force/energy exchange with external drivers
- Usable from C, C++, Python, or other languages via C bindings

## Source Organization
- Fortran modules providing callable interfaces to PW/CP engines
- Example driver programs

## Examples
`examples/` - Demonstrates coupling interface usage.

## Dependencies
- PW, CPV, Modules

## Build
- `make couple` from root (requires `pw` and `cp` first)
- CMake target: `couple`
