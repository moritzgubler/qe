# pseudo - Sample Pseudopotential Files

## Purpose
Repository of sample pseudopotential files in UPF (Unified Pseudopotential Format) used by examples and tests throughout QE.

## Format
UPF v2 (XML-based) and UPF v1 (plain text). Files have `.UPF` or `.upf` extension.

## Contents
Pseudopotentials for common elements used in QE's examples and test suite (Si, C, Al, Cu, Ni, O, H, etc.). These are for testing only - production calculations should use validated pseudopotential libraries.

## Pseudopotential Sources
Production pseudopotentials can be downloaded from:
- SSSP (Standard Solid State Pseudopotentials): https://www.materialscloud.org/discover/sssp
- PseudoDojo: http://www.pseudo-dojo.org/
- QE pseudopotential page: https://pseudopotentials.quantum-espresso.org/

## Cleanup
`clean_ps` script removes downloaded pseudopotentials, keeping only those in the repository.
