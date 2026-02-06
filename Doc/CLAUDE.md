# Doc - General Documentation

## Purpose
General documentation for the Quantum ESPRESSO suite. Package-specific documentation is in each package's own `Doc/` subdirectory.

## Key Documents

### User Guide
- `user_guide.tex` / `user_guide.pdf` - Main QE user guide covering installation, configuration, and parallelism

### Reference Materials
- `brillouin_zones.tex` / `brillouin_zones.pdf` - Brillouin zone labels and high-symmetry paths for all crystal systems
- `Hubbard_input.tex` - Guide for DFT+U/DFT+U+V input specification
- `ExternalForceFields.tex` - External force field documentation
- `constraints_HOWTO.tex` - Geometric constraints for MD/relaxation
- `plumed_quick_ref.tex` - PLUMED metadynamics integration

### Change Logs
- `release-notes` - Current release notes and bug fixes
- `ChangeLog.pw` - Historical PW changelog (pre-2004)
- `ChangeLog.cp` - Historical CP changelog (pre-2004)
- `ChangeLog.old` - Historical QE changelog (post-2004, no longer updated)

### Images
- PNG/PDF figures for crystal structures, Brillouin zones, and concepts

## Building Documentation
```bash
make doc    # from root (requires latex2html, pdflatex, tcl, tcllib, xsltproc)
```

## Package-Specific Documentation
Each package has its own `Doc/` with:
- `INPUT_<PROGRAM>.txt` - Plain-text input documentation
- `INPUT_<PROGRAM>.html` - HTML input documentation
- `INPUT_<PROGRAM>.def` - helpdoc definition file (source of truth)
- Package-specific user guides where applicable
