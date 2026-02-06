# install - Installation Support

## Purpose
Configuration support scripts, compiler detection, external library build helpers, and plugin management for the legacy Make-based build system.

## Key Files

### Configuration
- `configure` support files for autoconf-based configuration
- Compiler detection and flag generation
- Produces `make.inc` (included by all Makefiles)

### External Library Builders
- `extlibs_makefile` - Build rules for DeviceXlib, LibMBD, Wannier90
- `oldlibs_makefile` - Legacy LAPACK and FoX build rules
- `plugins_makefile` - Third-party plugin build/install rules (GIPAW, D3Q, WANT, Yambo)

### Dependency Generation
- `makedeps.sh` - Script to generate `make.depend` files throughout the source tree

### Platform-Specific
- `README.FX10` - Notes for Fujitsu FX10 platform

## Usage
Called by the root Makefile and configure script. Not typically modified directly.
