# dev-tools - Developer Utilities

## Purpose
Collection of developer tools for code maintenance, documentation generation, memory analysis, and code style management.

## Language
Mixed: Tcl, Python, Perl, Shell scripts.

## Key Tools

### Documentation Generation
- `helpdoc` - Generates HTML/text documentation from `.def` input definition files. Used to create `INPUT_*.html` and `INPUT_*.txt` files for all QE programs.
- `README.helpdoc` - Documentation for helpdoc usage

### Memory Analysis
- `mem_counter` - Tracks memory allocation/deallocation in Fortran code
- `mem_analyse.py` - Python script for detecting memory leaks from allocation logs

### Code Style
- `src-normal` - Fortran code style normalizer (standardizes formatting)

### Call Graph Analysis
- `callhtml.pl` - Generates HTML call graph from Fortran source
- `calltree.pl` - Generates call tree from Fortran source

### GUI Maintenance
- `check_gui` - Validate GUI input definitions against code
- `update_gui_help` - Update GUI help text from helpdoc output

### Emacs Support
- `gen-emacs-mode` - Generates Emacs major mode definitions for QE input files

## Usage
These tools are for QE developers, not end users. They are not built as part of normal compilation.
