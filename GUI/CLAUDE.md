# GUI - Graphical User Interfaces

## Purpose
Graphical tools for QE input file preparation, visualization, and editor support.

## Language
Tcl/Tk, shell scripts.

## Components

### PWgui (`PWgui/`)
Main graphical user interface for QE. Built on the Guib framework.
- Creates input files for pw.x, ph.x, pp.x, neb.x, projwfc.x, etc.
- Atomic coordinate editing
- k-point definition and visualization
- Input validation

### Guib (`Guib/`)
Generic GUI builder framework used by PWgui.
- Tcl-based framework for scientific application GUIs
- Widget system, input validation, help integration
- Reusable for other scientific codes

### QE-modes (`QE-modes/`)
Emacs major modes for editing QE input files.
- Syntax highlighting for QE input format
- Auto-indentation
- Keyword completion

## Building
```bash
make gui    # from root (requires Tcl/Tk)
```
PWgui is built in `GUI/PWgui/` with a link in `bin/pwgui`.

## Dependencies
- Tcl/Tk runtime (for PWgui and Guib)
- Emacs (for QE-modes)
