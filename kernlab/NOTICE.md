# Notices

This project is a modern Fortran translation of the computational code in
`kernlab` 0.9-33.

The upstream package declares `GPL-2` and identifies the original R and native
code authors in `original/kernlab-master/inst/COPYRIGHTS`. The translated
Fortran source is distributed under GPL-2.0-only to preserve that licensing.

The upstream source contains additional components with their own notices,
including string-kernel/suffix-array code. Those files are retained under
`original/`; the compiled Fortran library uses a new portable spectrum-string
kernel and does not compile the retained C/C++ implementation.
