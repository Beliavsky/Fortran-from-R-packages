# Third-party notice

`actuar-fortran` is a translation of the computational portions of the R
package `actuar` by Vincent Goulet and contributors.

The upstream package declares `License: GPL (>= 2)` in `DESCRIPTION`.
Accordingly this translated package is distributed under GPL-2.0-or-later. A
copy of GPL version 2 is included as `COPYING`.

v0.2 vendors `expint-fortran`, the previously produced Fortran translation of
the R package `expint`, because upstream `actuar` imports `expint` for extended
incomplete-gamma calculations. The vendored dependency retains its own
`COPYING`, upstream metadata and provenance files under
`dependencies/expint-fortran/`.

Selected unmodified upstream actuar R and C source files are included under
`upstream/` solely for provenance and algorithm comparison and retain their
original copyright notices.
