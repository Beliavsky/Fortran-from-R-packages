# License and provenance

This translation is derived from the R package `leaps` version 3.2, whose
DESCRIPTION declares `License: GPL (>= 2)`.  Accordingly, this translated
package is distributed under the GNU General Public License, version 2 or
(at your option) any later version (`GPL-2.0-or-later`).  A copy of GPL-2.0
is included in `licenses/GPL-2.0.txt`.

The numerical algorithms are based on Fortran code by Alan Miller, as stated
in the upstream package README and source headers.  The supplied free-format
Fortran 95 modules `lsq.f90` and `find_sub.f90` retain Alan Miller's original
attribution and revision history; translated/modernized copies are compiled
as `src/leaps_lsq.f90` and `src/leaps_find_subsets.f90`.

For traceability, the relevant original R and Fortran source files are kept
verbatim under `upstream/`.  Those files are provided under their original
upstream licensing terms; the package-level GPL-2.0-or-later terms apply to
this derivative distribution.
