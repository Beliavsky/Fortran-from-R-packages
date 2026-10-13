# License, attribution, and provenance

This translation is licensed under the GNU General Public License, version 2
or (at your option) any later version. The full GPL version 2 is in LICENSE.
There is no warranty, including no implied warranty of merchantability or
fitness for a particular purpose.

Upstream: exactRankTests 0.8-37, date 2026-05-02, CRAN publication
2026-05-04. Authors: Torsten Hothorn and Kurt Hornik. Maintainer:
Torsten.Hothorn@R-project.org. Source mirror:
https://github.com/cran/exactRankTests at commit
`dc41a275dfd99c46b43c5e9138dd040410343617`.
CRAN package: https://cran.r-project.org/package=exactRankTests

The original native distribution implementation carries
Copyright (C) 2000-2003 Torsten Hothorn; its simulation implementation carries
Copyright (C) 2003 Torsten Hothorn. The original irank implementation says it
was modified from R-1.5.1/src/main/sort.c and refers to that file's copyright;
Copyright (C) 1995, 1996 Robert Gentleman and Ross Ihaka and Copyright (C)
1998-2002 The R Development Core Team apply to that source. Its complete
original leading notice is retained in upstream/R-sort-COPYRIGHT.txt. Relevant original R and
C files, DESCRIPTION, and NAMESPACE are retained under upstream/ as provenance,
with their original notices intact. They are not compiled into the FPM library.
No upstream dependency package is copied.

Fortran translation and tests: OpenAI-assisted contributors, 2026, distributed
under GPL-2.0-or-later. Revisions include pure APIs, caller-managed uniform
variates, probability-scaled dynamic programming, GCD lattice compression,
explicit statuses and bounded allocation, validated sorting orders, and IEEE
handling instead of R object/warning interfaces.

Algorithm attribution: Bernd Streitberg and Joachim Roehmel's Shift-Algorithm,
including the references in upstream/src/permdist.c: Statistical Software
Newsletter 12(1), 10-17 (1986); EDV in Medizin und Biologie 18(1), 12-19 (1987).

Shared dependency: rfortran-core by Vivek Rao and contributors, MIT licensed,
used through ../rfortran-core at repository commit
`cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`. The existing dependency supplies its
own license and source; this package does not duplicate it. No BLAS/LAPACK,
ARPACK, rfortran-compat, or external binary dependencies are required.
