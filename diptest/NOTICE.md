# Notices and licensing

This package translates the computational code of R package **diptest 0.77-2**.
The upstream package declares **GPL (>= 2)**. The translated code, reproduced
diptest material, table data and newly authored tests are distributed under
**GNU GPL version 2 or, at your option, any later version**. The mixed oracle
fixture containing R datasets retains R's **GPL version 2 or version 3** terms;
see the component notice below. See LICENSE for the complete GPL version 2 text
and licenses/GPL-3.txt for version 3. No warranty is provided.

## Upstream authors and provenance

The upstream DESCRIPTION identifies:

- **Martin Maechler**, author and maintainer,
  <maechler@stat.math.ethz.ch>, ORCID 0000-0002-8685-9910; ported to R and C from
  Fortran and S-plus and substantially enhanced the implementation.
- **Dario Ringach**, author of the Fortran and S-plus interface.

The algorithm is **AS 217**, by **P. M. Hartigan** (1985), based on work by
**J. A. Hartigan and P. M. Hartigan**. Upstream `src/dip.c` records translation
from `dip.f` by f2c, version of 22 July 1992 22:54:52, and editing/extension by
Martin Maechler at ETH Seminar fuer Statistik, 8092 Zurich, Switzerland.
The complete historical header is preserved in `upstream/src/dip.c`, including:

1. Martin Maechler's July 30, 1994 fix to terminate the cycling for unimodal data.
2. The October 31, 2003 correction attributed to **Yong Lu**
   <lyongu+@cs.cmu.edu> of a misplaced parenthesis in the original Fortran.

Regression cases attributed upstream to **Ferenc Mechler**
(<fmechler@med.cornell.edu>, September 5, 2002) are retained in the tests and
example. This translation uses the corrected current C computation, not the
historical erroneous Fortran program.

The `qDiptab` quantile table is by **Martin Maechler**, in its earliest form from
August 1994. The current table derives from 1,000,001 uniform samples per sample
size. It is distributed as part of GPL-licensed diptest; no separate table
license or notice is supplied upstream. Both the exact original RDS and an
unrounded decimal Fortran representation are included. The source R
`dip.test` interpolation, including averaged tied knots, is translated.

## Translation changes

The 2026 translation introduces modern free-form Fortran, pure typed APIs,
explicit statuses and input contracts, local binary-scaled geometry for extreme
ranges, defined constant-sample diagnostics, caller-supplied null simulations,
an intrinsic-RNG streaming convenience procedure, documentation and deterministic
tests. Source headers and this notice identify the upstream origin.

The `rfortran-core` package is a sibling FPM dependency reused for `dp` and
sorting. Its separate MIT license and attribution remain in that package; no
copy of its source is embedded here. No BLAS, LAPACK, R runtime or other translated
R-package dependency is copied or linked.

## References

- P. M. Hartigan (1985), “Computation of the Dip Statistic to Test for
  Unimodality,” Applied Statistics 34, 320–325. Algorithm AS 217.
- J. A. Hartigan and P. M. Hartigan (1985), “The Dip Test of Unimodality,”
  Annals of Statistics 13, 70–84.
- R package source: <https://CRAN.R-project.org/package=diptest>
- Development source and history: <https://github.com/mmaechler/diptest>

Test data from R's datasets package are attributed separately in
`test/data/README.md`. R's distribution-wide copyright notice is reproduced
unchanged in `licenses/R-COPYRIGHTS.txt`; it states that R as a whole is under
GPL version 2 or 3. The datasets metadata identifying its R Core Team and
worldwide contributors is preserved in `test/data/R-datasets-DESCRIPTION`.
The full R copyright notice mentions components that are not included here;
its inclusion does not imply those sources are bundled. These dataset terms
are not changed by the GPL-2.0-or-later license of the Fortran implementation.
