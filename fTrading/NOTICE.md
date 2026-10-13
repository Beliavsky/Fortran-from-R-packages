# Notice and provenance

This package translates computational code from **fTrading 3042.79**,
"Rmetrics - Trading and Rebalancing Financial Instruments", dated 2017-11-12,
published by CRAN on 2017-11-15. Authors and contributors identified by upstream:
Diethelm Wuertz (author and numerical function implementer), Tobias Setz
(maintainer), and Yohan Chalabi (contributor). The benchmark functions identify
A. Trapletti's tseries package as their source; Yohan Chalabi's 2008 Sharpe
revision is retained. Original comments, licensing, and authorship notices
remain unchanged in `upstream/`. The upstream Rmetrics copyright is
(C) 1999-2005 Diethelm Wuertz; the complete multi-project notice is preserved
in `upstream/inst/COPYRIGHT.html`.

- Official package: https://cran.r-project.org/package=fTrading
- Official source: https://cran.r-project.org/src/contrib/fTrading_3042.79.tar.gz
- Upstream project: http://www.rmetrics.org
- Retrieved 2026-10-12 from the official CRAN source collection.
- SHA-256: 3e99963eb59891b14a1531f8b30692c59ae0e2a97a6d63c3ad3b37772d79ce76
- MD5: 2063faf646712006e77de73786ded556
- Audited target repository revision: cbd437d6d92f5481706fd8ea5a7829aaff7b6e24
  (https://github.com/Beliavsky/Fortran-from-R-packages).

Upstream DESCRIPTION declares **GPL (>= 2)**. Each substantive R source file
contains a **GNU Library General Public License version 2 or later** notice.
Both statements are preserved rather than silently normalized. `LICENSE`
contains GPL version 2 and `LICENSES/LGPL-2.0-or-later.txt` contains GNU Library
GPL version 2. Translated Fortran source follows the substantive source's
LGPL-2.0-or-later notice; this assembled package also retains the package-level
GPL-2.0-or-later declaration and uses the GPL-2.0-or-later TTR dependency.

Changes in this translation: modern free-form Fortran 2008 array API,
FPM sibling dependencies, explicit precision and purity, safe invalid-input
results, omission of R classes/dates/names/plotting, deterministic fixtures and
independent numerical tests. The README documents semantic quirks, adaptations,
and repaired vector-interface bugs. No warranty is provided.

Dependencies are referenced only, never copied: TTR 0.24.4 (Joshua M. Ulrich
and contributors, GPL-2.0-or-later) and rfortran-core (its existing repository
contributors, MIT). Their own notices and source remain in their sibling
packages. No BLAS, LAPACK, ARPACK, r.f90, r_mod.f90, or translated R dependency
implementation is embedded here.
