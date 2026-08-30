# Notices and attribution

This directory is a modern free-form Fortran translation of the computational
algorithms in the R package **changepoint**, version **2.3** (2024-11-02),
originally published on CRAN on 2024-11-04.

Upstream package metadata lists:

- Rebecca Killick — author and maintainer
- Kaylea Haynes — contributor
- Harjit Hullait — contributor
- Idris Eckley — thesis advisor
- Paul Fearnhead — contributor and thesis advisor
- Robin Long — contributor
- Jamie Lee — contractor

The upstream `DESCRIPTION` declares `License: GPL` without a version qualifier.
That declaration is preserved verbatim in `upstream/DESCRIPTION`. Standard GPL
version 2 and GPL version 3 texts are included as `LICENSE-GPL-2.txt` and
`LICENSE-GPL-3.txt`; see `LICENSE.note` for the handling of the unversioned
upstream declaration.

The original package homepage is <https://github.com/rkillick/changepoint/> and
the CRAN package page is <https://CRAN.R-project.org/package=changepoint>.

The upstream package requests citation of, among other method-specific work:

Rebecca Killick and Idris A. Eckley (2014), "changepoint: An R Package for
Changepoint Analysis", *Journal of Statistical Software*, 58(3), 1-19.
<https://www.jstatsoft.org/article/view/v058i03>

The Fortran translation preserves algorithmic provenance from the upstream R
and C sources while replacing R objects, `.Call`/`.C` interfaces, S4 methods,
plotting, and dynamic R evaluation with typed Fortran APIs. No source from
`rfortran-core`, `rfortran-linalg`, BLAS, LAPACK, or another translated R
package is copied into this package.

The Fortran translation was prepared for the
`Beliavsky/Fortran-from-R-packages` collection. Translation-specific code is
provided subject to the rights and obligations applicable to the upstream GPL
work; this notice does not remove or replace upstream copyright or attribution.
