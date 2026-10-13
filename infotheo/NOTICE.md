# License, attribution, and upstream provenance

This is a modern free-form Fortran translation of R package **infotheo 1.2.0.1**,
by **Patrick E. Meyer**, maintainer `software@meyerp.com`. Upstream DESCRIPTION
specifies **GPL (>= 3)**. The translated implementation and test/support work
are provided under **GPL-3.0-or-later**. The complete GNU GPL version 3 text is
in `LICENSE`; original source and documentation metadata are retained below.
No warranty is provided. This is a translation, not an upstream release.

The upstream discretize documentation credits **Patrick E. Meyer, Frederic
Lafitte, Gianluca Bontempi, and Korbinian Strimmer**. Those attributions and all
original DESCRIPTION, NAMESPACE, CHANGES, R functions, C++ sources, header, and
manual pages are preserved under `upstream/`. The upstream files contain no
additional standalone copyright notice; none has been replaced or invented.

- CRAN package: https://cran.r-project.org/package=infotheo
- Exact source archive: https://cran.r-project.org/src/contrib/infotheo_1.2.0.1.tar.gz
- Archive SHA-256: `c0fb8ec97ad3a49f231c4c993b5eee70c6a61c8c30dc4a46197867e4763a29d4`
- CRAN mirror source: https://github.com/cran/infotheo
- Mirror commit: `a0a3450b5ed66f49fadb98a0509baeb7e78f167d`
- Upstream package date: 2014-07-16; CRAN publication: 2022-04-08
- Exact archive, R sources and C++ sources verified during translation on
  2026-10-12; the preserved R/C++ files match the official archive byte-for-byte.
- Target repository: https://github.com/Beliavsky/Fortran-from-R-packages
- Dependency inventory inspected at repository commit
  `cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`.

The computational translation is in `src/infotheo.f90`. It translates the
probability counting, entropy estimators, information combinations, and
numerical discretization from `upstream/src/entropy.cpp`,
`upstream/src/discretize.cpp`, and `upstream/R/*.R`. The R registration,
data-frame machinery, and object presentation were not ported. Deliberate
numerical corrections, safe input handling, and a named SG compatibility mode
are described in README.md. No compiled upstream implementation is included.

## Shared dependency

Before implementing numerical helpers, the target repository was checked for
existing translations and shared APIs. No top-level `infotheo` or `entropy`
translation existed in the inspected repository. The compatible **rfortran-core
0.2.0** APIs are reused through the sibling path dependency `../rfortran-core`:

- `r_kinds: dp`
- `r_special: r_digamma`
- `r_ordering: r_order`
- `r_sorting: r_sort`

rfortran-core is MIT-licensed, copyright (c) 2026 Beliavsky, with manifest author
Vivek Rao and contributors. Its source and license are **not copied into this
package** and remain in that dependency's own directory. Preserve its MIT notice
when redistributing a combined dependency build. The GPL-3.0-or-later
translation may depend on the MIT library. No rfortran-compat, translated
entropy dependency, BLAS, LAPACK, or ARPACK implementation is copied, embedded,
or linked as a system library.
