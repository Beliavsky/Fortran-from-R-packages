# Provenance, attribution and license

This package translates computational code from the R package **entropy 1.3.2**,
package date 2025-04-08, authored by **Jean Hausser** and **Korbinian Strimmer**.
The package declares `GPL (>= 3)`. Translated code and new test/example code are
released under **GNU GPL version 3 or any later version**. The complete GPL v3
text is supplied in `LICENSE`. There is no warranty, including implied
merchantability or fitness for a particular purpose.

Pinned source:
- CRAN mirror repository: https://github.com/cran/entropy
- Commit: `8921581fd845025b5d56c9f27ed1c23c334ac092`
- Source tree: https://github.com/cran/entropy/tree/8921581fd845025b5d56c9f27ed1c23c334ac092
- Project: https://strimmerlab.github.io/software/entropy/
- Destination repository: https://github.com/Beliavsky/Fortran-from-R-packages
- Translation date: 2026-10-11

The unmodified upstream `DESCRIPTION`, `NAMESPACE`, `NEWS`, and all `R/*.R`
files are retained under `upstream/` as the source and reproducibility record.
Their original source notices are preserved in full. Applicable source-file
copyright notices are:

- `Gstat.R`: Copyright 2021 Korbinian Strimmer
- `KL.plugin.R`, `discretize.R`: Copyright 2013 Korbinian Strimmer
- `entropy.ChaoShen.R`, `entropy.MillerMadow.R`, `entropy.plugin.R`:
  Copyright 2008 Korbinian Strimmer
- `entropy.Dirichlet.R`, `entropy.R`, `entropy.empirical.R`,
  `entropy.shrink.R`, `mi.plugin.R`: Copyright 2008-13 Korbinian Strimmer
- `entropy.NSB.R`: Copyright 2006-2008 Jean Hausser

New Fortran translation, documentation and validation material:
Copyright 2026 entropy Fortran contributors, GPL-3.0-or-later.
The archival NSB interface and its Jean Hausser notice are retained despite NSB
being untranslated; no NSB executable or replacement algorithm is included.

## Shared dependency

`rfortran-core` is used through `../rfortran-core`, not copied into this package.
It supplies `r_kinds:dp`, `r_transforms:r_log1p/r_expm1` and
`r_distributions:r_pchisq`. It is a separately maintained MIT-licensed package by
Vivek Rao and contributors; consult its own LICENSE and provenance. No LAPACK,
BLAS, R runtime, translated dependency source, or external numerical executable
is bundled or linked.

## Numerical references cited upstream

The numerical operations follow the pinned R implementations: maximum-likelihood
and Dirichlet plug-in estimates, Miller–Madow (1955) bias correction,
Chao–Shen (2003) coverage correction, and Hausser–Strimmer (2009)
uniform-target shrinkage. See upstream package documentation and source for
original bibliographic context. Fortran adaptations and compatibility limits
are explicitly recorded in README.md and each coverage mapping in fpm.toml.
