# NOTICE and attribution

## Upstream R package

This project is a source translation of **ranger 0.18.0**, "A Fast Implementation of Random Forests".

Upstream package metadata lists:

- Marvin N. Wright — author and maintainer
- Stefan Wager — contributor
- Philipp Probst — contributor

The upstream R package declares `License: GPL-3`. The Fortran package therefore retains GPL-3 as its package license; see `LICENSE`.

Upstream citation:

Marvin N. Wright and Andreas Ziegler (2017), "ranger: A Fast Implementation of Random Forests for High Dimensional Data in C++ and R", *Journal of Statistical Software* 77(1), 1–17, DOI 10.18637/jss.v077.i01.

The original `DESCRIPTION`, `NAMESPACE`, `NEWS`, and `inst/CITATION` are retained under `upstream/` for provenance.

## Separate upstream C++ core notice

The upstream C++ core contains the following licensing statement in `src/globals.h`:

> Copyright (c) [2014-2018] [Marvin N. Wright]
>
> This software may be modified and distributed under the terms of the MIT license.
>
> Please note that the C++ core of ranger is distributed under MIT license and the R package "ranger" under GPL3 license.

The MIT license text for that upstream core is retained in `LICENSE.MIT-CORE`. The Fortran translation combines algorithmic material originating in that core with GPL-3 R-level algorithms and newly written Fortran glue, and is distributed as the GPL-3 package described above.

## Other algorithm provenance retained from upstream

The infinitesimal-jackknife code in upstream `R/infinitesimalJackknife.R` states that it is a ranger adaptation of `randomForestCI` by Wager et al. and identifies Stefan Wager as author of that implementation. This translation preserves the raw IJ, Monte Carlo bias-correction, and finite-sample formulas and records that provenance here.

The hierarchical-shrinkage interface cites Agarwal et al. (2022) in the upstream documentation; the translated recursion follows upstream `src/utilityRcpp.cpp`.

## Translation policy

No Rcpp, RcppEigen, Matrix, BLAS, LAPACK, or copied translated-R-package source is embedded here. `rfortran-core` is referenced as a sibling FPM dependency solely for the shared `r_kinds::dp` kind.
