# Provenance and attribution

This is a modified modern free-form Fortran translation of computational
operations from **posterior 1.7.0**, not an official upstream release.

- CRAN source mirror: https://github.com/cran/posterior
- Tested upstream revision: `0e696c89298974736164ffb6b738345b178cbf17`
- CRAN publication: 2026-04-01
- Upstream project: https://github.com/stan-dev/posterior
- Documentation: https://mc-stan.org/posterior/
- Target: https://github.com/Beliavsky/Fortran-from-R-packages
- Tested shared-repository revision: `cbd437d6d92f5481706fd8ea5a7829aaff7b6e24`
- Translation date: 2026-10-11

Upstream authors: Paul-Christian Bürkner, Jonah Gabry, Matthew Kay and Aki Vehtari.
Contributors: Måns Magnusson, Rok Češnovar, Ben Lambert, Ozan Adıgüzel,
Jacob Socolar, Noa Kallioinen and Teemu Säilynoja. Original author/maintainer
metadata, license notice, NAMESPACE and R source are retained verbatim under
`upstream/`; `upstream/SHA256SUMS` records their checksums. These R files are
provenance references, not runtime dependencies. No third-party dependency
source is embedded.

Upstream LICENSE identifies: YEAR 2021; COPYRIGHT HOLDER posterior package
authors; Stan Developers and their Assignees; Trustees of Columbia University.
The convergence source additionally preserves Copyright (C) 2012–2018 Trustees
of Columbia University and Copyright (C) 2018, 2019 Aki Vehtari, Paul Bürkner.
The unimplemented rstar provenance source retains Copyright (C) 2020 Ben Lambert,
Aki Vehtari. Full applicable BSD-3-Clause terms appear in `LICENSE`.

## Changes and dependency decisions

The Fortran translation uses native real arrays, explicit numerical outputs,
optional status integers and pure procedures. Rank normalization, ESS truncation,
Pareto-grid fitting and smoothing follow the pinned R algorithms. Autocovariance
uses the compatible direct shared implementation rather than R's FFT. Numerical
callbacks replace selected R expression interfaces; supplied uniforms replace
implicit R RNG in the implemented resampling/PIT paths. R objects, labels,
recycling beyond documented broadcasting, presentation and dynamic dispatch are
not reproduced. The coverage documents list all remaining computations.

Existing translations and shared modules were inspected before implementation.
The sibling `rfortran-core` provides the compatible dp kind, ordering, biased
autocovariance, type-7 finite quantiles, normal quantiles, incomplete beta and
stable log/exponential functions. It is an FPM path dependency with its own MIT
license, never copied into this directory. Endpoint handling around finite-core
quantiles is local because zero times infinity is not a valid exact-order-statistic
interpolation. The existing matrixStats translation was inspected; no dependency
on it is necessary because these operations are already provided by the focused
core API or Fortran intrinsics. There are no external BLAS/LAPACK links and no
vendored BLAS, LAPACK, ARPACK, R compatibility module or package dependency.

The algorithms cite Vehtari et al. (2021), Bayesian Analysis 16, 667–718,
doi:10.1214/20-BA1221; Zhang and Stephens (2009), Technometrics 51, 316–325;
Vehtari et al. (2024), Pareto-smoothed importance sampling; and the references
preserved in the original R source for nested Rhat and uniformity tests.

Optimized-build roundoff guard: scalar and vectorized EXP can differ by one
ULP at equal log-tail/cutoff values. Log-weight tail excesses are constrained
to their mathematically nonnegative domain before fitting, preventing a tied
cutoff from becoming a spurious negative excess. The pinned R smoothing
references pass at unchanged tolerances in debug and release builds.
