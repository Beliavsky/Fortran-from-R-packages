## Resubmission

This is a resubmission of 1.7-0. The previous submission was returned by the
auto-check service with one NOTE on r-devel-windows-x86_64:

    checking compiled code ... NOTE
    Error in ccE(lines, flags = new_flags, include = include) :
      'cc' is not on the path
    Calls: <Anonymous> ... lapply -> FUN -> lapply -> FUN -> getFunsHdr -> ccE
    Execution halted

We believe this is not a finding about the package but the check failing to
start: `tools:::ccE()` begins with `if (Sys.which("cc") == "") stop("'cc' is
not on the path")`, so it exited before inspecting any object file. The
Debian flavour of that same submission reported OK on every check, including
"checking compiled code", and the check completes normally on our own
machines. We could not find anything to change in the package in response to
it; if the maintainers would prefer a package-side change, we will gladly
make one.

This resubmission does contain unrelated improvements found during a further
audit, which is why the tarball differs from the previous one:

* The optimizer no longer reports a spurious convergence code 52 when the
  local polish step starts from an already-converged point (this affected the
  printed output of the documented examples).
* Standard errors are no longer withheld from legitimate high-curvature fits:
  the guard against a penalty-contaminated numerical Hessian now tests the
  finite-difference stencil directly instead of using a fixed magnitude
  cutoff, which does not scale with the sample size.
* Portfolio and forecasting inputs are validated (positive volatilities,
  positive-definite correlation matrices, finite future covariates), and
  `N` must be a whole number.
* New regression tests for all of the above.

## Submission

This is an update of the package already on CRAN (current version 1.1-2,
published 2025-09-03) to version 1.7-0. The main changes since 1.1-2 are:

* Datasets `greenbrown`, `mccc`, and `ff5ind` (used by the examples and the
  companion article).
* A C++ (Rcpp/RcppArmadillo) Hamilton filter, validated against a pure-R
  reference to ~1e-8.
* An `"rsdc_fit"` S3 interface (print/summary/coef/logLik/nobs/vcov/confint/
  predict/simulate/plot), numerical and bootstrap standard errors, multi-step
  forecasts, uncertainty bands, Viterbi decoding, and broom/ggplot2 methods.
* A canonical partial-correlation (Joe, 2006) reparameterization and a
  reparameterized global search that makes maximum-likelihood estimation
  feasible at higher cross-sectional dimensions, plus `rsdc_starts()` warm
  starts and a multi-start replication diagnostic.

See NEWS.md for the full history.

## Test environments
* Local: macOS Tahoe 26.5, R 4.5.2, aarch64 (source tarball, --as-cran):
  0 errors | 1 warning | 0 notes
* win-builder, R devel (2026-07-30 r90327 ucrt), x86_64-w64-mingw32:
  0 errors | 0 warnings | 1 note

## R CMD check results
No errors. The one warning and the one note are described below; neither is
raised by package code, and they do not overlap (each environment shows only
its own).

### win-builder R-devel: "checking compiled code ... NOTE"

The check did not report a finding about the package: it failed to run. The
note contains an R traceback ending in

    Error in ccE(lines, flags = new_flags, include = include) :
      'cc' is not on the path

`tools:::ccE()` begins with `if (Sys.which("cc") == "") stop(...)`, so this
is a missing prerequisite on that build machine rather than a diagnostic
about the compiled code. The same check completes and reports OK on the
local macOS build. Everything else on win-builder R-devel passed, including
`checking whether package 'RSDC' can be installed`, `checking compilation
flags in Makevars`, `checking for portable use of $(BLAS_LIBS) and
$(LAPACK_LIBS)`, the tests, and the vignette re-build.

### Local macOS: "checking whether package 'RSDC' can be installed ... WARNING"

The warning is not raised by package code. It comes from R's own header
`R_ext/Boolean.h`, which suppresses `-Wfixed-enum-extension` for Apple clang
via `#pragma clang diagnostic ignored`. On this machine the C++ compiler
(Apple clang 21.0.0) no longer recognises that warning group, whereas the R
installation was built with Apple clang 16.0.0, so the pragma itself triggers
`-Wunknown-warning-option`. The package's `src/` contains no pragmas and no
enum declarations; the diagnostic disappears when the compiler matches the
one R was built with, as on CRAN's build machines.

## Reverse dependencies
There are no reverse dependencies on CRAN.

## Notes for the maintainer
* The package contains compiled code (Rcpp/RcppArmadillo). `src/Makevars` and
  `src/Makevars.win` link BLAS/LAPACK explicitly, as RcppArmadillo requires;
  win-builder confirms the Windows build installs cleanly and that the use of
  `$(BLAS_LIBS)`/`$(LAPACK_LIBS)` is portable.
* Examples that estimate a model are wrapped in \donttest{} to keep check times
  short; the full suite runs locally and on the environments above.
