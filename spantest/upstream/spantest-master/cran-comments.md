## Submission

This is a feature and performance update of the 'spantest' package
(CRAN version 1.1-3 -> 1.3-0). See NEWS.md for the full list of changes; the
main ones are:

* New exported function `span_simulate()` that simulates benchmark and
  test-asset returns for mean-variance spanning size/power studies (normal,
  Student-t and skew-t innovations with i.i.d., AR(1), GARCH(1,1) or AR-GARCH
  dynamics).
* The Gungor-Luger Monte Carlo tests `span_gl_a()` / `span_gl_ad()` and the
  GARCH recursion used by `span_simulate()` now run in compiled code for speed;
  their results are unchanged (verified numerically). `span_as()` was also made
  faster through vectorised internal computations.

Because of the compiled code, the package now imports 'Rcpp' and links to
'RcppArmadillo' (LinkingTo). No user-facing function was removed and no
existing signature changed.

## Test environments

* local: macOS 15 (Apple silicon), R 4.5.2
* win-builder: R-devel and R-release
* macOS builder (mac-builder): R-release
* R-hub: Linux, Windows and macOS (R-devel)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies on CRAN.
