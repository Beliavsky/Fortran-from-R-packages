## Resubmission
This is a resubmission for 'optimflex' version 0.1.7.
The package was previously submitted as version 0.1.6.

## Changes since 0.1.6
* Fixed bugs in `dogleg()` and `double_dogleg()` so that when the user 
  supplies an exact Hessian, it is used directly in the optimization 
  rather than being overridden by the BFGS approximation.

## Test results
* R CMD check --as-cran: 0 errors | 0 warnings | 1 note
* Note: "unable to verify current time" — a local network issue, not a 
  package problem.
* Win-builder (R-release): OK
* Local Windows 11 (x86_64): OK