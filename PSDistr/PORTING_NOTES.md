# Porting notes

- Upstream PSDistr 0.0.1 is pure R and declares `License: GPL-3`.
- The only pracma primitive used by the package is `nthroot`. The supplied `pracma-fortran` tree is included under `vendor/` but the PSDistr library is intentionally self-contained; `signed_root(x,n)` handles the same real signed-root operation and also permits a real-valued shape exponent where the upstream formulas use one.
- `qen` is implemented by safeguarded bisection with adaptive bracketing rather than R's `uniroot` on a fixed large interval.
- `qeck` uses direct beta-quantile inversion instead of applying `uniroot` to the ECK CDF. This is algebraically equivalent and more stable near the support endpoints.
- SPC simulation avoids the upstream R expression that can transiently evaluate a fractional power of a negative normal variate before its branch is selected.
- Invalid parameter sets return IEEE quiet NaN rather than R character strings; this keeps the numerical API type-safe.
