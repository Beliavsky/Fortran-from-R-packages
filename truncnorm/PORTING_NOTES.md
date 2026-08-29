# Porting notes

## Translation boundary

`truncnorm` is almost entirely computational native code. The Fortran port
therefore retains every exported numerical routine and all native sampling and
root-finding algorithms. R's `.Call` registration, SEXP allocation/macros,
interrupt checks and R documentation objects are interface code and are not
recreated.

## `r_mod` reuse

The user-supplied `r_mod.f90` is retained verbatim under
`upstream/r_mod-original.f90`. The build copy `src/r_mod.F90` differs only in
free-form continuation/line wrapping required for standard line lengths. A
whitespace/continuation-normalized comparison is identical.

The port uses `r_mod` for Normal density/CDF/quantiles, uniform/Normal RNG and
seeding. It does not reimplement those helpers.

The upstream C implementation relies on R Mathlib's `pnorm(..., log=TRUE)` and
`logspace_sub`. `r_mod` has no log-tail Normal CDF/SF or `logspace_sub`, so a
small local tail-stable helper layer is genuinely necessary. This is used in
expectation and quantile calculations; ordinary CDF behavior remains
consistent with the upstream parameterization.

## Extreme-tail behavior

The upstream package deliberately changes behavior when a finite truncation
interval lies more than six standard deviations beyond the mean:

- `etruncnorm` returns the interval midpoint;
- `vtruncnorm` returns the variance of a uniform distribution on the interval;
- `dtruncnorm` can fall back to a uniform density if the normalizing CDF
  difference is no longer representable in ordinary floating point.

These behaviors are preserved rather than silently replacing them by a more
accurate asymptotic truncated-normal calculation. This matters for parity with
the package's own regression tests, including `(a,b,mean,sd)=(0,1,5,0.1)`.

The CDF/quantile implementation uses stable log interval probabilities when
possible. This only extends usability in very remote tails; it does not alter
the documented finite-region formulas or the explicit extreme mean/variance
fallbacks.

## Quantile root finding

`src/zeroin.c` comes from R's main distribution and NETLIB Brent code. One
important C detail is preserved literally: the apparent three-way swap uses
sequential assignments. Treating it as a simultaneous swap changes the
algorithm and caused endpoint sticking in early translation tests.

## Random generation

The sampler-selection logic and thresholds are direct translations of
`src/rtruncnorm.c`. R Mathlib's exponential RNG accepts a scale argument;
`r_mod`'s exponential helper uses a rate. The Fortran code therefore samples
`-log(1-U)/a` directly in the exponential-rejection branches, which is exactly
an exponential variate with rate `a` and corresponds to the upstream
`rexp(1/a)` scale call.
