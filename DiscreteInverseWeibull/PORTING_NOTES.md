# Porting notes

## Rsolnp

Method `M` uses the supplied `Rsolnp-fortran` package through its callback/data API, with the upstream bounds `q in [0,1]` and `beta in [2,100]`. The objective returns a finite large penalty at invalid/singular boundary points; mathematically the second moment requires `beta > 2`.

## Numerical stability

The PMF and survival differences are evaluated in log/stable-difference form to avoid cancellation for large observations.

## Upstream `heuristic` bug

The supplied R function begins by assigning `beta1 <- 1`, which makes its documented `beta1` argument ineffective. The Fortran translation honors the caller-supplied `beta1`. A maximum iteration count is also provided to prevent an unbounded loop on pathological samples.

## PP method

The probability-plot method is computational, not plotting code, so it is retained. The R `ecdf` and `lm` calls are replaced by direct empirical-CDF and ordinary least-squares calculations.

## RNG

The inverse-transform RNG is preserved but uses the Fortran intrinsic random generator, so seeded streams do not match R exactly.
