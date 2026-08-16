# Porting notes

## Scope

This port targets the computational/statistical implementation of
COMPoissonReg 0.8.2. R formula parsing, S3 printing/summary methods, Rcpp glue,
R data objects, and vignette/graphics infrastructure are intentionally omitted.

## Normalizing constant

The upstream hybrid rule is retained: use the asymptotic approximation when
`lambda^(-1/nu) < hybrid_tol`; otherwise truncate the infinite series using the
same Stirling/geometric tail bound. The Fortran summation advances terms by
recurrence rather than repeatedly evaluating `lgamma`, while retaining log-space
accumulation.

The exact `nu=0, lambda<1` geometric limit is handled explicitly for the
normalizer, CDF, quantile, mean, and variance. This is mathematically equivalent
to the limiting model and avoids an unnecessarily long truncation.

## Regression optimizer

R uses `stats::optim` (default `L-BFGS-B`) and its returned Hessian. The Fortran
port is standalone and therefore uses a BFGS maximizer with bounded trial-step
size. The regression score is analytical with respect to the link predictors:

- beta score uses `Y - E(Y)`;
- gamma score uses `nu * (E(log(Y!)) - log(y!))`;
- ZICMP zeta score uses the exact mixture score.

An independent central-difference Hessian of the log-likelihood is computed at
the final estimate and inverted to obtain covariance estimates. Fixed
coefficients and all three offset vectors are supported.

## Upstream source issues not propagated

Two issues are visible in the supplied R source:

1. `get.offset` constructs `list(x=x, s=s, w=s)` rather than using its `w`
   argument. The Fortran `cmp_offset_t` stores the requested `w` vector.
2. `NAMESPACE` registers `S3method(leverage,zicmpfit)`, and
   `deviance.zicmpfit` calls `leverage(object)`, but no
   `leverage.zicmpfit` implementation exists anywhere in the supplied source.
   The port therefore provides the implemented CMP leverage/deviance path and
   does not invent a ZICMP leverage formula.

ZIP/ZICMP densities also handle the exact `p=1` endpoint directly.

## Fisher information

The `fim.*` functions are translated even though upstream comments call them
potentially numerically unstable and do not export them. As upstream, exact FIM
variants use numerical derivatives of the log normalizing constant; Monte Carlo
variants use simulated scores.

## Vectorization

R automatically vectorizes and recycles scalar distribution parameters. The
Fortran scalar probability functions are designed to be called inside ordinary
or `do concurrent` loops. `rcmp_vec` and `rzicmp_vec` are provided for varying
per-observation parameters because those are important to regression/bootstrap.

## Random streams

Fortran `random_number` replaces R's RNG, so streams do not reproduce R seeds
bit-for-bit. Distributional tests are used instead.
