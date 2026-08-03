# Porting notes

## R vectors

The three deterministic distribution functions are elemental Fortran
functions. Array syntax therefore replaces R vector recycling. Inputs must be
scalar or conformable arrays; nonconformable R-style recycling is not used.

## Invalid parameters

The upstream R code replaces negative or missing scale and shape values with
zero and then evaluates formulas containing division by zero. The Fortran API
instead treats `alpha <= 0` or `beta <= 0` as invalid. Deterministic functions
return IEEE NaN, while random generation returns `gnorm_invalid_argument`.

## Gamma calculations

The CDF uses series and continued-fraction evaluations of the regularized
incomplete gamma function. Quantiles use a safeguarded bracketed inverse. This
is dependency-free and accurate across the tested heavy- and light-tailed
parameter ranges, but it is not a copy of R's platform math library.

## Random generation

The identity

```text
(|X-mu|/alpha)^beta ~ Gamma(shape=1/beta, scale=1)
```

is used with a Marsaglia-Tsang gamma sampler and a random sign. A Park-Miller
integer generator provides a portable stream. A supplied seed is deterministic.

## Corrected upper-log quantile branch

For `lower_tail=.false.` and `log_probability=.true.`, the upstream R source
computes `log(1-p)` even though `p` is already a log probability. The Fortran
code uses `1-exp(p)`, consistent with R's documented distribution conventions.

## Omitted material

The R plotting examples, roxygen/S3 packaging, and vignette rendering are not
part of the compiled library. The original files remain under `original/`.
