# Porting notes

## Scope

This port translates the full numerical content of ZeroOneDists 1.0.0:
all five d/p/q/r distribution sets and the numerical components embedded in the
five GAMLSS family constructors. R S3 objects, expression construction,
formula handling, and residual plotting are not translated.

## Numerical implementation

- Beta density/CDF support is implemented internally using `log_gamma` and a
  continued-fraction regularized incomplete beta calculation.
- BER quantiles use safeguarded Newton/bisection on the exact mixture CDF,
  replacing R `uniroot` while solving the same equation.
- Normal probabilities use the Fortran `erfc` intrinsic; inverse-normal values
  use bracketed inversion/Newton refinement.
- UMB quantiles use bisection on the translated CDF.
- UHLG and UPHN retain their closed-form quantiles.
- BER/BER2 family scores use central numerical derivatives with upstream step
  `1e-5`; UPHN uses the upstream step `0.01`.
- UHLG and UMB retain their analytic first derivatives. UMB also retains its
  analytic second derivative.
- For BER, BER2, UHLG, and UPHN the family working Hessian reproduces the
  upstream `-score_i*score_j` construction with the `-1e-15` floor.

## BER2 boundary correction

The upstream helper computes

```text
theta = nu * (1 - abs(2*mu - 1))
base_mu = (mu - theta/2) / (1 - theta)
```

At `mu=0.5, nu=1`, this gives `0/0`, although the mixture weight is one and the
limiting distribution is exactly Uniform(0,1). The Fortran port treats this
case explicitly as Uniform(0,1). This also preserves the stated BER2 property
`E(X)=mu` at that boundary.

## Boundary handling

Fortran numerical functions return IEEE NaN for invalid parameters instead of
raising an R exception. A few top-level R routines admit mathematically
singular boundary parameters (for example `sigma=0` in `dBER` because the R
check uses `< 0`); the Fortran routines require valid positive beta shapes.
For CDF/quantile routines, conventional support endpoints are handled directly
where that is numerically meaningful.

## Dependency translations

The user supplied:

- `gamlss-fortran-v0.9.0`
- `gamlss.dist-fortran-v0.3.0`

Both FPM manifests declare GPL-3.0-only. The upstream ZeroOneDists source is
MIT. The supplied GAMLSS fitter also uses a closed integer family registry, so
linking it would require modifying/forking that GPL dependency merely to add
five family IDs. The active ZeroOneDists FPM library therefore remains
self-contained and MIT licensed while reproducing the relevant GAMLSS link and
family callback semantics directly.
