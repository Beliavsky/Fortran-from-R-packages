# gamlss.dist-fortran

Modern Fortran 2018 translation of computational routines from the R package
`gamlss.dist` 6.1-1, packaged for the Fortran Package Manager (FPM).

The project concentrates on reusable probability kernels and likelihood fitting.
R plotting, `gamlss.family`/S3 objects, formulas, model frames, residual plotting,
and other R-specific infrastructure are intentionally omitted.

## v0.3.0 coverage

v0.3.0 retains the v0.1/v0.2 distribution catalog and completes the major
families identified as the previous release targets.

New continuous families and variants:

- `ST3C`
- skew normal `SN1`, `SN2`
- standardized skew-t `SST`
- generalized-t `GT`
- ex-Gaussian `exGAUS`
- `PARETO`, `PARETO1`, `PARETO2`, `PARETO2o`
- flexible gamma `GAF`

New discrete/count families and variants:

- double binomial `DBI`
- `PIG2`
- flexible negative binomial `NBF`
- zero-inflated flexible negative binomial `ZINBF`
- `ZIPIG`, `ZAPIG`
- `ZISICHEL`, `ZASICHEL`
- `ZIBB`, `ZABB`
- `ZIBNB`, `ZABNB`
- `ZAZIPF`

Most translated families provide `d*`, `p*`, `q*`, and `r*` routines.
`DBI` uses exact finite-support normalization rather than an asymptotic
normalizing approximation.

## GAMLSS-style likelihood fitting

`fit_gamlss` provides numerical maximum-likelihood regression with independent
design matrices for `mu`, `sigma`, `nu`, and `tau`. v0.3.0 exposes 62 family
constants through this generic interface.

v0.3 additionally adds fixed-denominator fitting routines where the denominator
cannot be represented by the generic four-parameter interface:

- `fit_dbi`
- `fit_zibb`
- `fit_zabb`

Results include coefficients, fitted distribution parameters, log likelihood,
AIC, convergence state, and numerical-Hessian covariance matrices.

## Example

```fortran
program v03_remaining
   use gamlss_dist
   implicit none

   print '(a,f12.6)', 'GT 90% quantile: ', &
      qGT(0.90_dp,mu=0.0_dp,sigma=1.0_dp,nu=4.0_dp,tau=1.5_dp)
   print '(a,f12.6)', 'Double-binomial P(Y=3): ', &
      dDBI(3.0_dp,mu=0.35_dp,sigma=0.6_dp,bd=8.0_dp)
end program v03_remaining
```

With FPM:

```text
fpm test
fpm run --example basic
fpm run --example v02_extended
fpm run --example v03_remaining
```

## Parameterization fidelity

The API keeps GAMLSS parameter meanings rather than substituting textbook
parameterizations. Examples include:

- `NBI`: negative-binomial size = `1/sigma`
- `NBII`: negative-binomial size = `mu/sigma`
- `NBF`: local overdispersion = `sigma * mu**(nu-2)`
- `GAF`: local CV parameter = `sigma * mu**(nu/2-1)`
- `WEI2`: Weibull scale = `mu**(-1/sigma)`, shape = `sigma`
- `BE`: beta shapes are reconstructed from GAMLSS `(mu,sigma)`
- `BCCG`, `BCT`, `BCPE`: positive-support truncation normalization is retained
- `SICHEL`/`SI`: upstream Bessel-ratio recurrences are retained
- `ST3C`: numerically equivalent to `ST3`; the upstream C acceleration layer is
  unnecessary in the native Fortran implementation

See `docs/API_MAP.md`, `docs/TRANSLATION_NOTES.md`, and
`docs/VALIDATION.md` for details.

## License and provenance

The upstream `gamlss.dist` package is licensed `GPL-2 | GPL-3`. This combined
Fortran work is GPL-3.0-only because reusable numerical support is adapted from
GPL-3 material in the earlier VGAM Fortran translation. Upstream metadata,
compiled reference code, and selected R sources used for translation checks are
retained under `upstream/`. See `THIRD_PARTY_NOTICE.md` and `LICENSES/`.
