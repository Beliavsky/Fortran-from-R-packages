# GeneralizedHyperbolic-fortran

Modern Fortran/FPM translation of the numerical core of the R package
`GeneralizedHyperbolic` 0.8-7 by David Scott.

## Implemented

- Generalized inverse Gaussian (GIG): `dgig`, `pgig`, `qgig`, `rgig`
- GIG moments, mean, variance, mode and four parameterization conversions
- Generalized hyperbolic (GH): `dghyp`, `pghyp`, `qghyp`, `rghyp`
- GH mean, variance, mode, five parameterization conversions and rescaling
- Hyperbolic specialization: `dhyperb`, `phyperb`, `qhyperb`, `rhyperb`
- Normal inverse Gaussian specialization: `dnig`, `pnig`, `qnig`, `rnig`
- Skew-Laplace: `dskewlap`, `pskewlap`, `qskewlap`, `rskewlap`
- Likelihood fitting: `gig_fit`, `hyperb_fit`, `nig_fit`
- Hyperbolic-error linear-model entry point: `hyperblm_fit`

The GIG random generator is a direct modern-Fortran translation of the
Dagpunar rejection algorithm used by the upstream R package. GH random
variates use the same normal variance-mean mixture representation as upstream.

## Build

```sh
fpm test
```

or compile directly with a Fortran 2018 compiler.

## Scope

R plotting, QQ plots, S3 print/summary methods, bundled datasets, and formula /
data-frame interfaces are not translated. Numerical functions are exposed
through the `generalized_hyperbolic` module.

## License

GPL-2.0-or-later, matching the upstream package's `GPL (>= 2)` declaration.
