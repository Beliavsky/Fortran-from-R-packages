# statmod-fortran

Modern Fortran/FPM translation of the computational core of the R package
**statmod 1.5.2**.

## Highlights

The port includes inverse-Gaussian DPQR, Gaussian quadrature, Digamma-family
mathematics, expected unit-deviance approximations, secure Gamma and
negative-binomial GLM fits, `fitNBP`, two-component mixed models,
heteroscedastic REML scoring, ELDA, forward selection, robust scale and exact
statistical tests, growth-curve permutation tests, GLM score tests, Tweedie
family kernels and randomized quantile residuals.

The full native expected-deviance Chebyshev tables are retained. The optional
`qres.tweedie` path is supported by a vendored translated Tweedie numerical
core. Plotting, printing, S3/formula/model-frame infrastructure and other
non-computational R UI code are intentionally omitted.

See `API_MAPPING.md` for the R-to-Fortran mapping and `PORTING_NOTES.md` for
algorithm and validation details.

## Build

```text
fpm build
fpm test
```

Linear algebra delegates to the local `rfortran-linalg` dependency and its
pinned pure-Fortran LAPACK backend. System BLAS and LAPACK are not required.

## Example

```fortran
program statmod_example
use statmod
use r_mod, only: dp
implicit none
type(quad_rule_t) :: g
real(dp) :: m, v

print *, pinvgauss(1.0_dp, mean=2.0_dp, dispersion=0.5_dp)
g = gauss_quad(5, 'legendre')
print *, g%nodes
call expected_deviance(2.0_dp, 'poisson', m=m, v=v)
print *, m, v
end program
```

## License

Code derived from `statmod` retains the upstream `GPL-2 | GPL-3` choice,
represented here as `GPL-2.0-only OR GPL-3.0-only`. The supplied `r_mod.F90`
remains MIT-licensed. The vendored Tweedie translation retains
GPL-2.0-or-later. See `LICENSE`, `LICENSES/`, and `NOTICE.md`.
