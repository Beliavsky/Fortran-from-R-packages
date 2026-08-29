# tweedie-fortran

Modern Fortran/FPM translation of the computational core of the R package
`tweedie` 3.1.0.

The port preserves the package's Dunn-Smyth series and Fourier-inversion
algorithms, stored interpolation grids, Tweedie DPQR operations, deviance and
parameter conversions, likelihood calculations, and the numerical core of
profile fitting. R plotting, formula parsing, lifecycle aliases, and model
presentation code are intentionally omitted.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

BLAS and LAPACK must be available to the linker because the supplied `r_mod`
helper contains linear-algebra facilities. They are declared in `fpm.toml`.

## Main API

Use the umbrella module:

```fortran
use tweedie
```

Important procedures include:

- `dtweedie`, `ptweedie`, `qtweedie`, `rtweedie`
- `rtweedie_vec_params` for R-like recycling of vector `mu`/`phi` in RNG calls
- `dtweedie_series`, `ptweedie_series`
- `dtweedie_inversion`, `ptweedie_inversion`
- `dtweedie_saddle`
- `tweedie_dev`, `tweedie_lambda`, `tweedie_convert`
- `tweedie_integrand_values`
- `tweedie_loglik`, `tweedie_aic`
- `dtweedie_dlogfdphi`, `dtweedie_dldphi`
- `tweedie_phi_mle`
- `tweedie_glm_fit`, `tweedie_profile_grid`

`API_MAPPING.md` gives the mapping from the R API to the Fortran API.

## Numerical implementation

The upstream package already ships its Fourier inversion implementation in
Fortran. Those numerical files are retained unchanged and called directly by
the new standalone wrappers. R's runtime-specific C registration and printing
shims are not required.

The R-level Dunn-Smyth series calculations and the full set of ten stored
Chebyshev interpolation grids are translated into native Fortran modules.
All 4,160 stored grid coefficients were checked element-by-element against the
attached R source.

For p=1, p=2, and p=3, the public routines use the same Poisson, gamma, and
inverse-Gaussian special cases as the package. For 1<p<2, random generation is
direct compound Poisson-gamma; for p>2 it uses inverse-CDF generation, matching
upstream.

## Licensing

The upstream `tweedie` package declares **GPL (>= 2)**. Consequently,
tweedie-derived code in this port is **GPL-2.0-or-later**.

The separately supplied `r_mod.f90` remains **MIT-licensed**. Its build copy
is semantically identical to the supplied file and differs only by line wraps
needed for standard free-form source limits.

See `LICENSE`, `LICENSES/`, `NOTICE.md`, and the retained upstream metadata.

## References

The underlying numerical methods are due to Peter K. Dunn and Gordon K. Smyth:

- Dunn & Smyth (2005), *Statistics and Computing* 15(4), 267-280 — series
  evaluation of Tweedie exponential dispersion models.
- Dunn & Smyth (2008), *Statistics and Computing* 18(1), 73-86 — Fourier
  inversion for Tweedie exponential dispersion models.

The original package `CITATION` file is retained verbatim under `upstream/`.

## Validation status

A clean GNU Fortran 14.2.0 build under Fortran 2018 with
`-Werror=implicit-interface` and `-fcheck=all` passes the included tests.
Representative interpolation points in all ten stored-grid regions agree with
direct Fourier inversion to relative errors from roughly 3e-11 to 3e-8.

FPM itself was not installed in the validation environment, so `fpm test` was
not executed there; the same source graph was compiled and linked directly.
