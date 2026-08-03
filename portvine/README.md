# portvine-fortran

Modern Fortran/FPM translation of the computational core of the R package
`portvine`.

The package estimates rolling portfolio risk by fitting an ARMA-GARCH marginal
model to each return series, transforming standardized innovations to copula
uniforms, fitting a parametric C-vine or D-vine, simulating joint innovations,
and transforming them back to portfolio returns. D-vines can be sampled
conditionally on one or two variables for stress testing.

## Main features

- rolling moving-window ARMA-GARCH marginal fits;
- normal, skew-normal, Student-t, skew-t, GED, skew-GED, JSU, NIG, GHYP, and
  GH skew-t marginal transformations through the vendored `rugarch` port;
- greedy D-vine ordering using partial correlations of normal scores;
- parametric independence, Gaussian, Student-t, Clayton, Gumbel, Frank, Joe,
  BB1, BB6, BB7, BB8, and Tawn pair-copulas through the vendored vine library;
- D-vine and C-vine fitting, simulation, and rolling refits;
- conditional D-vine sampling for one or two leading variables;
- empirical VaR, mean ES, median ES, and Monte Carlo integrated ES;
- time-varying portfolio weights and realized portfolio returns;
- typed result objects retaining fitted marginal and vine models.

## Build and run

```text
fpm build
fpm test
fpm run
```

The package is self-contained: its two numerical dependencies are vendored as
FPM path dependencies.

GNU Fortran users can also run:

```text
scripts/build_all.sh checked
scripts/build_all.sh optimized
```

On Windows:

```bat
scripts\build_all.bat checked
scripts\build_all.bat optimized
```

## Matrix orientation

Fortran arrays use variables by observations:

```fortran
returns(n_assets, n_observations)
weights(n_assets, n_vine_windows)
```

This is the transpose of the usual R data-frame representation.

## Minimal setup

```fortran
use portvine
use rugarch, only : garch_spec, dist_sstd
use rvinecopulib, only : bicop_indep, bicop_gaussian, bicop_student

type(garch_spec) :: spec
type(marginal_settings_type) :: ms
type(vine_settings_type) :: vs
type(portvine_roll_result) :: answer

spec = make_portvine_spec(1, 1, 1, 1, cond_dist=dist_sstd)
ms = make_marginal_settings(750, 50, n_assets, spec)
vs = make_vine_settings(250, 25, vine_dvine, &
   [bicop_indep, bicop_gaussian, bicop_student])

call estimate_risk_roll(returns, ms, vs, [0.01_dp, 0.05_dp], &
   [risk_var, risk_es_mean], 2000, answer)
```

For conditional estimation, supply `cond_indices` and `cond_u`. The weights of
conditioning assets must be zero, matching the R package.

## Important adaptations

The attached `rugarch` translation does not jointly optimize all ARMA and GARCH
parameters in one likelihood. This package estimates ARMA coefficients by
iterated conditional least squares, then fits the GARCH model to those
innovations and filters the combined fixed specification. This is a practical
numerical counterpart, not a bit-for-bit reproduction of `rugarch::ugarchroll`.

The R package's unrestricted regular-vine mode is represented by a C-vine in
this native dependency set. D-vine mode, including conditional sampling, is the
closest direct counterpart and is the recommended mode for portvine workflows.
See `PORTING_NOTES.md` and `API_MAP.md`.

## License

The original R package is MIT licensed and its notice is preserved. The
combined Fortran distribution links GPL-3.0-only `rugarch` and vine libraries,
so the FPM project is distributed under GPL-3.0-only. See `NOTICE.md` and
`THIRD_PARTY_LICENSES.md`.
