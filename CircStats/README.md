# CircStats-fortran

Modern Fortran 2018 computational port of the R package **CircStats 0.2-7**.
The project uses the Fortran Package Manager (FPM) layout and is self-contained.

## Scope

Translated numerical functionality includes:

- circular mean, resultant length, dispersion, range, summaries, correlations and trigonometric moments;
- concentration estimation (`A1`, `A1inv`, von Mises MLE and bias correction);
- circular change-point statistics;
- circular-on-circular trigonometric regression and higher-order tests;
- cardioid, triangular, von Mises, mixed von Mises, wrapped Cauchy and wrapped normal densities;
- von Mises cumulative probabilities;
- RNGs for cardioid, triangular, von Mises, mixed von Mises, wrapped Cauchy, wrapped normal and Levy alpha-stable laws;
- Rayleigh, V0, Kuiper, Watson uniform, Watson von Mises and Watson two-sample tests;
- Rao spacing and Rao homogeneity tests;
- Kent-Tyler wrapped-Cauchy maximum likelihood estimation;
- von Mises bootstrap confidence intervals;
- degree/radian conversion and combinatorial helpers.

The package's 43 x 4 Rao spacing critical-value table is embedded exactly from
`data/rao.table.rda`.

## Intentionally omitted

Graphics-only behavior is not translated: `circ.plot`, `rose.diag`, `plotedf`, and
the plotting side effects of `pp.plot`/`watson.two`. The numerical parameter estimates
returned by `pp.plot` are available as `pp_fit`.

## API notes

R names containing dots are mapped to Fortran-style underscore names. Examples:

- `circ.mean` -> `circ_mean`
- `est.rho` -> `est_rho`
- `r.test` -> `rayleigh_test`
- `v0.test` -> `v0_test`
- `wrpcauchy.ml` -> `wrpcauchy_ml`
- `vm.bootstrap.ci` -> `vm_bootstrap_ci`
- `I.0`, `I.1`, `I.p` -> `i0`, `i1`, `ip`

R data-frame/list returns are represented by derived result types.

## Build

```sh
fpm build
fpm test
fpm run --example basic
```

The code has also been validated directly with GNU Fortran using:

```sh
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

## License

CircStats declares `License: GPL-2`. This translation is distributed under
GPL-2.0-only, preserving the upstream license. See `COPYING`, `LICENSES.md`,
`UPSTREAM.md`, and the complete source snapshot in `upstream/`.
