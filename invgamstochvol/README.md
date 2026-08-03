# invgamstochvol-fortran

Modern Fortran translation of the computational code in the R package
`invgamstochvol` 1.0.0.

The library evaluates the closed-form marginalized likelihood of the stationary
inverse-gamma stochastic-volatility model described by Leon-Gonzalez and Majoni
(2023), and draws inverse volatilities from the exact smoothing distribution.

## Implemented API

- `ourgeo` - truncated Gauss hypergeometric function `2F1`
- `lik_clo` - exact closed-form log likelihood and smoothing recursion tables
- `draw_k0` - exact posterior draw of the inverse-volatility path

The likelihood result retains the coefficient and rising-factorial tables used
by the original package, so a likelihood evaluation can be passed directly to
`draw_k0` without recomputation.

## Build with FPM

```text
fpm build
fpm test
fpm run --example likelihood_example
fpm run --example posterior_example
fpm run demo_invgamstochvol
```

The project has no external dependencies.

## Minimal example

```fortran
use invgamstochvol

type(invgam_likelihood_result) :: fit
real(dp) :: residuals(8)
real(dp), allocatable :: inverse_volatility(:)
integer :: status

residuals = [0.20_dp, -0.10_dp, 0.35_dp, -0.25_dp, &
   0.05_dp, 0.40_dp, -0.30_dp, 0.15_dp]

call lik_clo(residuals, 0.7_dp, 4.1_dp, 0.85_dp, fit)
call draw_k0(fit, 4.1_dp, 0.85_dp, 0.7_dp, inverse_volatility, &
   status=status)
```

## Porting notes

R matrices, Rcpp lists, Armadillo containers, and OpenMP controls are replaced
with standard Fortran arrays and derived types. `nproc` and `nproc2` are
accepted for source-level familiarity but the current implementation is
serial. Posterior simulation uses a deterministic portable random-number
stream when a seed is supplied.

See `API.md`, `PORTING.md`, and `TRANSLATION_COVERAGE.md` for details.

## License

MIT. The original package sources and metadata are retained under `original/`.
