# greeks-fortran

A modern Fortran translation of the computational code in the R package
`greeks` 1.5.6. The project is an FPM library and has no external numerical
or plotting dependencies.

## Included functionality

- Exact Black-Scholes Greeks for European vanilla, cash-or-nothing, and
  asset-or-nothing options.
- Exact Black-Scholes Greeks for continuously monitored geometric Asian calls
  and puts.
- American call and put values and Greeks from a corrected CRR binomial tree.
- Malliavin Monte Carlo Greeks for European, arithmetic Asian, and geometric
  Asian options.
- A geometric-control-variate estimator for arithmetic Asian options.
- Black-Scholes Halley implied volatility and a general Newton implied
  volatility interface.
- Black-Scholes and jump-diffusion path simulation for Asian estimators.
- Brownian-path, trapezoidal-integration, and path-functional helpers.
- Standard errors for Monte Carlo estimates.
- Optional pure Fortran payoff callbacks for the general Malliavin routines.

The Shiny, ggplot2, plotly, tibble, and R S3/UI layers are not translated.

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example european_greeks
fpm run --example american_and_implied
fpm run --example asian_monte_carlo
```

The package version in `fpm.toml` is the plain semantic version `1.5.6`.

## Minimal example

```fortran
program example
  use greeks
  implicit none
  type(greek_result) :: result
  character(len=24) :: requested(3)

  requested = [character(len=24) :: 'fair_value', 'delta', 'gamma']
  call bs_european_greeks(100.0_dp, 100.0_dp, 0.05_dp, 1.0_dp, &
    0.20_dp, 0.0_dp, 'call', requested, result)
  print *, result%values
end program example
```

## API conventions

R dots are changed to underscores, and public R names are represented by
Fortran procedures such as `bs_european_greeks`,
`binomial_american_greeks`, and `implied_volatility`. The R dispatcher
`Greeks` is represented by `option_greeks`, because a procedure named
`greeks` would collide with the public module name.

Inputs are scalar. R's single-vectorized-parameter behavior is represented by
calling the Fortran procedure in a loop. Requested Greek names are passed as a
character array. Results use the `greek_result` derived type, which contains
names, values, Monte Carlo standard errors, a status code, and a message.

See `API.md`, `PORTING.md`, and `TRANSLATION_COVERAGE.md` for details.

## License

MIT, matching the original package. The original license, metadata, R sources,
and C++ computational sources are retained under `original/`.
