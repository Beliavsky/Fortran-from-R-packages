# derivmkts-fortran

A self-contained modern Fortran translation of the computational routines in
Robert L. McDonald's R package `derivmkts` 0.2.5.1.

The project uses the Fortran Package Manager (FPM), has no external numerical
dependencies, and preserves the original MIT license and attribution.

## Build and test

```text
fpm build
fpm test
fpm run derivmkts_demo
```

Examples:

```text
fpm run --example analytic_options
fpm run --example monte_carlo
fpm run --example simulation_and_tree
```

## Main functionality

- Black-Scholes calls, puts, asset binaries, and cash binaries
- Implied volatility and implied spot calculations
- Analytical Black-Scholes Greeks and callback-based numerical Greeks
- Bond price, yield, duration, and convexity
- Geometric-average Asian price and strike options
- Arithmetic and geometric Asian Monte Carlo
- Geometric-Asian control variates for arithmetic Asian calls
- Up/down and in/out barrier binaries and standard barrier options
- Immediate and deferred barrier rebates
- Perpetual American calls and puts
- Geske compound calls and puts
- Merton lognormal jump-diffusion prices
- European and American binomial trees with hedge trees and Greeks
- Single- and multi-asset lognormal path simulation with optional jumps
- Binomial distribution and quincunx simulation data

## Basic example

```fortran
program example
    use derivmkts, only: dp, bscall, bscallimpvol
    implicit none
    real(dp) :: price, implied_vol

    price = bscall(40.0_dp, 40.0_dp, 0.30_dp, 0.08_dp, 0.25_dp, 0.0_dp)
    implied_vol = bscallimpvol(40.0_dp, 40.0_dp, 0.08_dp, 0.25_dp, &
        0.0_dp, price)

    print '(a,f12.6)', 'price:       ', price
    print '(a,f12.6)', 'implied vol: ', implied_vol
end program example
```

## Typed results

R lists and data frames are represented by Fortran derived types:

- `option_pair`
- `greek_result`
- `perpetual_result`
- `compound_result`
- `binomial_result`
- `asian_mc_result`
- `simulation_result`
- `quincunx_result`

For example, `binomopt(..., returntrees=.true.)` returns stock, option,
probability, exercise, delta, and bond trees in a `binomial_result`.

## Numerical conventions

All real calculations use:

```fortran
integer, parameter :: dp = kind(1.0d0)
```

Rates and dividend yields are continuously compounded unless the bond routines
explicitly use periodic compounding. Volatility is annualized. Monte Carlo
standard-deviation fields reproduce the upstream convention: they are the
standard deviations of discounted payoff observations, not standard errors of
the estimated means.

## Scope

All exported numerical routines are translated. Plot-producing behavior from
`binomplot` and the animated display in `quincunx` are excluded. Their complete
underlying tree, probability, path, count, and expected-frequency data are
available to the caller.

R argument recycling, language introspection, list/data-frame formatting, and
random-seed restoration are replaced by explicit typed Fortran interfaces and
optional deterministic seeds.

See [COVERAGE.md](COVERAGE.md), [PORTING_NOTES.md](PORTING_NOTES.md), and
[VALIDATION.md](VALIDATION.md) for details.

## License

MIT. See [LICENSE](LICENSE) and [NOTICE](NOTICE). The original package is
retained in `original/derivmkts-0.2.5.1`.
