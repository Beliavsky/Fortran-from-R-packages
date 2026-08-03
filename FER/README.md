# FER Fortran

A modern Fortran translation of **FER 0.94** (Financial Engineering in R).
It provides self-contained implementations of standard option-pricing models
without external numerical-library dependencies.

## Features

- Bachelier price and implied volatility.
- Black-Scholes/Black-76 price and implied volatility.
- CEV option price and probability mass at zero.
- Hagan (2002) SABR implied volatility and price.
- Choi (2019) NSVh lambda-1 price.
- Margrabe exchange option formula.
- Kirk, Bjerksund-Stensland (2014), and Bachelier spread formulas.
- Fortran-native special functions and bracketed implied-volatility solvers.
- Demonstration, example, and regression-test programs.

## Directory layout

```text
app/       default demonstration program
example/   additional FPM example
original/  retained original R package source and documentation
src/       Fortran library modules
test/      regression tests
```

## Build with FPM

```text
fpm build
fpm test
fpm run
fpm run --example vanilla_options
```

The package uses standard Fortran 2018.

## License

The original package declares `GPL (>= 2)`. This translation is licensed under
**GPL-2.0-or-later**. GPL-2.0 and GPL-3.0 license texts are included, and the
original package files are retained under `original/` for provenance.
