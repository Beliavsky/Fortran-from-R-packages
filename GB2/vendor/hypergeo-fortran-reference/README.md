# hypergeo-fortran

Modern free-format Fortran translation of the computational code in the R package
`hypergeo` 1.2-14 by Robin K. S. Hankin.

The library evaluates the Gauss hypergeometric function 2F1 for real or complex
parameters and arguments and exposes the continuation algorithms that the R package
makes available individually.  It also implements generalized hypergeometric series.

## Highlights

- complex Lanczos gamma, log-gamma, factorial, and digamma support
- generalized hypergeometric power series and continued fractions
- Gauss 2F1 power-series evaluation and analytic continuation
- Abramowitz and Stegun 15.3.3--15.3.14 transformations
- the integer-parameter limiting formulas used by `hypergeo_cover1/2/3`
- Wolfram special formulas used by `hypergeo_cover3`
- Gosper evaluation near the two difficult critical points
- Buhring continuation
- Euler integral evaluation with endpoint-safe tanh-sinh quadrature
- Cauchy/residue evaluation through the translated `elliptic` package
- hypergeometric ODE continuation through the translated `deSolve` package
- generalized continued fractions through the translated `contfrac` package

All Fortran source in this package and its vendored Fortran dependencies is free format
`.f90`; there are no fixed-form `.f`, `.for`, or `.f77` files.

## Build

```text
fpm build
fpm test
fpm run --example hypergeo_demo
```

The dependencies are vendored under `dependencies/`, so no network access is required.

## Basic use

```fortran
use hypergeo_fortran, only : dp, hypergeo

complex(dp) :: value

value = hypergeo(cmplx(1.21_dp, 0.0_dp, dp), &
                 cmplx(1.443_dp, 0.0_dp, dp), &
                 cmplx(1.88_dp, 0.0_dp, dp), &
                 cmplx(2.13_dp, 0.68_dp, dp))
```

For the example above the result is approximately

```text
-0.713344984340046 + 0.596474790855054 i
```

`hypergeo_info` may be supplied to the scalar `hypergeo` call to inspect which
continuation route succeeded.

## API naming

R function names containing dots are mapped to valid Fortran identifiers by replacing
`.` with `_`.  For example, `f15.3.11()` becomes `f15_3_11()` and
`w07.23.06.0029.01()` becomes `w07_23_06_0029_01()`.

See `TRANSLATION_NOTES.md` for the detailed coverage map and numerical differences.

## License

The upstream package is GPL-2.  This translation is distributed under GPL-2.0-only.
The vendored dependencies retain their own GPL-compatible licenses and attribution.
