# smoof-fortran

Modern Fortran/FPM translation of the computational benchmark functions in
`smoof` 1.7.0.

The public API is array-oriented rather than an emulation of R's closure and
attribute system. Scalar objectives are pure functions such as `ackley(x)`,
`rosenbrock(x)`, and `hartmann(x)`. Multi-objective benchmarks are subroutines
such as `dtlz2(x, m, f)`, `wfg1(z, m, k, f)`, and `uf1(x, f)`.

## Modules

- `smoof_single`: ordinary single-objective benchmark functions.
- `smoof_multi`: DTLZ1-7, ZDT1/2/3/4/6, WFG1-9, MOP1-7, BK1,
  Viennet, Kursawe, Dent, and bi-sphere.
- `smoof_cec09`: CEC 2009 UF1-UF10.
- `smoof_cec2019`: SYMPART, OMNI, and MMF1-MMF15a families.
- `smoof_ed`: ED1 and ED2.
- `smoof_nk`: array-level NK-landscape evaluator.
- `smoof`: convenience module re-exporting all of the above.

## Build

```text
fpm build
fpm test
fpm run --example single_example
fpm run --example multi_example
```

The project is written to standard free-form Fortran 2018 and does not require
BLAS, LAPACK, R, Rcpp, or Armadillo.

See `TRANSLATION_COVERAGE.md` for exact coverage and deliberate omissions.
