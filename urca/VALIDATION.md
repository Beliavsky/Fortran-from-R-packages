# Validation

## Compiler configuration

Validated with GNU Fortran 14.2.0, BLAS, and LAPACK using:

```text
-std=f2018 -fcheck=all -Werror=implicit-interface -O0
```

No unlimited-line-length extension is required.

FPM was not installed in the validation environment, so the same source
dependency graph was compiled directly with `gfortran`. `fpm.toml` is
included for normal FPM builds.

## Regression suite

Six test programs are included:

- `test_univariate.f90`
- `test_cointegration.f90`
- `test_restrictions.f90`
- `test_breaks.f90`
- `test_mackinnon.f90`
- `test_smoke.f90`

The example program is also compiled and executed.

All six tests and the example pass under runtime bounds/checking.

## Deterministic reference dataset

`test/reference_data.txt` contains 240 observations generated from a fixed
NumPy random seed. The first series is stationary; the remaining three form a
cointegrated synthetic system. Expected values were reconstructed
independently with NumPy/SciPy matrix calculations rather than copied from the
Fortran output.

Representative checks include:

- ADF drift/fixed-lag tau: `-4.6794330253858885`
- ADF phi1: `10.997867559051583`
- KPSS mu: `0.0922438642502086`
- Phillips-Perron Z-tau: `-6.311792252486218`
- ERS constant DF-GLS: `-4.220338063518298`
- ERS trend DF-GLS: `-4.277327382260179`
- Zivot-Andrews both-break statistic: `-5.25698543661746`
- Zivot-Andrews selected break: `105`
- Phillips-Ouliaris Pu: `109.9797706037426`
- Johansen eigenvalues approximately
  `0.3155669409245894, 0.2180893796195130, 0.0078980956581318`
- Johansen trace statistics approximately
  `1.887209295673415, 60.43874149122232, 150.67987687768203`
- `blrtest`: `23.203235948811248`
- `alrtest`: `15.67438347881805`
- `ablrtest`: `23.461224117377938`
- `bh5lrtest`: `56.64603846973627`
- `bh6lrtest`: `56.63732560007072`
- `lttest`: `2.9542314288676184`
- `cajolst` break point: `192`

The small differences tolerated in a few tests reflect BLAS/LAPACK rounding
and alternative but mathematically equivalent matrix factorizations.

## MacKinnon checks

The translated response-surface evaluator was independently reproduced in
Python from the same upstream coefficient data. Quantiles and p-value
inversion agree to tight floating-point tolerances across deterministic
specifications, finite sample sizes, and asymptotic cases.

Examples of asymptotic 5% tau quantiles:

- no constant: `-1.9408467713813484`
- constant: `-2.8613704741924773`
- constant + trend: `-3.4098441424441335`
- constant + trend + quadratic trend: `-3.8320053632822271`

`test_mackinnon.f90` also verifies finite-sample and normalized-statistic
cases and checks that quantiles invert to approximately 0.05.
