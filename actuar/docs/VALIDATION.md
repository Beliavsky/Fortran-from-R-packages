# Validation

`actuar-fortran` v0.3.0 was rebuilt from an empty object/module directory with
GNU Fortran 14.2.0 using:

```text
-std=f2018 -Werror=implicit-interface -Werror=trampolines -fcheck=all -O0
```

The complete actuar source tree and vendored `expint-fortran` dependency were
compiled together.

Actuar tests:

```text
test_core: PASS
test_risk: PASS
test_v02: PASS
test_v03: PASS
```

Vendored expint tests:

```text
test_expint: PASS
test_gammainc: PASS
```

All three actuar examples compile and run. The v0.3 example reports:

```text
MDE exponential rate:   2.000000
Covered-loss CDF at 0.6:   0.958490
Covered-loss density at 0.6:   0.094341
Probability mass at payment limit:   0.004277
```

Representative v0.3 regression checks include:

- individual and grouped Cramer-von Mises minimum-distance fitting;
- modified chi-square recovery of an exponential rate from exact grouped
  probabilities;
- grouped layer-average-severity optimization;
- ordinary and franchise deductible coverage branches;
- per-loss and per-payment conditioning;
- zero-payment and policy-limit probability masses;
- coverage density Jacobian under coinsurance and inflation;
- exact zero-noise recovery of Hachemeister barycenter contract regressions;
- one-level `hierarc_exact_fit` agreement with iterative Buhlmann-Straub for
  process variance, between variance, collective mean and credibility factors;
- callback-driven hierarchical simulation of frequency/severity mixing,
  frequencies, claim allocation and terminal aggregates.

All `.f90` files in the package and vendored dependency are at most 132
characters per source line. `fpm.toml` parses with Python `tomllib` and reports
version `0.3.0` with the local `expint-fortran` dependency.

FPM itself was not installed in the validation environment, so the standard
FPM layout was validated by direct compilation rather than invoking FPM.
