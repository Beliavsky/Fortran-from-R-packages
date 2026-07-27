# Validation

## Environment

The release was validated with GNU Fortran 14.2.0.

Checked compilation flags:

```text
-std=f2018
-Wall
-Wextra
-Wconversion-extra
-Werror
-Wno-compare-reals
-fcheck=all
-fbacktrace
```

The explicit `-Wno-compare-reals` suppresses warnings for intentional
floating-point boundary and tie checks. Other enabled warnings remain errors.

## Test programs

```text
test_data: PASS
test_microstructure: PASS
test_models_jumps: PASS
test_realized: PASS
```

Coverage includes:

- log-return construction;
- fixed-time aggregation;
- previous-tick and refresh-time synchronization;
- trade direction, liquidity measures, cleanup masks, and OHLCV;
- realized covariance and higher moments;
- bipower, tripower, quadpower, minimum, median, pre-averaged, average, and
  threshold estimators;
- PSD preservation;
- Hayashi-Yoshida and two-scale calculations;
- ReMeDI output;
- lead-lag contrasts;
- kernel spot volatility and drift;
- HAR fitting and forecasting;
- HEAVY fitting, stationarity, and forecasting;
- BNS jump-test and integrated-variance inference output validity.

The deterministic realized-measure checks were independently calculated with
NumPy/Python formulas.

## Executables

The following also compile and run:

```text
highfrequency_demo
realized_measures
microstructure
```

## Optimized build

The complete suite is also compiled and run with `-O2`.

## FPM

The validation environment did not contain the FPM executable. The project uses
standard automatic FPM source, test, application, and example discovery. The
manifest is valid TOML and contains no compiler-specific build directives.

## Audits

The release validation also checks:

- all translated Fortran files contain an SPDX identifier;
- all program units use `implicit none`;
- translated source is ASCII;
- no translated free-form source line exceeds 132 characters;
- the original package is retained;
- original and translated file SHA-256 manifests match;
- the final ZIP can be extracted and rebuilt from a clean directory.
