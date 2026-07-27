# Validation

## Toolchain

The project was compiled with GNU Fortran 14.2.0 in Fortran 2018 mode.

Checked build flags:

```text
-O0 -g -std=f2018 -Wall -Wextra -Werror -Wconversion-extra
-fcheck=all -fbacktrace
```

Optimized build flags:

```text
-O2 -std=f2018 -Wall -Wextra -Werror -Wconversion-extra
```

## Test suites

- `test_stats`: fixed Sharpe, modified-Sharpe, alpha, and factor-loading references
- `test_inference`: alpha, HAC, asymptotic Sharpe, modified-Sharpe, and deterministic bootstrap tests
- `test_peer_ratios`: null-proportion adjustment, lambda selection, ratio ranges, and positive/negative split
- `test_screening`: within-group, cross-group, target, rolling, coefficient-level, and bootstrap screening

Expected output:

```text
test_inference: PASS
test_peer_ratios: PASS
test_screening: PASS
test_stats: PASS
validation: PASS
```

The demo and both examples are compiled and executed by the validation script.

## Portability checks

- checked and optimized builds pass
- warnings are promoted to errors
- bounds and allocation checks are enabled in the checked build
- source is free-form Fortran and remains within 132 columns
- every Fortran unit uses `implicit none`
- every translated Fortran file carries the preserved SPDX identifier and attribution
- the FPM manifest parses as TOML and uses automatic source, app, example, and test discovery

## FPM

The FPM executable was not installed in the validation environment. The source
was therefore compiled directly with the included GNU Fortran scripts. The
project layout and manifest are intended for:

```text
fpm build
fpm test
```
