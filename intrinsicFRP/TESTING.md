# Testing

## FPM

```text
fpm test
fpm run
fpm run --example basic_risk_premia
```

## Direct GNU Fortran validation

On Unix-like systems:

```text
./tools/test_gnu.sh debug
./tools/test_gnu.sh release
```

On Windows with GNU Fortran:

```text
tools\test_gnu.bat
```

The debug script uses Fortran 2018 conformance checks, warnings as errors,
bounds checking, runtime checks, and traps for invalid, divide-by-zero, and
overflow exceptions. `-Wno-maybe-uninitialized` suppresses a known GNU Fortran
false-positive class involving allocatable assignment; all other enabled
warnings remain errors.

The tests cover:

1. TFRP, FM/KRS FRP, FM/GKR SDF, GKR screening, HJ distance, and HAC covariance.
2. KP-style and Chen-Fang identification tests and Giglio-Xiu PCA premia.
3. Oracle GCV/CV/rolling tuning, soft thresholding, and FGX testing.
4. Invalid dimensions, nonfinite input, singular cases, and invalid controls.
