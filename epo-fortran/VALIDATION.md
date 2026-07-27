# Validation

## Compiler

The release was validated with GNU Fortran 14.2.0.

Checked build flags:

```text
-std=f2018 -O0 -g -fcheck=all -fbacktrace
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
```

Optimized build flags:

```text
-std=f2018 -O2
-Wall -Wextra -Wconversion-extra -Wimplicit-interface -Werror
```

## Test programs

- `test_statistics`: sample means, covariance, correlation, and exact symmetry.
- `test_simple`: raw and normalized simple EPO at zero, partial, and full
  shrinkage.
- `test_anchored`: endogenous and exogenous anchored EPO, endogenous gamma,
  zero-shrinkage equivalence, and full-shrinkage anchor recovery.
- `test_validation`: invalid method, invalid shrinkage, invalid risk aversion,
  missing anchor, singular covariance, and zero-sum normalization.

Fixed portfolio references were independently calculated with NumPy using
64-bit floating-point linear algebra.

## Results

```text
test_anchored: PASS
test_simple: PASS
test_statistics: PASS
test_validation: PASS
```

The demo and both examples also compiled and ran in checked and optimized
configurations.

## Source audits

- `fpm.toml` parses as TOML.
- All translated text is ASCII.
- All free-form Fortran lines are at most 132 columns.
- Every Fortran source has `implicit none`.
- Every translated Fortran source contains an MIT SPDX identifier.
- Original and translated file SHA-256 manifests are included.
- The exact final ZIP is extracted and rebuilt before release.

## FPM availability

The FPM executable was not available in the validation container. Direct GNU
Fortran scripts therefore perform the checked and optimized builds. The
manifest uses ordinary automatic `src`, `app`, `example`, and `test` target
discovery.
