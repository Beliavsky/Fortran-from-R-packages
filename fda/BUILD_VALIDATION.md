# Build validation

Validation environment:

- GNU Fortran 14.2.0
- Fortran 2018 mode (`-std=f2018 -pedantic`)
- runtime checking (`-fcheck=all -fbacktrace`)
- warnings (`-Wall -Wextra -Wimplicit-interface`) with implicit-interface diagnostics promoted to errors

## Direct strict validation performed

The maintained source was compiled from a clean directory against temporary,
non-packaged validation modules exposing the current `r_kinds` and `r_linalg`
interfaces used by `fda`.

The deterministic test program `test/test_fda.f90` then ran successfully and
printed:

```text
all fda tests passed
```

The example `example/smooth_sine.f90` also compiled and ran successfully,
printing:

```text
effective df:    7.73008
SSE:            1.6243E-04
```

The tests cover deterministic reference values/identities for:

- B-spline values and partition of unity;
- normalized Fourier values, derivatives, and Gram matrix;
- monomial, exponential, power, and polygonal basis derivatives/values;
- roughness matrices;
- trapezoidal and Simpson quadrature;
- polynomial interpolation;
- symmetric solves and generalized eigenvalues;
- functional-data evaluation and inner products;
- exact unpenalized polynomial smoothing and effective-df/lambda inversion;
- Fourier-basis FPCA variance proportions and self-CCA correlations;
- adaptive solution of `u'' + u = 0` through `pi/2`.

## FPM/fprettify environment limitation

`fpm` was not installed in the execution image.  An attempt to retrieve the
current FPM 0.13.0 Linux release binary failed.  `fprettify` was also absent;
an attempted pip installation failed because outbound DNS/network access is
unavailable in the environment.

Therefore the literal commands `fpm build`, `fpm test`, `fpm clean --all`, and
`fprettify` could not truthfully be executed here.  The package contains a
normal FPM manifest with sibling dependency paths, and the maintained source
was instead compiled and executed directly with strict GNU Fortran checking.
No claim is made that unavailable commands were run.
