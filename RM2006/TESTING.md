# Testing

The test target `test/test_rm2006.f90` covers:

- A complete two-asset, five-observation reference calculation.
- All six covariance slices against independently generated fixed values.
- Time-scale and normalized-weight reference values.
- The constant-return identity, for which every covariance slice must equal
  the same outer product.
- Samples shorter than the default `kmax=14`.
- Invalid parameter handling.
- Non-finite input handling.

The fixed reference matrices were generated independently from the published R
algorithm, not by calling the Fortran implementation.

## Strict compiler validation

The package was tested with GNU Fortran 14.2 using:

```text
-std=f2018 -pedantic -Wall -Wextra -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow
-fimplicit-none -O2
```

Both examples and the test program compiled and ran successfully under those
settings.

The validation environment did not provide the `fpm` executable. The FPM
manifest was parsed independently, and all source, application, example, and
test targets were compiled directly with `gfortran` using the dependency order
specified by the modules.
