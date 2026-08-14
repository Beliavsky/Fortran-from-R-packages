# Validation

The release was compiled from a clean build directory with GNU Fortran 14.2.0
using:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

and linked with system LAPACK/BLAS.

Permanent tests:

1. `test_irlba.f90`
   - dense IRLB vs full LAPACK SVD;
   - implicit centering/scaling vs explicitly transformed matrix;
   - native CSC sparse IRLB vs dense reference;
   - restart from an existing SVD result;
   - smallest-singular-value fallback;
   - rank-deficient matrix with exact zero tail singular values.
2. `test_randomized.f90`
   - 60 randomized dense problems;
   - centered/scaled variants for every case;
   - sparse variants for every third case;
   - maximum allowed relative leading-singular-value error: 2e-7.
3. `test_highlevel.f90`
   - truncated PCA;
   - symmetric partial eigenvalue recovery with positive/negative spectrum;
   - sparse SVD target support;
   - dense and CSC sparse randomized SVD.
4. `test_complex.f90`
   - complex dense SVD residual check.
5. `example/basic.f90`
   - example truncated dense SVD.

All tests pass under the flags above. Source-form audit also verifies no `.f90`
line exceeds 132 characters, source/docs are ASCII-clean, and `fpm.toml` parses
with Python's TOML parser.
