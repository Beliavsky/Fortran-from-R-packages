# Build validation

Validation performed in the translation environment:

- GNU Fortran 14.2.0, `-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all`.
- All library modules compiled from a clean directory.
- `test/test_tmb.f90` compiled, linked, and ran successfully: `All TMB Fortran tests passed.`
- `example/simple.f90` compiled, linked, and ran successfully.
- `audit.py` passed after cleanup.
- No BLAS, LAPACK, ARPACK, CHOLMOD, Eigen, CppAD, R, or other external library was linked.

The execution environment did not provide the Fortran Package Manager (`fpm`) binary, and its
network sandbox could not resolve external hosts to install it. Consequently, the requested literal
`fpm build`, `fpm test`, and `fpm clean --all` commands could not be executed here. The included
`fpm.toml` follows the current manifest specification, and the same source/test/example targets were
compiled directly with gfortran as the available verification route.
