# Validation

This translation was validated on 2026-09-22 with GNU Fortran 14.2.0.

The sandbox does not contain the `fpm` executable and cannot perform Git network fetches. Therefore the literal target commands `fpm build`, `fpm test`, `fpm run --example basic`, and `fpm clean --all` could not be completed here. They must be rerun on the target Windows/FPM installation.

For local numerical validation only, an out-of-tree adapter exposing the same six `la_lapack` generic names used by this package (`gees`, `geev`, `gges`, `ggev`, `trsen`, `tgsen`) was linked to the sandbox's system Reference LAPACK. That adapter and the system link flags are **not included** in the package and are not present in `fpm.toml`.

Validation performed on the maintained package sources:

- strict compilation with `-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace`;
- deterministic standard/generalized real/complex eigenvector residual tests;
- real/complex Schur and generalized-Schur reconstruction tests;
- real/complex TRSEN/TGSEN and high-level `ordqz` ordering tests;
- high-level `qz_decompose`, `qz_geigen`, `geigen`, `fda_geigen`, and Hermitian-transpose tests;
- basic example execution;
- independent optimized rebuild with `-O2` before packaging;
- static release audit for free-form/ASCII source, 132-column limit, dummy `INTENT` and FORD comments, semicolon-packed statements, self-comparison NaN tests, duplicate Fortran files, copied dependencies, unsafe math flags, and build artifacts.

The final FPM/Windows validation should run:

```text
fpm clean --all
fpm build
fpm test
fpm run --example basic
fpm clean --all
```
