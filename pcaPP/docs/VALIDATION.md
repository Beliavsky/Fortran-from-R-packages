# Validation

Validation was performed with GNU Fortran 14.2.0 on the translation sources.
No `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or system BLAS/LAPACK flags
were used.

## Direct compiler validation

The package was compiled twice outside the source tree so no object/module
files were retained in the deliverable.

Strict checked build:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -fcheck=all -fbacktrace -O0
```

Optimized warning-clean build:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -O2
```

In both configurations:

- `test/unit/main.f90` passed (`All pcaPP unit tests passed.`);
- `test/algorithms/main.f90` passed (`All pcaPP algorithm tests passed.`);
- `example/robust_pca.f90` ran successfully.

The deterministic tests cover Qn finite-sample behavior, Kendall tau-b with
ties, scaling, all seven exported spatial-median paths, PCAgrid, PCAproj,
sPCAgrid, covariance reconstruction, TPO/BIC tuning-grid selection, and a
seeded Zou-simulation shape/finiteness check.

## Static release audit

The final source-tree audit verifies:

- 19 of 19 exported computational R functions have exactly one mapping entry;
- all mapped Fortran procedures exist and are public in their declared module;
- all `r_source` and `fortran_source` coverage paths exist;
- README and `fpm.toml` both report 19/19 (100.0%) mapped;
- one package-local `dp = real64` definition is used consistently;
- all dummy arguments are declared individually with `INTENT` or `VALUE` and
  a trailing `!!` FORD documentation comment;
- maintained Fortran is ASCII, free-form, and no line exceeds 132 characters;
- no code semicolons, legacy `double precision`/`real*8`/D-exponent literals,
  or self-comparison NaN tests are present;
- no duplicate Fortran source files are present;
- no vendored numerical dependency source or system BLAS/LAPACK link is
  present; and
- no build products, module files, executables, caches, or nested ZIP archives
  are present in the release tree.

## FPM environment limitation

The requested literal FPM commands were attempted from the package root:

```text
fpm build
fpm test
fpm run --example robust_pca
fpm clean --all
```

All four returned exit status 127 with `fpm: command not found`, because this
execution environment does not provide an FPM executable.  These attempts are
not reported as successful FPM runs.  Manual cleanup was performed after the
failed `fpm clean --all` attempt, and the same source/test/example targets were
validated directly with GNU Fortran as described above.
