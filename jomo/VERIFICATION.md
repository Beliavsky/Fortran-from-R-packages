# Verification

The maintained Fortran source was verified with GNU Fortran 14.2.0 using
strict interface diagnostics and runtime checks. All 14 deterministic test
programs passed and both example programs completed successfully. The expanded
suite includes the full SMC proposal/rebuild driver, all four single-level
substantive families, multilevel binary/ordinal random effects, heterogeneous
level-1 covariance proposals, and level-2 continuous/categorical proposals.

Compiler verification flags included:

```text
-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all -fbacktrace
```

The final source audit checks:

- free-form line lengths;
- duplicate Fortran source content;
- multiple dummy arguments on one declaration;
- missing `INTENT`/`VALUE` attributes;
- missing trailing FORD `!!` dummy-argument documentation;
- semicolon-separated statements;
- self-comparison NaN idioms;
- `double precision`, `real*8`, `kind(0.0d0)`, and D-exponent literals;
- copied BLAS/LAPACK/ARPACK or translated dependency source;
- build products and caches.

All of those source/package audits passed across 29 maintained/test/example Fortran files.

## FPM availability in the translation environment

The execution environment used for this translation did not contain an `fpm`
executable.  The required commands were nevertheless invoked explicitly:

```text
fpm build
fpm test
fpm clean --all
```

Each invocation returned `command not found` (exit status 127).  An attempt to
obtain the current upstream FPM binary was blocked by the environment's binary
network-download path.  No substitute or fake FPM executable was used.

The package retains a normal `fpm.toml` layout (`src`, `test`, and `example`),
and the same complete source/test/example tree was compiled directly with
GNU Fortran as described above.  Before packaging, any temporary compiler
products were kept outside the package directory and the package itself was
checked to contain no build artifacts.
