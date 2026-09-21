# Validation notes

The release manifest defines one library, two deterministic test programs, and one example. The intended FPM release commands are:

```text
fpm build
fpm test
fpm run --example dense_fpca
fpm clean --all
```

## Validation performed for this translation

On 2026-09-20 the maintained Fortran sources were compiled from a clean external build directory with GNU Fortran 14.2.0 using both an optimized Fortran 2018 build and a runtime-checking build. Both deterministic test programs passed, and the dense-FPCA example completed successfully.

The optimized validation used `-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror=line-truncation`. The runtime-checking validation used `-std=f2018 -O0 -g -fcheck=all -fbacktrace -Werror=line-truncation`.

The validation environment did not contain an `fpm` executable, and network restrictions prevented installing one. The literal `fpm build`, `fpm test`, example, and `fpm clean --all` commands were attempted and returned `fpm: command not found`. Therefore this file does not claim that FPM itself was executed successfully in that environment. The manifest was parsed as TOML, the two tests use separate FPM source directories to avoid multiple-main collection, and all source dependencies were validated by direct dependency-ordered compilation.

## Release audits

Before packaging, automated checks verify the following:

- every maintained Fortran dummy argument has `INTENT` or `VALUE`, is declared separately, and has a trailing FORD `!!` comment;
- maintained Fortran source is free form, ASCII, and uses the single exported `dp` kind;
- no code statement contains disallowed semicolon separation;
- no self-comparison NaN test is used;
- there are no system BLAS/LAPACK link directives or copied BLAS/LAPACK/ARPACK/R compatibility sources;
- maintained Fortran files have unique names and unique content;
- the README and `fpm.toml` agree on 48 of 48 mapped computational exports; and
- no object, module, executable, archive, cache, or ZIP artifact is included in the package tree.

Do not enable finite-math/fast-math assumptions because documented behavior includes IEEE NaN handling.
