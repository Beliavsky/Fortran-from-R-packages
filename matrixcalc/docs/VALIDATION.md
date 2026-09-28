# Validation

This checkpoint was validated on 2026-09-22 with GNU Fortran 14.2.0. The maintained package sources were compiled against temporary out-of-tree interface implementations of the sibling `rfortran-core` and `rfortran-linalg` APIs. Those validation files are not part of this package or its ZIP.

## Compiler validation

Two independent builds passed:

- strict checked build: `-std=f2018 -O0 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -fbacktrace -ffpe-summary=all`;
- optimized build: `-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface`.

In both builds all deterministic tests passed:

- `test_core`
- `test_structured`
- `test_special`
- `test_linalg`

The `basic` example also compiled and ran successfully in both builds.

The temporary validation implementation of `rfortran-linalg` was linked outside the package for test execution. The package itself contains no BLAS/LAPACK source and its FPM manifest contains no system BLAS/LAPACK link directives.

## Static release audit

The final audit reported zero findings and checked:

- 60 mapping entries for 60 exported computational R functions;
- agreement of counts, fraction, mappings, and the empty untranslated list;
- existence of every mapped upstream R source and Fortran source;
- presence of every mapped public Fortran procedure/generic;
- 10 maintained Fortran source/test/example files;
- ASCII-only maintained Fortran and a 132-column maximum;
- explicit `INTENT` or `VALUE` and trailing FORD `!!` documentation for dummy declarations;
- no duplicate Fortran source files;
- no vendored dependency sources;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no semicolon-packed statements;
- no self-comparison NaN tests;
- no unsafe fast-math flags;
- no object/module/executable/library files, caches, or nested ZIP files in the release tree.

## FPM availability

The required literal FPM commands were attempted in this environment:

```text
fpm build
fpm test
fpm run --example basic
fpm clean --all
```

Each returned exit code 127 with `fpm: command not found` because the sandbox does not provide the FPM executable. Therefore this document does not claim a literal FPM run here. The package is laid out for those commands when extracted as a top-level repository directory alongside `rfortran-core` and `rfortran-linalg`.
