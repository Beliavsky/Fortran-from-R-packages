# Build validation

Validation was performed on 2026-09-02 with GNU Fortran 14.2.0.

## Dependency/API validation

Before implementation, the `Fortran-from-R-packages` repository was checked for
reusable translations. This package uses sibling path dependencies on
`rfortran-core`, `rfortran-linalg`, and `numDeriv`; no dependency source is
copied into this directory. The shared `rfortran-linalg` package owns the pinned
`fortran-lapack` dependency, so `pbkrtest/fpm.toml` contains no system BLAS or
LAPACK link directives.

The maintained source was compiled against explicit validation interfaces that
match the shared APIs used by this package. The temporary validation interfaces
and the temporary system-LAPACK link used to exercise them were kept outside the
package directory and are not included in the archive.

## Direct GNU Fortran validation

A fresh temporary build directory was used. All maintained modules were compiled
with:

```text
-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface \
-fcheck=all -fbacktrace
```

The deterministic test program completed with:

```text
pbkrtest tests passed
```

The example completed with:

```text
LRT =    6.00000, df = 2, p =   0.049787
Parametric-bootstrap p =   0.285714
Gamma moment-match p =   0.140513
```

Warnings emitted during this validation were confined to the external temporary
validation stubs (unused dummy arguments and one real equality comparison) and
to the ELF linker executable-stack note caused by an internal test callback.
The maintained package sources compiled without warnings under the flags above.

## Mechanical audits

The pre-package audit verifies that:

- `fpm.toml` parses and contains no `-lblas`, `-llapack`, or BLAS/LAPACK link list;
- all maintained Fortran is free-form and has no code line longer than 132 columns;
- there are no semicolon-separated maintained statements;
- there are no forbidden `double precision`, `real*8`, `kind(0.0d0)`, or `D`-literal spellings;
- there are no self-comparison NaN tests;
- every dummy argument has a standalone declaration with explicit `INTENT` or `VALUE` and a nonempty trailing FORD `!!` comment;
- there are no duplicate Fortran source files;
- no dependency source tree is copied into the package; and
- no object, module, executable, cache, archive, or build-product file is present in the package tree.

## FPM environment limitation

The validation environment did not contain an `fpm` executable. No conda/mamba
or other local FPM installation was available, and outbound DNS/download access
from the build container was unavailable, so FPM could not be installed during
this run. Therefore the literal commands `fpm build`, `fpm test`, and
`fpm clean --all` could not be executed here.

This is recorded explicitly rather than claiming an FPM run that did not occur.
The FPM manifest was parsed mechanically, and the package's maintained source,
test, and example targets were compiled and run directly with GNU Fortran as
described above. The package directory was then manually checked to contain no
build products before archiving.

`fprettify` was also not installed in this environment; source formatting was
checked directly for the maintained style constraints.
