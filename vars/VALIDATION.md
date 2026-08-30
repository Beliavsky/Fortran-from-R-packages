# Validation

Validation was performed on 2026-08-29 in the translation sandbox.

## Required FPM gate

The requested commands were attempted from this `vars/` directory immediately
before packaging:

```text
fpm build
fpm test
fpm clean --all
```

All three attempts returned exit status 127 because the sandbox does not contain
an `fpm` executable.  A bootstrap of the current FPM 0.13.0 release was also
attempted, but outbound downloads from the execution container are disabled and
the release-asset download service rejected the transfer.  No substitute or fake
`fpm` command was used.

`tools/release_check.py` intentionally treats a missing `fpm` as a release-check
failure.  On a normal checkout of `Fortran-from-R-packages` with FPM installed,
it runs `fpm build`, `fpm test`, the source/tree hygiene scans, and finally
`fpm clean --all`.

## Compiler validation performed here

As a secondary build validation, every maintained package source was compiled
with GNU Fortran using strict Fortran 2018 checks:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
-ffree-line-length-none
```

Temporary validation-only interface shims for the documented
`rfortran-core`/`rfortran-linalg` APIs were located outside this package and are
not included in the archive.  They link to the system BLAS/LAPACK installation
where needed.  The interfaces used by this translation were separately checked
against the current shared-module source in the target repository.

Results:

- all `src/*.f90` files compiled without warnings under the strict flags above;
- `test/test_vars.f90` compiled, linked, and printed `All vars tests passed.`;
- `example/fit_var.f90` compiled, linked, and ran successfully;
- `example/impulse_response.f90` compiled, linked, and ran successfully.

## Release hygiene

The final maintained Fortran tree contains 12 unique `.f90` files across
`src/`, `test/`, and `example/`.  The longest maintained source line is 112
characters.  Checks found no:

- duplicate Fortran source files;
- semicolon-separated Fortran statements;
- `double precision`, `real*8`, `kind(0.0d0)`, or `d0`/`D` real literals;
- unsuffixed real constants in maintained Fortran source;
- alternate maintained real-kind definitions;
- vendored `rfortran-*`, BLAS, LAPACK, ARPACK, `r.f90`, or `r_mod.f90` source;
- object/module/library/executable files, caches, build directories, or nested
  archives in the package tree.

All maintained Fortran files import the same `dp` kind from `r_kinds`; the
public `vars` module re-exports `dp`.
