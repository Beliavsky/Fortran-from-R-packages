# Build and verification record

This record describes validation performed on the maintained Fortran translation
of MCMCglmm 2.36 in this execution environment.

## Compiler validation

Compiler: GNU Fortran 14.2.0.

The complete maintained source, deterministic test program, and every example
were compiled in two modes using local validation representations of the
repository sibling dependencies. The validation representation of
`rfortran-linalg` matches its public interfaces and uses host LAPACK only as an
external test scaffold; no scaffold source or system BLAS/LAPACK link is present
in the MCMCglmm package.

Checked mode:

```text
-std=f2018 -O0 -fcheck=all -fbacktrace
-Wall -Wextra -Werror -Wimplicit-interface -Wsurprising -pedantic
```

Result: the complete deterministic test suite passed and all 14 examples
compiled and ran successfully.

Optimized mode:

```text
-std=f2018 -O2
-Wall -Wextra -Werror -Wimplicit-interface -Wsurprising -pedantic
```

Result: the complete deterministic test suite passed and all 14 examples
compiled and ran successfully.

No `-ffree-line-length-none` extension was used. All maintained Fortran source
lines are within the standard 132-column free-form limit.

## Shared dependency interface verification

The package manifest uses only sibling path dependencies:

```text
rfortran-core = { path = "../rfortran-core" }
rfortran-linalg = { path = "../rfortran-linalg" }
ape = { path = "../ape" }
```

The maintained source follows the repository's actual shared API contracts,
including the preallocated-output `solve_spd` interface and the impurity of
LAPACK-backed symmetric eigensolver wrappers. The package contains no local copy
of those shared modules.

## FPM and fprettify availability

Immediately before final packaging the requested commands are re-attempted from
the package root:

```text
fpm --version
fpm build
fpm test
fpm clean --all
fprettify --version
```

All four FPM commands returned status 127 (`fpm: command not found`) in this
environment. `fprettify --version` likewise returned status 127. Consequently
this record does not substitute direct gfortran validation for an actual FPM or
fprettify run.

## Source-policy audit

Before packaging, every maintained `.f90` file under `src/`, `test/`, and
`example/` is audited for:

- lines longer than 132 columns;
- semicolon-separated executable/declaration statements;
- legacy `double precision`, `real*8`, `kind(0.0d0)`, or `d0`/`D0` real forms;
- self-comparison NaN idioms;
- dummy arguments lacking `INTENT` or `VALUE`;
- multiple dummy arguments on one declaration;
- dummy declarations lacking meaningful trailing `!!` FORD comments;
- duplicate Fortran filenames or byte-identical duplicate source files;
- fixed-form maintained source;
- copied `r_kinds`, `r_linalg`, R compatibility, BLAS/LAPACK, ape, or CSparse
  dependency source;
- object/module/executable/library/cache/ZIP build artifacts inside the package;
- system `-lblas`/`-llapack` links or unsafe fast-math flags in build metadata.

Freeze audit result: 44 maintained `.f90` files (29 under `src/`, one under
`test/`, and 14 under `example/`) with zero policy violations. The manifest also
parsed successfully as TOML.

## Upstream metadata verification

The preserved files must remain byte-identical to the uploaded upstream archive:

- `upstream/DESCRIPTION` == upstream `DESCRIPTION`;
- `upstream/NAMESPACE` == upstream `NAMESPACE`;
- `upstream/CITATION` == upstream `inst/CITATION`.

The upstream archive does not contain a root `COPYING` file. The translation's
`COPYING` contains the GPL version 2 license text corresponding to upstream's
`GPL (>= 2)` declaration; `NOTICE.md` records the `-or-later` grant and source
provenance.

## Final archive validation

After creating the release ZIP, its exact extracted contents are rebuilt in the
same checked and optimized modes, all tests and examples are rerun, and the
source/provenance audit is repeated. The final response reports the outcome and
SHA-256 of that exact archive.
