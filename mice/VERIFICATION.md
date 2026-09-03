# Verification

## Environment

- GNU Fortran: 14.2.0
- Language mode: Fortran 2018
- Runtime checking: `-fcheck=all`
- Interface diagnostics: `-Werror=implicit-interface`
- Additional warnings: `-Wall -Wextra`
- `fprettify`: not installed in this execution environment

## Strict compiler regression

All maintained library modules were compiled from a clean temporary build
directory with:

```text
-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all -O0
```

The target repository's `rfortran-core` and `rfortran-linalg` directories were
not present in this container. For this compiler-only regression, minimal
interface-compatible `r_kinds`/`r_linalg` test stubs were kept **outside** the
package tree. Their signatures match the shared `dp`, `inverse_matrix`, and
`solve_spd` interfaces used by this translation. No stub or dependency source
is included in the deliverable.

Results:

```text
test_2l_norm_full  PASS
test_ampute        PASS
test_categorical   PASS
test_diagnostics   PASS
test_fcs           PASS
test_fcs_parity    PASS
test_native        PASS
test_parity_methods PASS
test_pmm           PASS
test_pooling       PASS
test_regression    PASS
test_twolevel      PASS
mice_example       PASS
mice_parity_example PASS
```

Total: **12/12 deterministic tests and 2/2 examples passed**.

The parity-specific tests cover the dedicated proportional-odds fit,
`midastouch`, `mpmm`, LDA, NARFCS offsets, heterogeneous `2l.norm`, and FCS
dispatch of the newly added univariate methods.

## FPM attempts

The required commands were attempted from the package root. This runtime does
not have the Fortran Package Manager executable installed, so each command
returned exit code 127 (`fpm: command not found`):

```text
fpm build
fpm test
fpm clean --all
```

The exact captured output is retained in `FPM_ATTEMPTS.txt`. No replacement or
mock `fpm` executable was used.

## Static/package audit

The final maintained-source audit checks library, tests, and examples for:

- free-form line length at or below 132 columns;
- duplicate Fortran source content;
- executable semicolon-separated statements;
- `double precision`, `real*8`, `kind(0.0d0)`, and D-exponent literals;
- self-comparison NaN idioms;
- unsafe finite-math/fast-math flags;
- each dummy argument declared separately with explicit `INTENT` or `VALUE`;
- a meaningful trailing FORD `!!` comment on every dummy declaration;
- copied `r.f90`, `r_mod.f90`, BLAS, LAPACK, ARPACK, `r_kinds`, or `r_linalg`
  dependency sources;
- object/module/executable/cache/ZIP build products in the package tree.

Result: **33 Fortran files audited; 0 source-policy violations**.

The retained upstream metadata and direct computational reference files were
also compared byte-for-byte with the supplied `mice-master.zip`: **31/31
reference files matched exactly**.
