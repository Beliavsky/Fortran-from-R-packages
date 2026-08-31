# Verification

## Compiler regression

The maintained Fortran library, every deterministic test, and the example were
compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror=implicit-interface -fcheck=all -O0 -g
```

Result:

- 7 deterministic test programs passed;
- 1 example program passed;
- no runtime bounds/checking failures occurred;
- the likelihood reference values were independently calculated before being
  encoded in the regression test.

The test set covers pooling/confidence intervals, D1-D4 MI tests, linear and
transformed constraints, cluster means and multilevel R-squared, convergence
and autocorrelation diagnostics, error/status paths, and Gaussian LM/LMM
likelihood helpers.

## Shared-module interface check

The target repository's current public APIs were checked before finalizing the
translation. The package imports `dp` from `r_kinds`, probability functions
from `r_distributions`, and `solve_system`, `solve_spd`, and
`signed_log_determinant` from `r_linalg`.

The repository's sibling `rfortran-core` and `rfortran-linalg` source trees are
not mounted in this execution environment. For the strict local compiler run,
temporary interface-compatible stand-ins were kept outside the package tree;
their procedure signatures match the current repository APIs. No stand-in or
dependency source is present in the deliverable.

## Static source audit

The final maintained source audit covers 17 `.f90` files under `src/`, `test/`,
and `example/` and reports zero issues. Checks include:

- standard free-form line length (132 columns or fewer);
- no duplicate Fortran source contents;
- no semicolon-separated executable statements;
- no `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals;
- no self-comparison NaN idioms;
- each dummy argument declared separately with `INTENT` or `VALUE`;
- every dummy declaration has a nonempty trailing FORD `!!` comment;
- no object/module/executable/cache/build products in the package tree;
- no copied `r.f90`, `r_mod.f90`, BLAS, LAPACK, ARPACK, or dependency source.

The manifest parses successfully as TOML and has only sibling path dependencies
on `../rfortran-core` and `../rfortran-linalg`; it contains no system BLAS or
LAPACK link directive.

## Upstream reference integrity

The retained `upstream/DESCRIPTION`, `NAMESPACE`, `README.md`, and the eight R
files listed in `PROVENANCE.md` were compared byte-for-byte with the uploaded
source archive and match exactly.

## FPM and fprettify availability

A real `fpm` executable is not installed in this runtime. Explicit attempts to
run `fpm build`, `fpm test`, and `fpm clean --all` therefore returned exit 127.
After confirming FPM 0.13.0 as the current stable release, a bootstrap download
was attempted, but outbound DNS is disabled in the container and curl returned
exit 6 (`Could not resolve host: github.com`). The exact command transcript is
in `FPM_ATTEMPTS.txt`; no fake FPM executable was substituted.

Neither the `fprettify` command nor its Python module is installed. The source
was nevertheless kept in fprettify-compatible modern free-form style and was
checked for the formatting constraints required by this translation.
