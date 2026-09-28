# Validation record

The maintained Fortran source is validated with GNU Fortran using two independent local builds against a temporary out-of-tree interface adapter to the system LAPACK. That adapter is a validation aid only and is **not** part of the package or ZIP.

Checked build flags:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Optimized build flags:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface
```

Deterministic tests cover:

- real symmetric-definite and general generalized eigenproblems;
- complex Hermitian-definite generalized eigenproblems;
- real and complex generalized Schur reconstruction;
- QZ sorting and generalized-eigenvalue extraction;
- real, complex, and wide rectangular GSVD identities;
- `gsvd.R`, `gsvd.oR`, `gsvd.D1`, and `gsvd.D2` helpers;
- the basic example.

The GSVD tests validate `U^H A Q = D1 [0 R]`, `V^H B Q = D2 [0 R]`, and `alpha^2 + beta^2 = 1` on the active rank.

## FPM limitation in this sandbox

The sandbox used to prepare this checkpoint has GNU Fortran but no `fpm` executable. Literal attempts to run `fpm build`, `fpm test`, `fpm run --example basic`, and `fpm clean --all` therefore return exit code 127 (`fpm: command not found`). The package manifest is nevertheless configured to use the Git `fortran-lapack` dependency rather than system `-llapack/-lblas` links.

On a normal checkout with FPM installed, run:

```text
fpm clean --all
fpm build
fpm test
fpm run --example basic
```
