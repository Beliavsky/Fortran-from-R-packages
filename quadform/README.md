# quadform-fortran

Modern free-format Fortran translation of the computational code in Robin K. S. Hankin's R package `quadform` 0.0-4.

The upstream package is a compact collection of efficient quadratic-form identities. This translation preserves the same real/complex conjugation semantics and multiplication-order variants while exposing Fortran-friendly names.

## Features

- Hermitian transpose: `ht`
- conjugate cross-products: `cprod` / `cp`
- conjugate transpose-cross-products: `tcprod` / `tcp`
- quadratic forms: `quad_form` / `qf`
- inverse quadratic forms: `quad_form_inv` / `qfi`
- three-argument forms: `quad3_form`, `quad3_form_inv`, `quad3_tform`
- transposed-orientation forms: `quad_tform`, `quad_tform_inv`
- diagonal-only and trace-only evaluation
- Cholesky-factor quadratic form
- real and complex double-precision generic interfaces
- self-contained pivoted linear solve for inverse forms

The R names containing dots are mapped to underscore-separated Fortran names. The terse aliases exported by the R package (`qf`, `qfi`, `q3`, `q3i`, etc.) are also available.

## Build

```text
fpm build
fpm test
fpm run --example demo_quadform
```

The package is dependency-free.

## Numerical conventions

For complex matrices, `quad_form(M,X)` computes

```text
conjg(transpose(X)) * M * X
```

and does **not** conjugate `M`, matching the upstream R package.

Inverse routines return an empty `0 x 0` matrix and a nonzero optional `info` value when the coefficient matrix is singular or dimensionally invalid.

## Upstream

The complete source snapshot supplied for translation is retained under `upstream/quadform-master/`.
