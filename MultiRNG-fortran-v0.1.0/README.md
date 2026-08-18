# MultiRNG-fortran

Modern, self-contained, free-format Fortran translation of the computational code in
R package **MultiRNG 1.2.4**.

## Requirements

- Fortran 2018 compiler
- FPM (Fortran Package Manager), for normal package use

There are no C/C++ or external numerical-library dependencies.

## Build

```text
fpm build
fpm test
fpm run --example example_multirng
```

All compiled Fortran sources are free-format `.f90`, use `implicit none`, and are kept
within the standard free-form line length.

## Main API

The public `multirng` module exports:

- `seed_rng`
- `draw_d_variate_normal`
- `draw_d_variate_t`
- `draw_d_variate_uniform`
- `draw_correlated_binary`
- `draw_dirichlet`
- `draw_multinomial`
- `draw_dirichlet_multinomial`
- `draw_multivariate_hypergeometric`
- `draw_multivariate_laplace`
- `draw_wishart`, `draw_wishart_flat`
- `draw_inv_wishart`, `draw_inv_wishart_flat`
- `draw_inv_wishart_legacy`
- `generate_point_in_sphere`
- `loc_min`

The Wishart routines return a natural Fortran rank-3 array `(replicate,d,d)`; `*_flat`
routines provide the upstream row-flattened representation.

## Inverse-Wishart note

The upstream R implementation of `draw.inv.wishart()` computes a Wishart draw using
the inverse of `inv.sigma` and returns that matrix without inverting the draw. That does
not generate the inverse-Wishart density documented in the package manual.

`draw_inv_wishart()` implements the documented distribution. For strict reproduction
of the upstream R code, use `draw_inv_wishart_legacy()`.

See `TRANSLATION_NOTES.md` for other behavioral details preserved from the R code.
