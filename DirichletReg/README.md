# DirichletReg modern Fortran port

This is a computational port of R package **DirichletReg 0.7-2** to modern
Fortran with an FPM project layout. The upstream package is by Marco Johannes
Maier and is licensed under GPL (>= 2). The Fortran port uses
`GPL-2.0-or-later`; the original metadata and source are retained under
`upstream/`.

## Implemented computational functionality

- Dirichlet density / log-density with either a common alpha vector or a
  per-observation alpha matrix.
- Dirichlet random generation using gamma variates.
- Composition preparation corresponding to `DR_data`, including normalization,
  boundary compression, beta-to-two-component conversion, and subcomposition
  collapsing.
- Common Dirichlet regression:
  `alpha_ij = exp(X_j beta_j)`.
- Alternative mean/precision regression:
  softmax mean model with a reference component and
  `phi_i = exp(Z gamma)`, `alpha_ij = mu_ij phi_i`.
- Analytic scores and Hessians for both model parameterizations.
- Starting-value construction based on the strategies used by the R package.
- Two-stage optimization: BFGS followed by Newton-Raphson with backtracking.
- Covariance matrix, standard errors, AIC, BIC, Wald intervals, coefficient
  z-tests, and likelihood-ratio tests.
- Fitted/predicted `alpha`, `mu`, and `phi` values.
- Raw, standardized, and composite residuals.
- Ternary and quaternary simplex coordinate transformations.

Plotting code is intentionally omitted.

## R-specific infrastructure intentionally not ported

The Fortran API works directly with response and design matrices. It does not
attempt to emulate R's `Formula`, model-frame evaluation, S3 classes/printing,
`update()`, or formula-driven `drop1()`. The numerical operation underlying
model comparison is provided as `likelihood_ratio_test`.

For a common model, create one `design_block` per response component. Blocks
may have different numbers of columns. For the alternative model, pass the
common mean-model matrix `X`, the precision-model matrix `Z`, and the base
component number.

## Build and test

```text
fpm test
fpm run --example intercept_model
```

The source is Fortran 2018 compatible and does not require external numerical
libraries.

## Minimal common-model example

```fortran
program demo
  use dirichletreg, only : dp, design_block, dirichletreg_model, fit_common
  implicit none
  real(dp) :: y(5,3)
  type(design_block) :: x(3)
  type(dirichletreg_model) :: model
  integer :: j

  y = reshape([ &
    0.20_dp,0.25_dp,0.22_dp,0.18_dp,0.24_dp, &
    0.30_dp,0.35_dp,0.28_dp,0.32_dp,0.31_dp, &
    0.50_dp,0.40_dp,0.50_dp,0.50_dp,0.45_dp], [5,3])

  do j = 1, 3
    allocate(x(j)%x(5,1))
    x(j)%x = 1.0_dp
  end do

  call fit_common(y, x, model)
  print *, model%loglik
  print *, exp(model%coefficients)
end program demo
```

## Main modules

- `dirichletreg`: umbrella public API.
- `dirichletreg_distribution`: density and RNG.
- `dirichletreg_data`: composition preparation.
- `dirichletreg_common`: common-model likelihood, score, Hessian, prediction.
- `dirichletreg_alternative`: alternative-model likelihood, score, Hessian,
  prediction.
- `dirichletreg_fit`: model fitting.
- `dirichletreg_inference`: residuals and inference utilities.
- `dirichletreg_geometry`: simplex-coordinate transforms.

See `PORTING_NOTES.md` for the source mapping and intentional numerical
corrections.
