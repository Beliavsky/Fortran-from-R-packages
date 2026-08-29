# Porting notes

## Scope

The source package has no compiled code. Its computational core is a small set of matrix utilities plus matrix-normal density, rectangle probability, and random generation. The Fortran port therefore consists of two numerical modules and an umbrella module:

* `matrix_normal_utils` -- trace, vectorization, half-vectorization, special matrices, Kronecker products, and symmetry/definiteness checks.
* `matrix_normal_distribution` -- validation, density, rectangle probabilities, and simulation.
* `matrixNormal` -- public umbrella module.

R documentation generation, data-frame coercion, dimnames/names, warnings, package-version checks, and other language-specific presentation logic are not translated.

## Matrix-normal convention

The package documentation states

`vec(A) ~ N(vec(M), V kron U)`

for `A ~ MatNorm(M,U,V)`, with R/Fortran column-major `vec`. The density formula and the corrected upstream `rmatnorm()` implementation are consistent with this convention.

### `pmatnorm()` covariance-order correction

The current R source nevertheless calls `mvtnorm::pmvnorm(..., sigma=kronecker(U,V))`. This conflicts with the documented distribution and with `rmatnorm()`, whose source explicitly marks `kronecker(U,V)` as incorrect and uses `kronecker(V,U)`.

The Fortran `pmatnorm()` therefore uses `V kron U` by default. To reproduce the current R implementation exactly for regression work, pass

`legacy_covariance_order=.true.`

which uses `U kron V`.

## `vech()` correction

The R source computes

```
stack <- A[lower.tri(A, diag=TRUE)]
vech.A <- full[stack]
```

so the numerical *values* in the lower triangle are used as indices into `vec(A)`. That only gives the intended result accidentally for some matrices. The Fortran routine implements the documented mathematical operation directly: it stacks the lower triangle, including the diagonal, column by column.

## Density implementation

The R code explicitly forms `solve(U)` and `solve(V)` and evaluates

`tr(U^{-1}(A-M)V^{-1}(A-M)^T)`.

The Fortran implementation preserves that expression numerically but avoids explicit inverse matrices. It uses SPD solves and Cholesky log determinants from the supplied mvtnorm Fortran dependency. This is algebraically equivalent and more stable.

## Random generation

The upstream R implementation forms `V kron U`, vectorizes `M`, calls `mvtnorm::rmvnorm`, and reshapes the result. The Fortran implementation follows the same construction through the attached `mvtnorm` port. The dependency port automatically uses Cholesky decomposition with an eigen fallback, so the R string selector `method="chol"/"eigen"/"svd"` is not reproduced.

The upstream argument `s` is documented as a placeholder and currently has no effect. The Fortran generic keeps a one-draw form and additionally makes the integer-`s` form useful by returning `s` draws in a rank-3 array.

## Probability result

The attached mvtnorm port returns a `probability_result` object rather than R attributes. `pmatnorm()` exposes that object directly, preserving probability value, estimated numerical error, status code, evaluation count, and diagnostic message.

## Dependency and licensing note

The matrixNormal-derived sources are GPL-3.0-only. The user-supplied mvtnorm Fortran dependency identifies itself as GPL-2.0-only and is retained unchanged under `vendor/mvtnorm-fortran` with its own license and notices. This archive does not relicense that dependency. The two codebases remain separately identified; users redistributing linked artifacts should review the applicable license terms for their use case.
