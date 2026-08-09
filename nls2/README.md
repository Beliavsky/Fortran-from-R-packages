# nls2-fortran

Modern Fortran/FPM translation of the computational core of the R package **nls2 0.3-4**.

`nls2` extends nonlinear least squares with multiple starting values and search-based start generation. This translation keeps that computational role while replacing R formulas/model objects with typed Fortran callbacks.

## Included

- Gauss-Newton nonlinear least squares with step halving.
- Numerical or analytical Jacobians.
- Optional observation weights.
- Lower/upper bounds for the bounded compatibility path.
- Brute-force/grid start evaluation.
- Uniform random start generation.
- Latin-hypercube start generation.
- Multiple explicit starting rows.
- Best-start selection plus all per-start results.
- Partially-linear variable-projection fitting.
- Partially-linear brute-force/random/LHS start searches.
- Singular-Jacobian start-point evaluation (the useful `nls2` feature that ordinary R `nls` can reject).
- Covariance, residual standard error, log-likelihood, residual df, and Pearson residual helpers.
- Reproducible RNG seeding.

## Build

```text
fpm build
fpm test
```

or, with GNU Fortran and the strict flags used for validation:

```text
scripts/test_gfortran.sh
```

On Windows:

```text
scripts\test_gfortran.bat
```

## Callback model

Instead of an R formula, supply a typed procedure:

```fortran
subroutine model(x, par, yhat, ierr)
    use nls2, only : dp
    real(dp), intent(in) :: x(:,:), par(:)
    real(dp), intent(out) :: yhat(:)
    integer, intent(out) :: ierr
end subroutine model
```

Then, for example:

```fortran
call fit_nls(model, x, y, start, result)
```

or use `nls2_fit` with a matrix of starting values/bounds and one of the search algorithms.

See `API.md` and `TRANSLATION_COVERAGE.md` for details.
