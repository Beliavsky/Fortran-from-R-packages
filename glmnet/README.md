# glmnet-fortran

Modern Fortran/FPM implementation of the computational core of the R package
`glmnet` 5.0, **Lasso and Elastic-Net Regularized Generalized Linear Models**.

The library fits regularization paths for:

- Gaussian linear regression
- Binomial logistic regression
- Poisson log-linear regression
- Multinomial logistic regression, with grouped or ungrouped penalties
- Multiresponse Gaussian regression with a grouped elastic-net penalty
- Cox proportional-hazards regression for right-censored or `(start, stop]` data

It also provides prediction and coefficient interpolation, cross-validation,
relaxed refitting, model assessment, ROC/AUC and confusion tables, Cox
partial-likelihood gradients and C-index, missing-value utilities, deterministic
fold generation and multinomial sampling, and dense/CSC conversion.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

The package uses a plain semantic version (`5.0.0`) accepted by FPM.

## Minimal example

```fortran
program example
   use glmnet, only : dp, glmnet_control_type, glmnet_path_result, fit_glmnet
   implicit none
   real(dp) :: x(100, 5), y(100)
   type(glmnet_control_type) :: control
   type(glmnet_path_result) :: fit

   ! Fill x and y.
   control%alpha = 0.8_dp
   control%nlambda = 50
   call fit_glmnet(x, y, 'gaussian', fit, control)

   print *, fit%lambda(fit%nlambda)
   print *, fit%beta(:, 1, fit%nlambda)
end program example
```

See `example/` and `app/demo_glmnet.f90` for complete programs.

## Design

The public API is array based. R lists, S3 classes, formula processing, factors,
and sparse Matrix classes are represented by Fortran derived types and explicit
procedures. Coefficients are always returned on the original predictor scale.

The main result type, `glmnet_path_result`, contains:

- `lambda`
- `intercept(nout, nlambda)`
- `beta(nvars, nout, nlambda)`
- `dev_ratio`, `objective`, and `df`
- iteration and convergence diagnostics
- standardization metadata

IEEE NaNs are not accepted by fitting routines. Use `na_replace` or `prepare_x`
first when data contain missing values.

## Numerical scope

This port preserves the elastic-net objectives, warm-started paths, IRLS for
binomial and Poisson models, grouped penalties for multiresponse and
multinomial models, Breslow/Efron Cox partial likelihoods, coefficient bounds,
penalty factors, exclusions, offsets, observation weights, and the usual
lambda sequence.

The implementation is intentionally self-contained and portable. It is not a
line-for-line translation of the current templated C++ engine and does not
claim its sparse-data speed. Important adaptations are listed in `PORTING.md`.

## Validation

Run either:

```text
scripts/test_gfortran.sh
scripts/test_gfortran_optimized.sh
```

The strict configuration uses Fortran 2018, warnings as errors, bounds and
runtime checking, floating-point traps, and backtraces. The optimized
configuration uses `-O3` with warnings as errors.

## License

The upstream package declares `GPL-2`. This translation is distributed under
**GPL-2.0-only**. The license text is in `LICENSE`; retained upstream material
is under `original/glmnet-master/`.
