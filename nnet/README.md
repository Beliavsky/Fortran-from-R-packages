# nnet-fortran

Modern Fortran/FPM translation of the computational core of R package `nnet`
7.3-21 by W. N. Venables and B. D. Ripley.

The port provides the single-hidden-layer feed-forward neural-network engine and
the multinomial-log-linear estimator that `nnet::multinom` builds on top of it.
It is matrix oriented and does not attempt to reproduce R formula/model-frame or
S3 presentation machinery.

## Implemented numerical functionality

- Standard `nnet` topology with one hidden layer, bias weights, and optional
  input-to-output skip connections.
- Logistic hidden/output units and linear output units.
- Squared-error, entropy, softmax cross-entropy, and censored-softmax losses.
- Case weights, per-weight decay, initial weights, fixed/free weight masks,
  `rang`, iteration/tolerance controls, and the upstream `MaxNWts` guard.
- Exact analytic back-propagation gradient translated from `nnet.c`.
- Exact Ripley analytic Hessian translated from `nnet.c`.
- Raw and class-index prediction.
- Binary and multiclass multinomial-log-linear fitting from labels or a response
  count/proportion matrix.
- Multinomial offsets with the same fixed-unit-coefficient construction used by
  the R package.
- Censored multinomial fitting.
- Multinomial Fisher information, SVD generalized-inverse covariance, deviance,
  AIC, log-likelihood, rank, and effective degrees of freedom.
- `class.ind`, randomized tie-breaking `which.is.max`, and the native `summ2`
  duplicate-design-row summation kernel.

The supplied MIT-licensed `r_mod.f90` is used for the R-compatible BFGS
optimizer, RNG helpers, and QR rank calculation.  No independent replacement
of those helpers is included.

## Build

```text
fpm build
fpm test
fpm run --example basic_nnet
```

The generalized inverse used for `vcov.multinom` delegates to the local
`rfortran-linalg` dependency and its pinned pure-Fortran LAPACK backend.
System BLAS and LAPACK are not required.

## Minimal example

```fortran
use nnet, only: dp, nnet_model_t, nnet_fit, nnet_predict

type(nnet_model_t) :: fit
real(dp) :: x(8,1), y(8,1), decay(1)
real(dp), allocatable :: pred(:,:)

! fill x and y ...
decay = 0.0_dp
call nnet_fit(fit, x, y, hidden_size=0, linout=.true., skip=.true., &
   decay=decay, maxit=200)
pred = nnet_predict(fit, x)
```

For `multinom_fit_labels`, the `x` argument is a design matrix.  Include an
intercept column if an intercept is desired, just as R's `multinom()` receives
an intercept column from `model.matrix` while fixing the network bias weight.

See `API_MAPPING.md` and `PORTING_NOTES.md` for exact scope and differences.
