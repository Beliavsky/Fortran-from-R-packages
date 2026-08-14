# Validation

## Strict build

The full dependency and mlr source tree was compiled directly with gfortran
using the equivalent of:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

The FPM manifest explicitly disables implicit typing and implicit external
procedures.

## Permanent tests

Seven test programs plus the example cover:

1. regression/classification measures, AUC, Brier/log loss;
2. exact linear regression, logistic regression, k-means and k-NN;
3. standardization, missing-value imputation, k-fold partitions and SMOTE;
4. grid/random tuning and exhaustive/forward feature selection;
5. Cox fitting and concordance through survival-fortran;
6. callback-driven 5-fold regression evaluation;
7. binary classification threshold tuning.

## Independent numerical cross-checks

A separate validation driver was compared with NumPy/SciPy on deterministic
random problems generated with NumPy's PCG64:

- 50 ordinary linear-regression problems (20-79 rows, 1-5 predictors): maximum
  absolute coefficient difference from `numpy.linalg.lstsq` was about
  `3.55e-15`.
- 30 nonseparable binary logistic-regression problems (80-179 rows, 1-4
  predictors): all Fortran Newton fits converged. Against independent SciPy
  BFGS minimization of the binomial negative log likelihood, the maximum
  objective difference was about `5.68e-14`; maximum coefficient difference
  was about `8.02e-8`.

The exact linear example and the callback-driven CV example recover zero RMSE
up to floating-point roundoff.
