# intrinsicFRP-fortran

A modern Fortran/FPM implementation of the computational core of the R package
`intrinsicFRP` 2.1.0 for evaluating linear factor asset-pricing models.

## Implemented functionality

- Tradable factor risk premia (TFRP), with Newey-West HAC standard errors.
- Fama-MacBeth and Kan-Robotti-Shanken factor risk premia.
- Fama-MacBeth and Gospodinov-Kan-Robotti SDF coefficients.
- GKR iterative factor screening.
- Oracle TFRP with GCV, deterministic k-fold CV, or rolling validation.
- Feng-Giglio-Xiu three-pass factor tests with self-contained Lasso selection.
- Hansen-Jagannathan/Kan-Robotti misspecification distance and interval.
- Kleibergen-Paap-style and Chen-Fang rank tests.
- Giglio-Xiu PCA risk-premium estimation and PCA-count selectors.
- Newey-West HAC covariance, variance, and standard errors, with optional
  VAR(1)/AR(1) prewhitening.

The public API is array based: observations are rows, assets or factors are
columns, and outputs use typed result objects with explicit status codes.

## Build with FPM

```text
fpm build
fpm test
fpm run
```

Examples can be run with:

```text
fpm run --example basic_risk_premia
fpm run --example oracle_and_screening
fpm run --example identification_tests
```

The project has no external FPM dependencies. Dense linear algebra, probability
functions, HAC estimation, deterministic random numbers, and coordinate-descent
Lasso are implemented in Fortran.

## Important adaptations

The upstream package relies on R lists/data frames, RcppArmadillo, and `glmnet`.
This port replaces those interfaces with derived types and self-contained dense
algorithms. The Kleibergen-Paap entry point uses a documented singular-value
asymptotic rank test rather than the complete covariance-weighted upstream test.
FGX Lasso tuning is deterministic and may select a different model from `glmnet`.
See `PORTING.md` for details.

## License

The original package is licensed under GPL version 3 or later. This translation
is distributed under **GPL-3.0-or-later**. See `LICENSE.md` and `NOTICE.md`.
