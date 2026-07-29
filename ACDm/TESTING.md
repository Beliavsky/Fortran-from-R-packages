# Testing

## Test programs

`test/test_distributions.f90`

- fixed independent SciPy density, CDF, and quantile references
- CDF/quantile inversion for all distribution families
- direct component-distribution identities

`test/test_models.f90`

- simulation, filtering, and finite likelihood for all 13 model recursions
- threshold and spline model branches
- exogenous regressors
- daily recursion resets
- fixed ACD recurrence references

`test/test_fit.f90`

- simulated ACD(1,1) maximum-likelihood recovery
- Hessian, covariance, standard errors, and robust covariance
- AIC/BIC and deterministic forecasting

`test/test_data_diagnostics.f90`

- trade, price, and volume duration construction
- duplicate-timestamp aggregation
- cubic spline, penalized spline, adaptive super smoother, and FFF adjustment
- grouped diurnal adjustment
- PIT and Cox-Snell residuals
- ACF, KDE, QQ, summary, and rolling statistics
- robust and nonrobust remaining-ACD, STACD, and TVACD tests

`test/test_profiles.f90`

- nonparametric and fitted hazard coordinates
- likelihood-profile grids

## GNU Fortran configurations

The release script supports two clean configurations.

Strict:

```text
-std=f2018 -O0 -g
-Wall -Wextra -Werror -Wimplicit-interface -Wconversion-extra
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

Optimized:

```text
-std=f2018 -O3 -march=x86-64
-Wall -Wextra -Werror -Wimplicit-interface
```

Run:

```text
./run_gfortran_tests.sh strict
./run_gfortran_tests.sh optimized
```

The script builds the library from a clean directory, compiles all tests,
applications, and examples, runs every target, and stops on the first failure.
