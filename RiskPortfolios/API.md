# Public API

All public entities are available through:

```fortran
use riskportfolios
```

All real-valued APIs use `real(dp)`.

## Mean estimation

```fortran
subroutine mean_estimation(rets, mu, method, lambda, info)
```

Arguments:

- `rets(:, :)`: observations by assets
- `mu(:)`: allocated output vector
- `method`: optional method constant
- `lambda`: optional EWMA decay, default 0.94
- `info`: optional status, zero on success

Method constants:

- `MEAN_NAIVE`
- `MEAN_EWMA`
- `MEAN_BAYES_STEIN`
- `MEAN_MARTELLINI`

Direct procedures are also public:

- `naive_mean`
- `ewma_mean`
- `bayes_stein_mean`
- `martellini_mean`

## Covariance estimation

```fortran
subroutine covariance_estimation(rets, sigma, method, lambda, n_factors, info)
```

Method constants:

- `COV_NAIVE`
- `COV_EWMA`
- `COV_LEDOIT_WOLF`
- `COV_FACTOR`
- `COV_CONSTANT`
- `COV_COR_SHRINKAGE`
- `COV_ONE_PARAMETER`
- `COV_DIAGONAL`
- `COV_LARGE`
- `COV_BAYES_STEIN`

`lambda` defaults to 0.94. `n_factors` defaults to 1.

Direct estimator procedures are public under corresponding descriptive names.

## Semideviation estimation

```fortran
subroutine semideviation_estimation(rets, semidev, method, lambda, info)
```

Method constants:

- `SEMIDEV_NAIVE`
- `SEMIDEV_EWMA`

Direct procedures:

- `naive_semideviation`
- `ewma_semideviation`

The observations must be ordered from oldest to newest for EWMA methods.

## Portfolio optimization

```fortran
subroutine optimal_portfolio(sigma, weights, portfolio_type, mu, semidev, &
   control, info)
```

Portfolio constants:

- `PORT_MEAN_VARIANCE`
- `PORT_MINIMUM_VARIANCE`
- `PORT_INVERSE_VOLATILITY`
- `PORT_EQUAL_RISK_CONTRIBUTION`
- `PORT_MAXIMUM_DIVERSIFICATION`
- `PORT_RISK_EFFICIENT`
- `PORT_MAXIMUM_DECORRELATION`

Constraint constants:

- `CONSTRAINT_NONE`
- `CONSTRAINT_LONG_ONLY`
- `CONSTRAINT_GROSS`
- `CONSTRAINT_USER`

Control type:

```fortran
type(portfolio_control) :: control
```

Fields:

- `constraint`: constraint constant
- `gross_limit`: L1 gross-exposure limit, default 1.6
- `gamma`: mean-variance risk-aversion parameter, default 0.8773
- `max_iterations`: optimizer iteration limit, default 4000
- `tolerance`: optimizer tolerance, default 1.0e-10
- `lower_bounds(:)`: optional asset lower bounds
- `upper_bounds(:)`: optional asset upper bounds
- `initial_weights(:)`: optional starting portfolio

`mu` is required by mean-variance optimization. `semidev` is required by the
risk-efficient method.

## Risk contributions

```fortran
rc = portfolio_risk_contributions(sigma, weights, proportional)
```

The default returns proportional variance contributions summing to one.
With `proportional=.false.`, it returns volatility contributions.

## Linear algebra and statistical utilities

The following utility routines are public because they are useful when
validating inputs and outputs:

- `solve_linear`
- `symmetric_eigen`
- `is_positive_definite`
- `sample_covariance`
- `population_covariance`
- `standard_deviations`
- `covariance_to_correlation`
- `quantile_type7`
- `median_value`
