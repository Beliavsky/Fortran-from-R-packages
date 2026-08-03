# Porting notes

## Array layout

Return matrices use `(observation, asset)`. Weight and expected-return vectors
have length `nassets`; covariance matrices are `nassets x nassets`. Random
portfolio and frontier weight collections use `(asset, portfolio)`.

## Constraint handling

`repair_weights` exactly projects onto box bounds and a selected total-weight
level. Cardinality constraints are applied by retaining the largest absolute
positions. Group, factor, return, leverage, turnover, and transaction-cost
conditions are evaluated explicitly and enter generic optimization through a
quadratic violation penalty. Differential evolution is therefore preferred for
strongly nonlinear or nonconvex specifications.

## Optimization

`opt_auto` selects projected gradient for minimum variance, maximum return, and
quadratic utility, and differential evolution for tail-risk, ratio, drawdown,
and risk-budget objectives. Solver controls are contained in
`portfolio_options`. The result always reports both convergence and final
constraint feasibility.

## Risk conventions

Historical VaR and ES are positive loss numbers computed from the lower return
tail. `alpha` is the lower-tail probability, such as `0.05`. Volatility risk
contributions sum to portfolio standard deviation, subject to floating-point
roundoff.

## Entropy pooling

The entropy solver uses the exponential-family dual and a Newton active-set
method. A probability-sum equality is inserted automatically when absent.
Inequalities use the upstream `A p <= b` convention.

## Higher moments

Coskewness is flattened as `(asset_i, (asset_j-1)*n + asset_k)` and cokurtosis
as `(asset_i, ((asset_j-1)*n + asset_k-1)*n + asset_l)`, matching the matrix
shape used by PortfolioAnalytics/PerformanceAnalytics.

## Numerical scope

The project is self-contained and intentionally avoids mandatory BLAS/LAPACK or
external solver linkage. Its dense linear algebra and population optimizers are
suited to portfolio research and moderate asset counts; large institutional
problems should connect the typed model layer to a dedicated QP/conic/MIP
backend.
