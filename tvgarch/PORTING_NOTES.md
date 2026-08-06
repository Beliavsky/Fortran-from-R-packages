# Porting notes

## Parameter order

The short-term GARCH-X vector follows the supplied dependency:

```text
intercept, ARCH coefficients, GARCH coefficients,
asymmetry coefficients, variance-regressor coefficients
```

The R `order.h` convention is retained by `make_garch_order_spec`:

```text
[p = GARCH order, q = ARCH order, r = asymmetry order]
```

The TV vector is:

```text
intercept.g, sizes(1:s), speeds(1:s), locations(1:sum(order.g))
```

For `speed_option=2`, speed parameters are stored on the log scale, matching
the R package.

## Estimation

`fit_tvgarch` follows the package's maximization-by-parts design:

1. Fit the TV component with `h=1`.
2. Hold the TV intercept fixed.
3. Alternate GARCH-X fitting on `y/sqrt(g)` and TV fitting conditional on `h`.
4. Refit the GARCH-X component once more.
5. Optionally jointly refine all free TV and GARCH parameters.

The R package uses linear constraints in `constrOptim`. The Fortran port uses
box-constrained Nelder-Mead and rejects candidates that violate positivity,
location ordering, or positive TV regimes. This is numerically equivalent in
model intent but can produce different local optima.

## Multivariate implementation

Marginal variances are fitted equation by equation. Constant correlations are
sample correlations of standardized residuals. DCC follows the package
recursion

```text
Q(t) = (1-a-b) Rbar + a z(t-1)z(t-1)' + b Q(t-1)
```

with the same first-100-observation initialization when available.

The spillover routine iterates lagged squared returns standardized by the TV
component. It accepts a logical equation-by-source inclusion matrix rather than
R's flattened binary vector.

## Known differences

- Numeric arrays replace `zoo` and R names/indexes.
- Gaussian random streams differ from R.
- Optimizer endpoints and finite-difference standard errors are not bitwise R-equivalent.
- The multivariate API uses arrays of typed marginal specifications instead of padded parameter matrices.
- `tvgarchTest` returns numeric tables rather than an S3 object.
