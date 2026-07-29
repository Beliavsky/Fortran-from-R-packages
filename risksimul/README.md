# risksimul-fortran

Modern Fortran 2018 / FPM translation of the numerical algorithms in the R
package `riskSimul` 0.1.2.

The library estimates portfolio tail-loss probabilities and conditional excess
under a Student-t copula with either Student-t or generalized-hyperbolic (GH)
marginals. It provides both naive Monte Carlo and adaptive stratified
importance sampling.

## Licensing

The upstream package is distributed under `GPL-2 | GPL-3`. This translation
preserves that choice as:

```text
GPL-2.0-only OR GPL-3.0-only
```

See `LICENSE`, `LICENSE-GPL-2.0`, `LICENSE-GPL-3.0`, and `NOTICE`.

## Build

```text
fpm build
fpm test
fpm run risksimul_demo
fpm run --example t_copula_risk
fpm run --example stratified_sampling
```

No external numerical library is required.

## Main API

```fortran
use risksimul
```

Create a t-copula portfolio with Student-t marginals:

```fortran
real(dp) :: correlation(2,2), parameters(2,3)
type(portfolio_model) :: portfolio

correlation = reshape([ &
   1.0_dp, 0.4_dp, &
   0.4_dp, 1.0_dp  &
], [2,2])

! Each row is: location, scale, degrees of freedom.
parameters = reshape([ &
   0.0_dp, 0.0_dp, &
   0.02_dp, 0.03_dp, &
   6.0_dp, 8.0_dp  &
], [2,3])

portfolio = new_portfolio( &
   8.0_dp, correlation, 't', parameters, &
   weight=[0.5_dp,0.5_dp] &
)
```

Naive simulation:

```fortran
type(simulation_result) :: result

result = NVTCopula( &
   50000, portfolio, [0.94_dp,0.97_dp], 12345_i8 &
)
```

Stratified importance sampling with the upstream-style wrapper:

```fortran
result = SISTCopula( &
   10000, [1000,3000], portfolio, &
   [0.94_dp,0.97_dp], &
   stratasize=[6,6], &
   CEopt=.false., &
   beta=0.75_dp, &
   mintype=objective_msre, &
   seed=12345_i8 &
)
```

The more flexible typed interface is `stratified_copula(portfolio,
thresholds, control)`, where `control` is a `sis_control` object.

## Returned statistics

For each threshold, `simulation_result` returns:

- tail-loss probability;
- unconditional excess;
- conditional excess;
- estimated variance;
- normal 95% confidence half-width and interval;
- relative confidence half-width in percent.

The conditional-excess point estimator preserves the finite-sample ratio
correction used by the upstream R package.

## Numerical implementation

The project includes:

- Student-t CDFs and quantiles;
- gamma CDFs and quantiles;
- generalized-hyperbolic density, CDF, quantile, moments, and supporting GIG
  calculations;
- Cholesky factorization and orthogonal-basis completion;
- deterministic optional random seeds;
- rare-event boundary root solving;
- derivative-free direction optimization;
- adaptive multiresponse allocation;
- MSE, relative-MSE, maximum-error, maximum-relative-error, and selected-target
  allocation objectives.

GH inverse CDFs are tabulated once per marginal and then interpolated during
simulation. This replaces the upstream `Runuran::pinvd.new` dependency.

## Source layout

- `src/`: Fortran library
- `test/`: numerical tests
- `app/`: demonstration program
- `example/`: focused usage examples
- `original/`: unmodified upstream package
- `provenance/`: supplied archive and checksums
- `scripts/`: direct GNU Fortran validation scripts

See `COVERAGE.md`, `PORTING_NOTES.md`, and `VALIDATION.md` for details.
