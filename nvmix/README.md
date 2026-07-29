# nvmix-fortran

A dependency-free modern Fortran/FPM numerical translation of **nvmix 0.1-2**,
a package for multivariate normal variance mixtures, grouped mixtures, Student
distributions and copulas, risk measures, dependence measures and fitting.

## Build

```text
fpm build
fpm test
fpm run nvmix_demo
fpm run --example distributions_and_risk
fpm run --example grouped_mixtures
```

The source is standard free-form Fortran 2018 and uses `implicit none`, double
precision (`dp = real64`) and typed status/result objects.

## Main capabilities

- General and grouped normal variance-mixture models.
- Constant, inverse-gamma, Pareto and mean-one gamma mixing distributions.
- Multivariate normal and Student density, probability and simulation wrappers.
- Grouped Student distributions and copulas with group-specific degrees of freedom.
- Normal-variance-mixture copula density, rectangle probabilities and simulation.
- Gamma-mixture laws for squared Mahalanobis distances.
- Skew-t density, CDF, quantiles, simulation and copula calculations.
- Normal and Student value-at-risk and expected shortfall formulas.
- Generic simulated VaR and ES for one-dimensional mixtures.
- Grouped-mixture correlation, Kendall's tau, Spearman's rho and Student tail dependence.
- Normal, Student and one-parameter variance-mixture fitting.
- Student and grouped-Student copula fitting.
- Numerical Mahalanobis QQ data.

## Typed model construction

```fortran
use nvmix

type(nvmix_model) :: model
real(dp) :: loc(2), scale(2,2)

loc = 0.0_dp
scale = reshape([1.0_dp, 0.5_dp, 0.5_dp, 1.0_dp], [2,2])

model = make_nvmix_model(loc, scale, mix_inverse_gamma, 7.0_dp)

print *, dnvmix([0.0_dp, 0.0_dp], model)
```

A grouped model is created with `make_grouped_model`. Each coordinate is mapped
to a mixing group, and each group has its own family and parameter.

```fortran
integer :: groupings(3), families(2)

groupings = [1, 1, 2]
families = mix_inverse_gamma

model = make_grouped_model(loc3, scale3, groupings, families, [6.0_dp, 14.0_dp])
```

## Integration controls

Probability and non-closed-form density calculations accept an
`integration_control` object:

```fortran
type(integration_control) :: control

control%samples = 32768
control%batches = 16
control%seed = 12345_i8
```

The Fortran implementation uses deterministic shifted Halton integration and
batch error estimates. It does not require `qrng`, `mvtnorm`, `mnormt`,
`Matrix`, `copula`, `pcaPP`, `ADGofTest` or `pracma`.

## Compatibility names

The library provides familiar names including:

- `dnvmix`, `pnvmix`, `rnvmix`, `qnvmix`
- `dgnvmix`, `pgnvmix`, `rgnvmix`
- `dNorm`, `pNorm`, `rNorm`, `fitNorm`
- `dStudent`, `pStudent`, `rStudent`, `fitStudent`
- `dgStudent`, `pgStudent`, `rgStudent`
- Student and grouped-Student copula wrappers
- `dgammamix`, `pgammamix`, `qgammamix`, `rgammamix`
- `VaR_nvmix`, `ES_nvmix`, `corgnvmix`, `lambda_gStudent`

Fortran is case-insensitive, so capitalization is documentary.

## Scope differences

The original package uses adaptive randomized quasi-Monte Carlo, Sobol and
generalized Halton point sets, specialized C integrands, extensive R function
callbacks and third-party fitting/diagnostic packages. This port uses a single
self-contained deterministic Halton engine. Results should agree within the
reported integration accuracy but are not expected to be bit-for-bit identical.

The grouped-Student tail-dependence routine uses the exact formula when the two
degrees of freedom agree and a symmetric effective-degrees-of-freedom
approximation otherwise. Exact grouped tail dependence in the R package uses
RQMC.

The copula fitting routines are dependency-free numerical estimators rather
than replicas of every R optimization option. See `COVERAGE.md` and
`PORTING_NOTES.md` for details.

## License

GPL-3.0-or-later. The original authors and source are retained for provenance.
