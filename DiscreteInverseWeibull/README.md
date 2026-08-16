# DiscreteInverseWeibull-fortran

Modern Fortran/FPM translation of the computational code in CRAN package `DiscreteInverseWeibull` 1.0.2.

The library provides the discrete inverse Weibull PMF/CDF/quantile/RNG, hazard and accumulated hazard, moment calculations, likelihood/loss functions, the heuristic estimator, and all four upstream estimation methods (`P`, `M`, `H`, `PP`). The method-of-moments fit uses the supplied `Rsolnp-fortran` translation as an FPM dependency.

```fortran
use discrete_inverse_weibull
real(dp) :: x(100)
type(diw_estimate) :: fit
call rdiweibull(x,0.5_dp,2.5_dp)
fit = estdiweibull(x,'H')
```

See `PORTING_NOTES.md` and `API_MAP.md` for compatibility details.
