# AdequacyModel-fortran

Modern Fortran 2018 translation of the computational code in the R package
**AdequacyModel 2.0.0** by Pedro Rafael Diniz Marinho, Marcelo Bourguignon,
and Cicero Rafael Barros Dias.

The project uses the Fortran Package Manager (FPM) layout and is licensed
under GPL-2.0-or-later, matching the upstream package.

## Scope

Translated computational areas:

- `pso()` -> `pso_optimize`
- `descriptive()` -> `descriptive`
- numerical part of `TTT()` -> `ttt_curve`
- `goodness.fit()` -> `goodness_fit` and `goodness_from_mle`
- local optimization choices used through R `optim()`:
  - BFGS
  - Nelder-Mead
  - nonlinear conjugate gradient
  - simulated annealing
- one-sample Kolmogorov-Smirnov statistic and finite-sample p-value
- numerical Hessians and Hessian-based standard errors
- PDF normalization checks on finite and infinite domains

The plotting performed by `TTT()` and the histogram graphics side effect of
`descriptive()` are intentionally omitted. `ttt_curve` returns the curve
coordinates so a caller can plot them using any graphics library.

## Build

```sh
fpm build
fpm test
fpm run --example example_normal
```

The code has also been validated directly with GNU Fortran 14 using:

```sh
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

## Public API

The facade module is:

```fortran
use adequacy_model
```

### Descriptive statistics

```fortran
type(descriptive_result) :: s
s = descriptive(x)
```

`descriptive_result` contains mean, median, modes, sample variance,
skewness, excess kurtosis, minimum, maximum, and sample size.

### Total-time-on-test curve

```fortran
real(dp), allocatable :: r(:), t(:)
call ttt_curve(x, r, t)
```

This computes exactly the numerical curve used by the upstream `TTT()`
function but performs no plotting.

### Particle swarm optimization

```fortran
type(optimize_result) :: opt
call pso_optimize(objective, data, lower, upper, opt)
```

The objective callback has the form:

```fortran
function objective(par, data) result(value)
    real(dp), intent(in) :: par(:), data(:)
    real(dp) :: value
end function objective
```

### Goodness of fit

For supplied parameters:

```fortran
type(goodness_result) :: g
call goodness_from_mle(pdf, cdf, par, data, g, domain)
```

For estimation plus diagnostics:

```fortran
call goodness_fit(pdf, cdf, starts, data, g, method='BFGS', domain=domain)
```

For `method='PSO'`, `lower=` and `upper=` bounds are required.

`goodness_result` contains the corrected Cramer-von Mises statistic, corrected
Anderson-Darling statistic, KS statistic and p-value, MLE, standard errors,
log likelihood, AIC, AICc, BIC, HQIC, objective value, convergence code, and
PDF/CDF validation information.

## Numerical notes

See `PORTING_NOTES.md` for deliberate corrections and differences from the R
implementation. The original package is preserved under `upstream/` for
license and source auditing.
