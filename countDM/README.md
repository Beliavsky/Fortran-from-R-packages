# countDM-fortran

Modern Fortran 2018/FPM translation of the computational code in the R package
`countDM` 0.1.0.

The upstream package estimates count-data models based on the Bell, Borel,
Poisson, Bell--Touchard, zero-inflated, and zero/one-inflated distributions.
This translation preserves the numerical/statistical functionality while
removing dependencies on R, `maxLik`, `lamW`, `numbers`, and `miscTools`.

## Main API

Use the umbrella module:

```fortran
use countdm
```

Distribution functions:

- `tp` / `touchard_polynomial`
- `bell_number`
- `dbell`, `dborel`, `dpoisson_count`
- `dbellt`, `pbellt`, `qbellt`, `rbellt`
- `dzibellt`, `pzibellt`, `qzibellt`, `rzibellt`
- `dzip`, `dzibell`, `dzoip`, `dzoibell`

Maximum-likelihood fitting:

- `bell_mle` / `bell_mle_closed`
- `mle_bell`
- `mle_borel`
- `mle_poisson`
- `mle_bt`
- `mle_zip`
- `mle_zibell`
- `mle_zibellt`
- `mle_zoip`
- `mle_zoibell`

`mle_result_t` contains the estimates, standard errors, covariance matrix,
maximized log likelihood, AIC, convergence flag, and iteration count.

Bundled data are available with either `data_criminal()` / `data_sbirth()` or
the allocation subroutines `get_data_criminal` / `get_data_sbirth`.

## Numerical implementation

- Touchard polynomials use the Stirling-number recurrence.
- Bell numbers are `T_n(1)`.
- The closed-form Bell MLE uses the principal Lambert-W solution
  `theta = W_0(mean(x))`.
- General MLEs use a self-contained BFGS optimizer, central numerical
  derivatives, and a numerical observed-information Hessian.
- Positive parameters are optimized on the log scale.
- Unit-interval parameters use the logit scale.
- Zero/one-inflation weights use a three-category softmax so
  `alpha >= 0`, `beta >= 0`, and `alpha + beta < 1` automatically.

## Upstream corrections

Two clear implementation defects in the original R source are corrected while
retaining the intended documented model:

1. `mle_zibellt()` defines its likelihood parameter order as `(pi, lambda,
   theta)` but supplies starting values as `(lambda, theta, pi)`. The Fortran
   routine uses the documented function order `(lambda, theta, pi)` consistently.
2. `qbellt()` in R ignores `lower.tail` and, when `log.p=TRUE`, logs the returned
   integer quantile instead of interpreting the input probability on the log
   scale. The Fortran quantile routines implement conventional probability
   semantics.

The upstream zero/one-inflated Bell likelihood is otherwise preserved as
written: masses `alpha` and `beta` are assigned directly at 0 and 1 and the
remaining Bell law is conditioned on values greater than 1.

## Validation

The source is validated with GNU Fortran using:

```text
gfortran -std=f2018 -Werror=implicit-interface -Werror=trampolines \
         -fcheck=all -O0
```

Tests cover known Touchard/Bell values, Bell--Touchard normalization and
quantiles, zero-inflated identities, bundled datasets, closed-form MLEs, all
numerical MLE routines, standard errors/AIC, and independent SciPy likelihood
optima for the inflated models.

## License and provenance

Upstream `countDM` is GPL (>= 2). This translation is distributed under
GPL-2.0-or-later. The original `DESCRIPTION`, `NAMESPACE`, R source, and Rd
manual files are retained under `upstream/` for provenance.
