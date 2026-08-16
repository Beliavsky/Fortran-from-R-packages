# DiscreteDists-fortran

Modern Fortran 2018/FPM computational port of the R package `DiscreteDists`
(version 1.1.2).

## Scope

The port provides the 14 discrete distribution families from the upstream
package:

- BerG
- COMPO
- COMPO2
- DBH
- DGEII
- DIKUM
- DLD
- DMOLBE
- DPERKS
- DsPA
- GGEO
- HYPERPO
- HYPERPO2
- POISXL

Each family has density/PMF, CDF, quantile and random-generation routines.  The
GAMLSS family constructors are represented by `discrete_family_t`, which
provides links, score/curvature contributions, deviance increments,
initialization, moments, validation and randomized quantile residuals without
requiring R or GAMLSS at runtime.

The exported likelihood/start-value helpers, `F11`, `AR`, `add`, `stopping`,
`obtaining_lambda`, `mean_var_hp`, `mean_var_hp2`, and `simulate_hp` are also
ported.  The only NAMESPACE export omitted is `plot_discrete_cdf`, because it
is plotting-only.

## Dependency

`COMPoissonReg-fortran v0.1.0` is vendored as an FPM dependency and supplies
the CMP normalizer/distribution machinery used by COMPO and COMPO2.

The upstream R package imports `nleqslv`, but the current source only mentions
it in commented-out initialization code, so it is not a runtime dependency of
this port.  The supplied GAMLSS translations were useful for checking family
semantics but are likewise not runtime dependencies.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

## Minimal example

```fortran
use discretedists
implicit none
real(dp) :: p
integer :: q

p = pcompo(4.0_dp, 3.0_dp, 0.8_dp)
q = qcompo(0.75_dp, 3.0_dp, 0.8_dp)
```

## Validation

The test suite covers distribution/CDF identities and all 14 quantile APIs,
starting-value estimators, GAMLSS-compatible family methods, and random-number
checks.  It is run both with bounds/runtime checking and with optimization.

See `PORTING_NOTES.md` for compatibility details and `API_MAP.md` for the full
R-to-Fortran mapping.
