# API

All floating-point calculations use `dp = kind(1.0d0)`. Public names are
available through `use fixedincome`.

## Types

- `term_t`: vector of term values with per-element day/month/year units.
- `daycount_t`: day-count specification and annual base.
- `spot_rate_t`: rate vector plus compounding, day count, and calendar.
- `spot_rate_curve_t`: rates, day terms, conventions, reference date, and an
  optional interpolation model.
- `forward_rate_t`: forward rates and their interval lengths.
- `interpolation_t`: interpolation method, parameters, and propagation flag.
- `fit_result_t`: fitted model, objective, iteration count, and status.

## Terms, dates, and day count

- `term(value, units)` and `term_units(values, units)`
- `parse_term("6 months")`
- `daycount("actual/365")`, `dib(dc)`
- `todays(dc,t)`, `tomonths(dc,t)`, `toyears(dc,t)`, `convert_term`
- `shift`, `term_difference`
- `gregorian_to_ordinal`, `ordinal_to_gregorian`, `offset_date`

The day-count parser accepts any `name/base` spelling with a positive numeric
base. Unit conversion follows the upstream package's fixed conversion map:
months contain `base/12` days.

## Compounding

- `compounding("simple"|"discrete"|"continuous")`
- `compounding_name(method)`
- `compound(method,time,rate)`
- `implied_rate(method,time,factor)`
- `discount(method,time,rate)`

Methods are also available as constants `COMPOUND_SIMPLE`,
`COMPOUND_DISCRETE`, and `COMPOUND_CONTINUOUS`.

## Rates and curves

- `spotrate`
- `spotratecurve`
- `forwardrate`
- `spot_compound`, `spot_discount`
- `curve_compound`, `curve_discount`, `forward_compound`
- `as_spotrate`, `as_forwardrate`, `as_spotratecurve`
- `forwardrate_from_curve`, `forwardrate_between`
- `curve_at_terms`, `insert_curve_points`
- `first`, `last`, `closest`, `maturities`

Constructors accept character convention names and return a status through an
optional integer argument. Curve terms are normalized to days and sorted.

## Interpolation

Constructors:

- `interp_flatforward`
- `interp_linear`
- `interp_loglinear`
- `interp_naturalspline`
- `interp_hermitespline`
- `interp_monotonespline`
- `interp_nelsonsiegel`
- `interp_nelsonsiegelsvensson`

Operations:

- `set_interpolation(curve,model)`
- `interpolation(curve)`
- `prepare_interpolation(model,curve)`
- `interpolate(curve,query_days)`
- `interpolation_error(curve)`
- `parameters(model)` and `interpolation_name(model)`

Flat-forward interpolation linearly interpolates the logarithm of accumulation
factors and then converts the result back to the curve's compounding regime.
Queries outside the data range return IEEE quiet NaNs and status
`FI_OUT_OF_RANGE`, matching R's default non-extrapolating behavior.

## Parametric curves

- elemental `nelson_siegel`
- elemental `nelson_siegel_svensson`
- `fit_interpolation(initial_model,curve)`

Fitting uses a bounded Nelder-Mead implementation with the same parameter
bounds as the upstream R code. The result status is `FI_OK` on convergence or
`FI_NO_CONVERGENCE` when the iteration limit is reached; the best parameters
found are returned in either case.

## Status codes

- `FI_OK`
- `FI_INVALID_ARGUMENT`
- `FI_SIZE_MISMATCH`
- `FI_OUT_OF_RANGE`
- `FI_NOT_CONFIGURED`
- `FI_NO_CONVERGENCE`
- `FI_UNSUPPORTED_CALENDAR`
