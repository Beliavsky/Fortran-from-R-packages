# Porting notes

## Array orientation

Fortran curve arrays are ordinary one-dimensional arrays. Historical yield input to `curve_calculation` is shaped `number_of_terms x number_of_dates`, so each column is one observation date.

## Dates and calendars

`qbc_date` avoids compiler- and locale-dependent date parsing. Month addition follows the same end-of-month clipping convention used by `lubridate`'s `%m+%` and `%m-%` operators.

Business-day adjustment recognizes weekends. The upstream package delegates named BOG, LDN, NY, and NYLDN holiday calendars to `quantdates`; those holiday datasets are not included. Applications needing exact settlement calendars should adjust the returned nominal dates with their own holiday service or extend `qbc_dates`.

## Coupon schedules

The four schedule constants are:

- `qbc_schedule_short_first`
- `qbc_schedule_long_first`
- `qbc_schedule_short_last`
- `qbc_schedule_long_last`

A zero payment frequency represents a single maturity cash flow.

## Curves

`qbc_curve%approximation` is `1` for piecewise constant and `2` for piecewise linear. Extrapolation is flat at the first and final rates, matching the practical behavior used by the upstream calibration code.

The direct-yield mode converts periodically compounded market rates to continuous rates before interpolation. Price calibration works directly in continuous zero rates, because this gives a stable unconstrained representation for discount factors.

## Optimization

`bootstrap_curve` and `basis_curve` use bounded Nelder-Mead. The objective is a weighted residual sum of squares. Bounds and starting values are explicit, and the returned `qbc_calibration_result` contains the objective, iteration count, status, and message.

## Swap valuation

Floating legs use the standard par-notional approximation plus the discounted spread leg. This is the same computational reduction used by the upstream package for variable legs. Cross-currency values are reported as local leg minus exchange-rate-adjusted foreign leg.

## Error handling

Most procedures accept an optional integer `status`. Public status constants are `qbc_success`, `qbc_invalid_argument`, `qbc_size_mismatch`, `qbc_no_convergence`, `qbc_singular`, and `qbc_infeasible`.
