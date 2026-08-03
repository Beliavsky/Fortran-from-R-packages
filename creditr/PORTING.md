# Porting notes

## R-to-Fortran mapping

| R function | Fortran equivalent |
|---|---|
| `CDS` | `price_cds`, `cds` |
| `add_dates` | `add_dates` |
| `add_conventions` | `add_conventions` |
| `spread_to_upfront` | `spread_to_upfront` |
| `upfront_to_spread` | `upfront_to_spread` |
| `spread_DV01` | `spread_dv01` |
| `CS10` | `cs10` |
| `IR_DV01` | `price_cds` with market quotes; `ir_dv01` for a zero-rate bump |
| `rec_risk_01` | `rec_risk_01` |
| `spread_to_pd` | `spread_to_pd` |
| `pd_to_spread` | `pd_to_spread` |
| `implied_RR` | `implied_rr` |
| `PV01` | `pv01` |
| `get_rates`/`build_rates` | `read_rate_quotes_csv`; explicit quotes |
| `check_inputs` | typed arguments and status validation |
| `adj_next_bus_day` | internal following-business-day routine |
| `separate_YMD` | `date_t` fields |
| `call_ISDA` | native valuation and sensitivity routines |

## ISDA numerical translation

The port translates the algorithms needed by `creditr` rather than linking the original C objects. In particular it implements:

1. money-market discount factors;
2. par-swap bootstrapping with modified-following coupons;
3. logarithmic flat-forward discount interpolation;
4. quarterly premium schedules;
5. ISDA-style accrued-on-default integration;
6. protection-leg integration over discount-curve timelines;
7. clean-spread hazard calibration; and
8. clean/dirty upfront and spread conversion.

The representative 2014 USD curve reproduces the retained C implementation's first six money-market discount factors at machine precision and later swap nodes within a few parts in `1e12`.

## Intentional interface changes

- R data frames and S4 objects are replaced by typed derived types.
- Market curves are explicit inputs. The Fortran library does not download XML/HTML data.
- The bundled serialized R object `rates.RData` is retained only for provenance.
- Vectorized data-frame valuation becomes an ordinary Fortran loop over contracts.
- Invalid inputs return status codes; undefined scalar formulas return IEEE NaN.

## Calendars

`add_dates` preserves the package's behavior:

- value date advances by weekdays only, as explicitly noted in the R source;
- base date uses the package's finite USD and JPY holiday lists;
- coupon calendars use weekend adjustment when the R calls pass `NONE`.

The finite historic holiday tables are not a general-purpose modern holiday service. Applications valuing other dates should supply already appropriate base dates and curves or extend the calendar module.

## Interest-rate DV01

The R package adds one basis point to each input market quote and rebuilds the curve. `price_cds(contract, curve, status, quotes)` does the same. If quotes are unavailable, `price_cds` and standalone `ir_dv01` use a parallel continuous-zero-rate bump, which is close but not identical.

## Source behavior retained

- `spread_to_pd` uses the package's 360-day tenor convention.
- `implied_RR` returns a percentage, not a decimal.
- `PV01` retains the package formula based on principal and the spread/coupon difference.
- Sensitivities are forward differences, not centered derivatives.
