# API map

| Upstream R routine or object | Modern Fortran equivalent | Notes |
|---|---|---|
| `fWilson` | `wilson_function` | Pure elemental scalar function; matrix construction is separate. |
| `fCreateKernelMatrix` | `create_kernel_matrix` | Constructs the symmetric Wilson matrix for calibration times. |
| `fFitKernelWeights` | `fit_kernel_weights` | Uses self-contained pivoted Gaussian elimination. |
| `fFitSmithWilsonYieldCurve` | `fit_smith_wilson_curve` | Returns `type(smith_wilson_curve)`. |
| `fFitSmithWilsonYieldCurveToInstruments` | `fit_smith_wilson_curve_to_instruments` | Accepts an array of `type(market_instrument)`. |
| `fCreateTimeVector` | `create_time_vector` | Sorted unique union of all payment times. |
| `fCreateCashflowMatrix` | `create_cashflow_matrix` | Instruments in rows and payment times in columns. |
| `fGetTimesLibor`, `fGetTimesSwap`, `fGetTimesBond` | `get_instrument_times` | Dispatches on `instrument_type`. |
| `fGetCashflowsLibor`, `fGetCashflowsSwap`, `fGetCashflowsBond` | `get_instrument_cashflows` | Returns both schedule times and cashflows. |
| R instrument data frame | `type(market_instrument)` | Integer type constants are `sw_libor`, `sw_swap`, and `sw_bond`. |
| R curve list element `P` | `curve%discount(...)` | Scalar and vector overloads. |
| R curve list element `K` | `curve%compound_kernel(...)` | Returns instrument-by-query-time values for vector input. |
| R curve list element `xi` | `curve%xi` | Public allocatable component. |
| Repricing through `P` and cashflows | `curve%repriced_values()` | Returns one fitted value per calibration instrument. |
| Plot/lines/points methods | Omitted | Plotting is outside the computational port. |
| Generic internal `fFitYieldCurve` | Smith-Wilson fit pipeline plus `fit_kernel_weights` | The R closure-based generic was not retained as a stored procedure closure. |

## Error handling

R exceptions and warnings are represented by an integer `info` and optional
character `message` output. Public status constants are:

- `sw_success`
- `sw_invalid_argument`
- `sw_dimension_error`
- `sw_singular_system`
- `sw_unknown_instrument`

As in the R source, a negative UFR does not prevent fitting. The Fortran
routine succeeds and returns a warning in `message`.
