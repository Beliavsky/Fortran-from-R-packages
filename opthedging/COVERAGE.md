# Computational coverage

## Package summary

- Original package: `OptHedging` 1.0
- Original author: Bruno Remillard
- Declared license: GPL version 2 or later
- Original numerical implementation: R wrappers plus one C source file

## Routine mapping

| Original routine | Fortran implementation | Status |
|---|---|---|
| `interpol1d` | `interpolation1d` | Complete |
| C `interpolation1d` | `linear_interpolate_uniform` | Complete |
| `hedging.iid` | `hedging_iid` | Complete |
| C `HedgingIID` | `hedging_iid` | Complete |
| C `Cn` | `call_payoff`, `put_payoff` | Complete and generalized |
| C `xbar` | `mean_value` | Complete |
| C `x2bar` | `mean_square` | Complete |

## Additional typed interfaces

The original R function returned an untyped list. The Fortran port returns a
`hedging_result` with the same numerical arrays and adds these methods:

- `option_value_at`
- `auxiliary_at`
- `shares_at`
- `initial_hedge_at`

The supplied `rho` and `phi1` fields correspond directly to the original R
output.

## Excluded infrastructure

No numerical functionality is excluded. The following R-specific mechanisms
are unnecessary in the compiled library:

- `.C` foreign-function calls.
- R memory allocation and list construction.
- R matrix reshaping.
- R documentation and namespace machinery at runtime.
- Plotting in the package example.

The original R source and documentation remain available under `original`.
