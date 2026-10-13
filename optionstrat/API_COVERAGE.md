# API coverage

Upstream optionstrat 1.4.1 has 30 explicit exports. Six plotting-only exports are excluded:
`plotbearcall`, `plotbearput`, `plotbullcall`, `plotbullput`, `plotdv`, `plotvertical`.
The remaining 24 exported computational functions are all substantially mapped.
Internal helpers and newly added corrected routines do not increase the denominator or numerator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `calldelta` | `optionstrat` | `calldelta` | substantial |
| `calleval` | `optionstrat` | `calleval` | substantial |
| `callgreek` | `optionstrat` | `callgreek` | substantial |
| `callpremium` | `optionstrat` | `callpremium` | substantial |
| `callrho` | `optionstrat` | `callrho` | substantial |
| `calltheta` | `optionstrat` | `calltheta` | substantial |
| `dv` | `optionstrat` | `dv` | substantial |
| `iv.calc` | `optionstrat` | `iv_calc` | substantial |
| `lambda` | `optionstrat` | `lambda` | substantial |
| `opteval` | `optionstrat` | `opteval` | substantial |
| `optiongamma` | `optionstrat` | `optiongamma` | substantial |
| `optionvega` | `optionstrat` | `optionvega` | substantial |
| `prob.above` | `optionstrat` | `prob_above`, `prob_above_quantile` | substantial |
| `prob.below` | `optionstrat` | `prob_below`, `prob_below_quantile` | substantial |
| `prob.btwn` | `optionstrat` | `prob_btwn`, `prob_btwn_quantile` | substantial |
| `putdelta` | `optionstrat` | `putdelta` | substantial |
| `puteval` | `optionstrat` | `puteval` | substantial |
| `putgreek` | `optionstrat` | `putgreek` | substantial |
| `putpremium` | `optionstrat` | `putpremium` | substantial |
| `putrho` | `optionstrat` | `putrho` | substantial |
| `puttheta` | `optionstrat` | `puttheta` | substantial |
| `r.cont` | `optionstrat` | `r_cont` | substantial |
| `tdiff` | `optionstrat` | `tdiff` | substantial |
| `vertical` | `optionstrat` | `vertical` | substantial |

## Mapping details

- `calldelta`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `calleval`: All six numerical fields are translated to option_result. R data.frame construction and recycling are omitted.
- `callgreek`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `callpremium`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `callrho`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `calltheta`: Preserves the upstream full 1 / t multiplier; this differs from a standard daily theta. calltheta_standard is a separately documented extension.
- `dv`: All ten numeric source fields are translated. Derived type replaces transposed data.frame; two-decimal rounding is omitted and invalid strike ordering returns NaN.
- `iv.calc`: Preserves the 0:0.001:5 nearest grid and its forward-ATM first-grid-point defect. Invalid inputs return NaN. The continuous implied_volatility extension is not counted separately.
- `lambda`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `opteval`: Both six-field option evaluations are translated to dual_option_result. R data.frame and character-label construction are omitted.
- `optiongamma`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `optionvega`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `prob.above`: Normal arithmetic-return probability and upstream whole-quantile negation are preserved. Quantile output is a price scalar; R data.frame and percent-change presentation are omitted. Corrected quantile extension is separate.
- `prob.below`: Normal arithmetic-return CDF and quantile branches are translated. Quantile output is a price scalar; R data.frame and percent-change presentation are omitted.
- `prob.btwn`: Normal arithmetic-return probability and the upstream sign-mirrored quantile range are translated into a price_range result. Nonzero-mean source convention is preserved; corrected extension is separate.
- `putdelta`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `puteval`: All six numerical fields are translated to option_result. R data.frame construction and recycling are omitted.
- `putgreek`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `putpremium`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `putrho`: Numerical operation translated as a pure elemental function. Decimal inputs, source units, and tiny-tenor convention are retained. Explicit zero-volatility limits and NaN validation replace R warnings, recycling, and undefined arithmetic.
- `puttheta`: Preserves the upstream full 1 / t multiplier; this differs from a standard daily theta. puttheta_standard is a separately documented extension.
- `r.cont`: Computational rate conversion uses the equivalent stable n * log1p(r / n); nonpositive frequency and nonpositive compounding base return NaN.
- `tdiff`: Date difference and 365-day years are translated for valid Gregorian ISO dates YYYY-MM-DD. General R datetime parsing, time zones, subday values, and fallback period coercion are omitted.
- `vertical`: All four call/put debit/credit branches and source probabilities are translated. Derived type replaces transposed data.frame; two-decimal rounding is omitted and invalid equal strikes return NaN.
