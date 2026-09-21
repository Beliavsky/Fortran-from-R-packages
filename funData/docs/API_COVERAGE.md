# API coverage

## Coverage basis

The coverage denominator is the set of distinct exported computational R functions in the upstream NAMESPACE. Plotting, printing, summaries, names/accessors/setters, constructors whose role is S4 object assembly, data-frame presentation conversion, and `fd2funData`/`funData2fd` wrappers that delegate the substantive conversion to the external `fda` package are excluded.

The resulting basis contains 20 functions. All 20 have meaningful numerical Fortran mappings, so mapped-function coverage is 100.0%. The package status remains `substantial` because mapped numerical operations do not imply complete R/S4 interface parity.

## Native representations

- `fun_data`: regular functional observations with support grids, an R-order dimension vector `[N,M1,...,Md]`, and flat column-major `real(dp)` data.
- `irreg_fun_data`: a collection of one-dimensional irregular curves.
- `multi_fun_data`: a collection of regular functional components with common observation count.
- `basis_spec`: component specification used by multivariate simulation.

## Material compatibility differences

- S4 dispatch, object names, formula-like/non-standard evaluation, warning/error text, and R attributes are not reproduced.
- Regular extraction uses resolved integer support indices; R callers normally provide support values that R resolves through membership tests.
- Simulation/sparsification/error routines use Fortran RNGs rather than R RNG state and are statistically compatible rather than stream-identical.
- `approxNA` directly performs linear interpolation of interior NaN runs instead of dispatching to `zoo::na.approx`.
- Integration is intentionally limited to support dimension 1-3, matching upstream's implemented numerical scope.
- One-point irregular full-domain extrapolation follows the documented constant-function-value rule rather than the apparent `rep(x,3)` typo in upstream `.extrapolate`.
- `funData2fd` and `fd2funData` are excluded because their substantive behavior is conversion to/from the external R package `fda`, not an independent funData numerical algorithm.

The authoritative per-function mapping is in `fpm.toml`; the README table is generated from the same 20-entry mapping used to write that manifest.
