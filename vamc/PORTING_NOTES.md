# Porting notes

## Representation

R data frames are replaced by derived types:

- `date_type`
- `yield_curve_type`
- `policy_type`
- `portfolio_type`
- `mortality_table_type`
- result and status types

Scenario arrays retain the upstream logical ordering
`(scenario, time step, index or fund)`.

## Exact numerical parity

The swap bootstrap follows the upstream secant iteration and date conventions.
The published `buildCurve` test vector is reproduced to better than `1e-9` in the
checked test suite.

Random streams are deterministic after `seed_rng` or the optional `seed` arguments,
but they are not R's Mersenne-Twister stream and therefore are not bit-identical to R.

## Upstream behaviors preserved by default

Several source behaviors appear accidental but are retained to make the translation
auditable:

1. Inside the range of a user holiday vector, `isBusinessDay` returns true for dates
   listed as holidays and false for unlisted dates. Pass
   `source_compatible_holidays=.false.` for conventional holiday exclusion.
2. `genIndexScen` uses row sums of the squared upper Cholesky factor for the drift
   correction. Pass `source_compatible_drift=.false.` to use covariance diagonals.
3. Anniversary-benefit policies reset issue/maturity dates every projection step,
   and the elapsed number of days is assigned as a number of months. Set
   `source_compatible_ab_renewal=.false.` to renew only at a maturity event using
   the original term in months.
4. Maturity-benefit routines apply `max` to the entire preallocated account-value
   vector, including future zeros. Set `source_compatible_maturity_vector=.false.`
   to use the current account value only.
5. The upstream maturity test is evaluated before the monthly date increment while
   the horizon contains only the month difference, so a maturity exactly at the end
   of the horizon is normally not paid. Set `source_compatible_timing=.false.` for
   end-of-period maturity recognition.
6. Historical aging uses zero-based month offsets as R indices. Set
   `source_compatible_aging_index=.false.` for the natural one-based month mapping.

The upstream withdrawal routines may reference `dWA` before its first assignment.
The Fortran implementation initializes it to zero so results remain deterministic
and memory-safe.

## Mortality

The uniform-distribution-of-deaths formula is preserved. The Fortran implementation
extends mortality beyond the supplied maximum age with probability one, preventing
the out-of-bounds access possible in the R loop if survival remains above the cutoff.

## Scope omissions

- R data-frame conversion and printing
- `.rda` deserialization
- CSV I/O in the private `readMortTable` helper
- R warning/error class behavior

All computational valuation, curve, scenario, mortality, policy-generation, and aging
routines are implemented.
