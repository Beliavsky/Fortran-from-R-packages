# credule

Modern, free-form Fortran translation of CRAN **credule 0.1.4** by Bertrand Le
Nezet. Prices credit default swaps and bootstraps piecewise-constant hazard-rate
curves from par spreads. All numerical library procedures are pure and use the
`dp` kind imported from the sibling `rfortran-core` package. There is no mutable
global solver state, R runtime dependency, or external BLAS/LAPACK dependency.

This is the upstream discretized model, not an ISDA standard-model implementation.
Times are year fractions; zero rates, spreads and hazard rates are annual decimals.
For example, 100 basis points is `0.01_dp`. No calendar, settlement, day-count,
upfront-payment, or business-day conventions are inferred.

## Build and run

Extract this directory as `credule/` beside an existing `rfortran-core/` at the
root of [Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages).
The only dependency is declared as `rfortran-core = { path = "../rfortran-core" }`.
Dependency source is not bundled here. No shared package provides a compatible
CDS or bracketed Brent API; the small upstream Brent component is translated
separately with its original license. No `rfortran-compat` is used.

```sh
cd credule
fpm build
fpm test
fpm run --example --all
```

FPM and gfortran are sufficient, including on Windows without separately installed
BLAS or LAPACK. Validation was performed on Linux with gfortran; Windows execution
has not been claimed. Do not compile with `-ffast-math`, `-Ofast`, or
`-ffinite-math-only`: validation and error contracts depend on IEEE NaNs/infinities.

## Translation coverage

Package status: **substantial**. **2 of 2 (100%)** exported computational R
functions are mapped. This fraction counts mapped functions, not complete R
compatibility. Both exported numerical operations are implemented; data-frame
construction, `.C` marshaling, coercion, warning/stop presentation, and vector-valued
R defaults are replaced by typed Fortran arrays, scalar options, and status codes.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `priceCDS` | `credule` | `price_cds` | substantial |
| `bootstrapCDS` | `credule` | `bootstrap_cds` | substantial |

The exact mappings and notes are also in `fpm.toml`; see `API_COVERAGE.md` for the
scope of internal helpers and deliberately changed edge behavior.

## API

```fortran
use credule, only : dp, price_cds, bootstrap_cds, credule_success
! Caller-allocated arrays have matching extents; every call returns status.
call price_cds(yield_tenors, yield_rates, credit_tenors, credit_sp, &
               maturities, recovery, spreads, status)
call bootstrap_cds(yield_tenors, yield_rates, maturities, spreads, &
                   recovery, survival, hazards, status)
```

Both accept optional scalar arguments `premium_frequency` (default 4),
`default_frequency` (default 12), and `accrued_premium` (default `.true.`).
Any positive integer frequency is supported, including the R interface's common
1, 2, 4 and 12 choices. Outputs must be allocated by the caller; procedures do
not allocate, resize, print, or stop. Failed calls fill **all supplied output
array elements with quiet NaNs**, so partial calibration cannot look successful.

- Yield tenors: finite, positive and strictly increasing. Rates: matching finite
  continuously compounded zero rates, which may be negative.
- Credit tenors: finite, positive and strictly increasing. Survival probabilities:
  matching, finite, nonincreasing values in `(0,1]`; zero is rejected because the
  log-hazard interpolation would be undefined.
- Pricing maturities: finite, positive and in any order. Bootstrap maturities must
  be strictly increasing. Each must span at least one premium and default interval.
  An empty maturity/output vector is accepted if the other inputs are valid.
- Recovery: finite in `[0,1]` for pricing and `[0,1)` for calibration. Full recovery
  has zero protection value and cannot identify a calibrated hazard rate.
- Bootstrap spreads: finite nonnegative decimals. Hazards are annual **interval**
  hazards, not cumulative average hazards. The retained calibration bracket is
  `[0,30]` per year. Nonmonotone or other inconsistent spread curves may fail.

Public supporting procedures:

- `discount_factor(tenors, rates, time)` returns the interpolated discount factor.
- `survival_probability(tenors, probabilities, time)` returns interpolated survival.
- `cds_legs(yield_tenors, yield_rates, credit_tenors, credit_sp, maturity, recovery,`
  `premium_annuity, default_leg, status, premium_frequency, default_frequency, accrued_premium)`
  returns the unit-spread premium PV and unit-notional protection PV.

The scalar interpolation helpers return quiet NaN on invalid inputs. Time zero
returns 1 when a valid nonempty curve is supplied. `cds_legs` validates inputs and
shares exactly the conventions of `price_cds`. All public procedures have explicit,
one-argument-per-line, trailing FORD documentation in the source.

| Status constant | Value | Meaning |
| --- | --- | --- |
| `credule_success` | 0 | All outputs calculated |
| `credule_invalid_input` | 1 | Invalid shape, curve, contract, frequency or spread |
| `credule_no_bracket` | 2 | No sign-changing hazard bracket on `[0,30]` |
| `credule_no_convergence` | 3 | Brent exceeded 256 iterations |
| `credule_numerical_failure` | 4 | Nonfinite/degenerate legs or survival underflow |

## Numerical conventions and compatibility

1. Zero rates are interpolated linearly; the first/last zero rate is held constant
   outside the yield knots. Discounting is `exp(-rate(time) * time)`.
2. Survival is interpolated by a constant hazard on each interval, starting from
   survival 1 at time zero. The standalone survival helper extends the final
   interval hazard beyond the last knot (or the origin-to-knot hazard for one knot).
3. Premium cash flows are summed at right endpoints on a regular grid. Default
   payments are likewise discounted at the **end** of each default interval.
   Accrued premiums, when enabled, add one half of a premium period times the
   discounted survival decrement. These are the exact upstream approximations.
4. Each leg uses `int(maturity * frequency)` complete intervals. Fractional final
   intervals are omitted, even when the two leg grids end at different dates.
   Floating-point multiplication and truncation are not replaced by rounding.
5. An important upstream peculiarity is deliberately retained: when **pricing past
   the last credit knot**, the unknown tail hazard is zero. Each leg independently
   anchors at `floor(last_credit_tenor * its_frequency) / its_frequency`. Survival
   remains constant after that anchor. This differs from the standalone survival
   helper's terminal-hazard extrapolation. Use a supplied credit knot at or after
   the requested maturity if continued hazard exposure is intended.
6. Bootstrap proceeds sequentially. Known points enter interpolation; each new
   hazard uses that same **leg-specific grid anchor** after the previous credit
   knot. The stored survival, however, advances over the actual tenor difference,
   exactly as upstream. For non-grid credit knots, this can prevent perfect
   `price_cds`/`bootstrap_cds` round trips. Grid-aligned knots are recommended;
   the reference suite checks non-grid R parity instead of asserting false inversion.
7. Pricing uses protection PV divided by unit-spread premium PV. This replaces
   upstream's unnecessary linear spread root solve and removes its artificial
   `[0,30]` spread bracket; ordinary valid inputs agree to floating-point tolerance.
   Calibration retains Brent's original machine-epsilon convergence criterion,
   with explicit bracket checks and a bounded iteration count.
8. Exact endpoint hazard solutions, notably a zero-spread/default-free curve, are
   valid here. Upstream flags endpoints as failures. Invalid inputs and unbracketed
   roots receive explicit statuses rather than the original undefined access,
   warning, or potentially misleading root. Extreme numerical underflow/overflow
   is reported rather than presented as a valid curve.

## Verification

`fpm test` checks 102 genuine R pricing cases and 78 genuine R bootstrap cases,
with frozen expectations generated by the installed, unmodified CRAN credule
0.1.4 package. `test/reference/generate_reference.R` regenerates them; R is needed
only for regeneration. References cover zero/negative/positive/shaped yields,
single- and multi-knot credit curves, low and distressed hazard rates, recoveries
0/0.4/0.8, annual through monthly premiums, 4/12/52/365 default grids, both accrual
choices, extrapolation, unsorted pricing maturities and fractional-grid knots.

Independent checks use closed-form geometric-series legs under constant rates,
24 varied-frequency/accrual pricing-calibration inversions, repricing and
no-arbitrage identities, explicit interpolation/extrapolation values, endpoint
roots, malformed inputs, nonfinite values, shape mismatch, and transactional error
outputs. R-reference tolerances are `2e-12 * max(1,abs(reference))` for prices and
`3e-12 * max(1,abs(reference))` for calibration. The example verifies its own market
repricing residual. Tests are deterministic and have no downloads or RNG calls.

## Licensing and provenance

The original package is **MIT**, copyright 2015 Bertrand Le Nezet. Its bundled
Brent solver is separately **GNU LGPL** and is credited to Richard Brent and
John Burkardt, with the 11 February 2013 correction acknowledgement to Thomas
Secretin. That component is retained as `src/credule_brent.f90`; this distribution
elects LGPL 2.1 or later from the upstream version-unspecified LGPL grant.
The combined manifest therefore says `MIT AND LGPL-2.1-or-later`.

Read `LICENSE`, `NOTICE`, and `licenses/`. Unmodified upstream R/C source,
metadata and license files are retained under `upstream/` for attribution and
verification only; FPM compiles only the maintained `.f90` sources. No third-party
Fortran dependency has been copied into this package.
