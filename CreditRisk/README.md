# CreditRisk

Modern free-form Fortran translation of the computational exports of CRAN
**CreditRisk 0.1.7**, by Alessandro Cimarelli and Nicolò Manca. It includes
structural survival models, geometric-Brownian firm-value simulation, CDS
valuation, interpolation, and four bounded calibrators.

## Important compatibility notes

This is a translation of source numerics, not a validation of financial models.
Two material upstream formula quirks are deliberately preserved:

- `Merton` survival uses `log(V0/L) + r - sigma^2*t/2` in its numerator.
  The usual textbook expression has `r*t`, not `r`. They agree when `t=1`
  (or `r=0`), but generally give different probabilities. Equity/debt still use
  the upstream generalized Black–Scholes formula.
- With `simplified=.false.`, `cds` subtracts accrued premium on default from
  the scheduled premium annuity. The conventional accrued-premium leg adds
  that integral. This package preserves the upstream subtraction, including
  the possibility of nonpositive premium annuities for extreme hazards.
  `simplified=.true.` instead uses discrete endpoint legs. The upstream R
  prose reverses these labels in one place; the implemented code is followed.

One upstream reporting bug is corrected: `calibrate_black_cox` returns actual
fitted `L` and `K`. R `calibrate.BlackCox` returns their starting values, although
its reported volatility, barrier-discount rate and error come from the fit.
The tests record both the unmodified R return and the actual R optimizer output.

Further intentional differences:

- Typed results replace R data frames/lists. No recycling, formula processing,
  plotting, data download, or R-specific interface is provided. The bundled
  `cdsdata` dataset is excluded from computational coverage and from the port.
- Invalid domains/shapes return `credit_invalid`; scalar option errors return
  quiet NaN. Curve outputs must be finite to have `credit_ok`. Check status
  before using results. A finite but mathematically unusual result is not a
  guarantee of economic plausibility. Overflow, underflowed survival, and
  degenerate leg ratios can make status invalid.
- Interpolation requires sorted, unique finite knots. R `approx` sorting,
  missing-value removal and tie averaging are omitted. Constant extrapolation
  is retained, and singleton input grids are additionally supported.
- Expiry/zero-volatility options, zero-variance AT1P, single-date structural
  curves, and zero combined hazard/discount rate are handled by explicit limits.
  Some corresponding R calls divide by zero or fail because of `2:n` loops.
- `Merton.sim` is available with supplied standard-normal shocks (for exact
  paired-path comparisons), or a local seeded Park–Miller/Box–Muller stream.
  The latter does not reproduce R's RNG or change global RNG state.
- Calibration uses the existing sibling L-BFGS-B 3.0 package, rather than R
  `optim`'s implementation. Defaults retain R starts, zero lower bounds,
  memory 5, `factr=1e7`, `pgtol=0`, and finite-difference step 0.001.
  This dependency scales that step by `max(1,abs(parameter))`, unlike R's
  absolute default step. Invalid model trials receive a finite penalty of
  `max(1,100*initial_MSE)`; R instead can stop with a nonfinite-objective
  error. The moderate objective-scaled penalty avoids artificial line-search
  stagnation caused by a huge penalty. Results include convergence diagnostics,
  initial/final errors, RMS and maximum quote residuals, and a warning when
  the maximum residual exceeds one basis point. Parameter-identifiability and
  local minima remain the caller's responsibility; optimizer trajectories need
  not match R.
- `calibrate_cds` interpolates *rates*, as R does. `cds2` instead computes a
  ratio of interpolated legs. These are different for off-grid quote dates.
- Upstream `calibrate.sbtv` silently ignores its ellipsis. The translated
  routine therefore exposes the same default CDS settings only.

See [API.md](API.md) for argument order, columns, parameter order and domains.

## Build and run

Place this directory at the root of `Fortran-from-R-packages`, alongside the
existing `rfortran-core` and `lbfgsb3` directories. They are reused through
sibling FPM path dependencies and are not included here.

```sh
fpm build
fpm test
fpm run --example --all
```

The default build needs FPM and gfortran, including on Windows. No separately
installed BLAS/LAPACK, C compiler, R runtime, or downloaded dependency is
required for these commands. Do not use `-ffast-math`, `-Ofast`,
`-ffinite-math-only`, or equivalent flags: invalid-result handling depends on
IEEE NaNs and infinities. Maintained real declarations and constants use the
shared `dp` imported from `r_kinds`. Public model procedures are pure; scalar
CDF/option routines are pure elemental. Calibration itself is impure because
it calls the shared optimizer.

```fortran
use creditrisk, only: dp, credit_curve_t, at1p, cds, credit_ok
implicit none
type(credit_curve_t) :: model, swap
real(dp) :: t(2), r(2)
t = [1.0_dp, 2.0_dp]
r = [0.02_dp, 0.025_dp]
model = at1p(1.0_dp, 0.6_dp, 0.4_dp, [0.2_dp, 0.23_dp], r, t)
if (model%status /= credit_ok) error stop 'Invalid structural curve'
swap = cds(t, model%values(:,5), r)
if (swap%status /= credit_ok) error stop 'Invalid CDS curve'
print *, swap%values(:,5)
```

## Testing and provenance

`reference/generate_reference.R` runs the actual installed CRAN package and
records numerical outputs plus its R session information. It uses a read-only
`trace` on `stats::optim` to observe optimized BlackCox parameters hidden by
the reporting bug; it does not modify the objective or optimizer. Numeric
fixtures are embedded in `test/reference_values.f90`, so normal tests need no
R installation. Simulation comparisons use innovations from the real R run.

Tests cover all 13 exports, including all four nontrivial calibrations, R
parameter/rate/objective agreement, reconstruction of fitted objectives,
put–call parity, equity plus debt, reflection-formula survival, scenario
mixtures, hazard-to-survival inversion, analytic CDS integrals, boundary
limits, deterministic RNG, interpolation and invalid shapes/domains. Separate
analytic assertions demonstrate the preserved Merton and CDS discrepancies
rather than mistaking agreement with R for textbook validation. Numerical
calibration comparisons allow optimizer-dependent differences: maximum
absolute errors of `1e-5` for parameters, `2e-7` for fitted rates and `1e-10`
for objective values on the original recorded cases. A further independent
SBTV case exercises invalid-trial recovery: its fit is assessed by objective
quality and rates rather than matching nonunique parameters. Regression tests
also require a residual warning for a converged but inadequate fit and reject
unrepresentable objectives before entering the optimizer.

## Translation coverage

Package status: **substantial**. **13 of 13 (100%)** exported computational R
functions have a mapping. This fraction measures mapped computational
functions, not complete R compatibility. The dataset and presentation-only R
code are outside the denominator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `BlackCox` | `creditrisk` | `black_cox` | complete |
| `Merton` | `creditrisk` | `merton` | complete |
| `Merton.sim` | `creditrisk` | `merton_sim`, `merton_sim_from_shocks` | substantial |
| `at1p` | `creditrisk` | `at1p` | complete |
| `calibrate.BlackCox` | `creditrisk_calibration` | `calibrate_black_cox` | substantial |
| `calibrate.at1p` | `creditrisk_calibration` | `calibrate_at1p` | substantial |
| `calibrate.cds` | `creditrisk_calibration` | `calibrate_cds` | substantial |
| `calibrate.sbtv` | `creditrisk_calibration` | `calibrate_sbtv` | substantial |
| `cds` | `creditrisk` | `cds` | complete |
| `cds2` | `creditrisk` | `cds2` | substantial |
| `cum_normal_density` | `creditrisk` | `cum_normal_density` | complete |
| `generalized_black_scholes` | `creditrisk` | `generalized_black_scholes` | complete |
| `sbtv` | `creditrisk` | `sbtv` | complete |

The exact per-function source mappings are in `fpm.toml`; [API.md](API.md)
summarizes callable interfaces and their limitations.

## License

The translated package and upstream CreditRisk source are MIT licensed; see
[LICENSE](LICENSE), [NOTICE.md](NOTICE.md), and the original upstream license
metadata. The separately reused `rfortran-core` dependency is MIT licensed.
The separately reused `lbfgsb3` dependency declares **GPL-2.0-only**: distributing
a combined linked program must comply with that dependency's GPL terms. The
MIT label for this package does not relicense that dependency or the combined
binary. No dependency source is copied or vendored.
