# ineq: modern Fortran numerical translation

Numerical translation of CRAN **ineq 0.2-13**, by Achim Zeileis with contributions
from Christian Kleiber. It computes inequality, concentration, poverty, empirical
and theoretical Lorenz curves, Mehran bounds, asymmetry, and majorization.
All computational procedures are pure; scalar theoretical curves without array
parameters are also elemental. All maintained Fortran uses `dp` from `r_kinds`.

## Build and run

Place this directory at the repository root next to `rfortran-core`:

```text
Fortran-from-R-packages/
  ineq/
  rfortran-core/
```

From `ineq`, with FPM and gfortran available:

```sh
fpm build
fpm test
fpm run --example --all
```

There are no separately installed BLAS/LAPACK requirements and no system-library
link flags. The sibling dependency is used through FPM's relative path support,
including on Windows. Verification was performed with Linux gfortran and FPM;
**a native Windows build has not been run**. IEEE NaN/infinity semantics are part
of the interface: do not use `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or
similar flags. Expected IEEE invalid/divide-by-zero values are tested; trapping
those exceptions is not a supported runtime configuration.

## Quick example

```fortran
program example
   use ineq_api, only : dp, gini, foster, lc, lorenz_curve, ineq_ok
   implicit none
   real(dp), parameter :: x(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 10.0_dp]
   type(lorenz_curve) :: curve
   print *, gini(x)
   print *, foster(x, 5.0_dp) ! Headcount below the poverty threshold.
   curve = lc(x)
   if (curve%status == ineq_ok) print *, curve%p, curve%l
end program example
```

Import all public numerical operations from `ineq_api`, or use `ineq_indices`
and `ineq_lorenz` directly. Public argument and result details, examples, and
edge policies are in [API.md](API.md).

## Translation coverage

Package status: **substantial**. **26 of 26 (100%)** exported computational R
functions have meaningful Fortran mappings. This fraction measures mapped
computational functions, not complete R compatibility or identical floating-point
rounding. All 26 mappings are substantial. `Pen` is a plotting-only export and is
excluded. Plot/lines S3 methods and the Ilocos dataset are excluded. `theorLc`
is counted and mapped as a numerical evaluator, although R closure creation is
not reproduced. No computational export is omitted from the denominator.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `ineq` | `ineq_indices` | `ineq` | substantial |
| `Gini` | `ineq_indices` | `gini` | substantial |
| `RS` | `ineq_indices` | `rs` | substantial |
| `Atkinson` | `ineq_indices` | `atkinson` | substantial |
| `Theil` | `ineq_indices` | `theil` | substantial |
| `Kolm` | `ineq_indices` | `kolm` | substantial |
| `var.coeff` | `ineq_indices` | `var_coeff` | substantial |
| `entropy` | `ineq_indices` | `entropy` | substantial |
| `conc` | `ineq_indices` | `conc` | substantial |
| `Herfindahl` | `ineq_indices` | `herfindahl` | substantial |
| `Rosenbluth` | `ineq_indices` | `rosenbluth` | substantial |
| `pov` | `ineq_indices` | `pov` | substantial |
| `Watts` | `ineq_indices` | `watts` | substantial |
| `Sen` | `ineq_indices` | `sen` | substantial |
| `SST` | `ineq_indices` | `sst` | substantial |
| `Foster` | `ineq_indices` | `foster` | substantial |
| `Lc` | `ineq_lorenz` | `lc` | substantial |
| `Lc.mehran` | `ineq_lorenz` | `lc_mehran` | substantial |
| `Lasym` | `ineq_lorenz` | `lasym` | substantial |
| `theorLc` | `ineq_lorenz` | `theor_lc` | substantial |
| `Lc.dagum` | `ineq_lorenz` | `lc_dagum` | substantial |
| `Lc.exp` | `ineq_lorenz` | `lc_exp` | substantial |
| `Lc.lognorm` | `ineq_lorenz` | `lc_lognorm` | substantial |
| `Lc.pareto` | `ineq_lorenz` | `lc_pareto` | substantial |
| `Lc.singh` | `ineq_lorenz` | `lc_singh` | substantial |
| `major` | `ineq_lorenz` | `major` | substantial |

The same inventory is recorded in `fpm.toml` and summarized in
[COVERAGE.md](COVERAGE.md). Formula/coercion handling, S3 classes, warning text,
R argument recycling, partial-name matching, graphics, and closure construction
are not implemented.

## Numerical and compatibility notes

- Missing inputs use IEEE NaN. Sample indices remove NaNs by default; pass
  `na_rm=.false.` to propagate them. Infinite empirical observations are rejected
  with NaN or a nonzero curve/status result, rather than reproducing every R
  special-value accident. Real arrays and scalar parameters replace R coercion
  and recycling. Frequencies must have the same length as observations.
- Signed finite observations follow the upstream arithmetic wherever defined.
  Economically meaningful inequality and poverty ordinarily require nonnegative
  incomes. Logarithms and fractional powers of negative quantities can yield NaN.
- Preserve the nonstandard upstream details: `Theil(parameter=0)` removes zeros
  before computing its mean; `Foster` uses exponent `parameter-1`; strict poverty
  comparisons omit observations equal to the threshold; `Sen` can return NaN
  when every poor income is zero; `Lasym` preserves its historical missing-value
  weight alignment and its interpolation's addition of one frequency unit.
- Stable scaling, log-domain algebra, and `log1p`/`expm1` improve overflow,
  underflow, and small-parameter behavior. A positive series stabilizes the
  exponential Lorenz curve near zero. Ordinary Lorenz and Mehran arithmetic
  is preserved before resorting to scaling, since exact slope ties matter.
  Kolm remains finite in many cases where upstream exponentiation overflows.
- `Lc.exp(1)` deliberately returns the mathematical limit **1**, rather than
  ineq 0.2-13's `0 * log(Inf)` NaN. Theoretical probabilities outside [0,1] yield
  NaN instead of extrapolation. Finite inadmissible theoretical parameters use
  the upstream fallback where its guards specify one, without emitting warnings.
- Undefined empirical curves (empty retained input, zero total income or
  frequency, nonfinite observations, mismatched weights) return empty coordinate
  arrays and `ineq_invalid`. R sometimes returns partial NaN coordinates or
  errors for these inputs. Mehran requires at least two retained classes.
- `major` returns status for incomparable lengths/totals. It uses exact total
  comparisons, not a tolerance, and does not claim bit-for-bit equality with R's
  extended-precision accumulators. Compensating rounding can affect very close
  noninteger cases in either implementation.

## Validation

`fpm test` runs frozen, actual-R references and independent mathematical checks,
without requiring R at build time:

- 966 assertions against R ineq 0.2-13: all measure/dispatcher mappings, positive,
  zero, signed, negative, tied, singleton and missing observations, corrected
  Gini, parameter branches and limits, poverty thresholds, weighted curves,
  Mehran examples, asymmetry intervals, beta parameters and theoretical tails.
- 477 independent assertions: pairwise and Lorenz-integral Gini identities,
  concentration/variance relations, scale invariance over 600 orders of
  magnitude, poverty-gap identities, small-parameter limits, overflow-safe
  income-frequency products, monotonic theoretical curves, degenerate beta
  cases, invalid-input statuses, and exact-mean/majorization regressions.

Reference comparisons use `2e-12 * max(1, abs(expected))`, except explicit
near-singular scalar parameters use `5e-9` because the original R power and
exponential formulas lose precision. Independent limit tests separately check
the stabilized answers. `validation/R-reference.tsv` contains the numerical
oracle, and `validation/R-session-info.txt` records the runtime. Recreate it,
from this directory with R and the exact ineq 0.2-13 package installed, using:

```sh
Rscript validation/generate_reference.R
fpm test
```

The generator intentionally rewrites the reference test and oracle. Default,
checked-debug, and release build/test/all-example runs are required for release.
[VALIDATION.md](VALIDATION.md) gives commands and the verification record.

## Dependencies and licensing

`rfortran-core` provides the common real kind, stable sorting, normal CDF and
quantile, regularized beta, logarithmic transforms and log-sum-exp helpers. A
repository survey found these compatible APIs, so no translated stats package,
BLAS, LAPACK, `r.f90`, or dependency source is embedded here.

The original package is licensed **GPL-2.0-only OR GPL-3.0-only**. This translation
uses the same choice; both complete license texts are included. See [NOTICE](NOTICE)
and [PROVENANCE.md](PROVENANCE.md) for upstream attribution, archive hash,
repository dependency pin, and retained source-file provenance.
