# diptest: corrected Hartigan dip test in Fortran

Modern, free-form Fortran translation of **diptest 0.77-2**, distributed under
GPL-2.0-or-later. It includes the corrected AS 217 dip statistic, sorted data,
modal interval, complete GCM/LCM diagnostics, the original 21-by-26 simulation
quantile table, and tabulated or Monte Carlo p-values.

The numerical source is the current corrected **C** implementation, not the old
Fortran code containing the misplaced-parenthesis defect. The 1994 termination
fix and Yong Lu's 2003 correction are retained.

## Build and run

Place this directory at the root of
[Beliavsky/Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages),
next to `rfortran-core`, then run:

```sh
cd diptest
fpm build
fpm test
fpm run --example --all
```

FPM and gfortran are required. R, a C compiler, BLAS and LAPACK are **not** required
for normal builds, tests or examples. The sibling FPM dependency provides the
single `dp` kind and pure sorting. Its source is not copied here. No system
library links or platform-specific source are used. Linux gfortran 14.2/FPM
0.13.0 were tested; a native Windows build was not executed.

```fortran
use diptest_api, only : dp, dip_test, dip_test_result, dip_ok
type(dip_test_result) :: result
call dip_test([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], result)
if (result%status /= dip_ok) error stop 'invalid sample'
print *, result%statistic, result%p_value
```

## API

All public procedures except the intrinsic-RNG convenience routine are pure.
All real values use `dp` imported from `r_kinds` and re-exported by `diptest_api`.
See [API.md](API.md) for argument and result details.

- `dip(x [, min_is_zero])`: scalar statistic, sorting a copy of finite `x`.
- `dip_compute(x, result [, min_is_zero, sorted])`: statistic and complete
  `dip_result` diagnostics. `sorted=.true.` checks the presorted contract.
- `dip_test(x, result [, uniform_samples])`: tabulated p-value by default, or
  reproducible Monte Carlo p-value using the supplied `n`-by-`B` null samples.
- `dip_test_simulate(x, result [, replicates])`: streaming intrinsic Fortran RNG
  simulation with `B=2000` by default; callers can seed with `random_seed`.
- `dip_p_value(statistic, n, p_value [, status, asymptotic])`: calibrate an
  already-computed statistic against the original table.

`dip_result` exposes sorted `x`, `statistic`, `n`, `status`, the 1-based `low` and
`high` indices, original-coordinate `lower` and `upper`, modal `gcm` and `lcm`,
full-sample `full_gcm` and `full_lcm`, and all `mn`/`mj` predecessor/successor
indices. GCM indices are descending; LCM indices are ascending.

## Statistical and numerical compatibility

- The default classical minimum is `1/(2*n)`, including constant samples and
  singletons. `min_is_zero=.true.` matches R's alternative convention.
- `dip` rejects NaNs and infinities. `dip_test` removes NaNs like R's
  `complete.cases` and rejects infinities. Empty retained samples return an
  explicit error, rather than propagating the R wrapper's ignored C fault.
- Tests with `n<=3` return p=1 and ignore simulation options, as R does.
- Table interpolation is linear between rows after multiplying by `sqrt(n)`,
  then linear in the statistic with constant boundary extrapolation. Repeated
  quantile knots are averaged as by `stats::approx(ties=mean)`. Consequently a
  classical minimum dip at `n=4` has p=0.782, not 1; this is upstream behavior.
- At `n=72000` the final table row is used. Larger `n` uses that same scaled row
  and sets `asymptotic=.true.`; no message is printed.
- Simulation returns `count(D_observed <= D_null)/B`, without a plus-one
  correction. Caller-supplied uniform matrices reproduce a particular R run's
  null samples. Fortran `random_number` does not reproduce R's RNG or promise
  bitwise reproducibility between compilers. The intrinsic convenience routine
  uses O(n) temporary memory rather than allocating all `n*B` samples.
- The implementation uses **one real kind throughout**. Upstream C uses
  `long double` for selected distances. Last-bit differences and alternate
  equally valid modal/hull selections near exact ties are therefore possible;
  exact cross-platform diagnostic equality is not promised by upstream either.
- Local differences are represented by binary fraction/exponent pairs, and
  ratios are evaluated before multiplication by index counts. This avoids
  overflowing ranges and steep-segment slopes without globally rescaling and
  erasing small clusters beside huge outliers. Subnormal input spacing is
  preserved when gradual underflow is supported. No procedure can recover
  distinctions already rounded away in the caller's data.
- For constant data, defined singleton/endpoint hulls replace upstream's
  uninitialized diagnostic arrays. This deliberately makes full diagnostics
  safe for a degenerate sample.
- R calls, classes, printing, plotting, debug output and language-level argument
  matching are omitted. Numeric failures use status codes and NaN outputs.

Do not compile with `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or flush-to-zero
assumptions: validation and extreme-scale behavior rely on IEEE semantics.

## Translation coverage

Package status: **substantial**. **2 of 2 (100%)** exported computational R
functions are mapped. This fraction measures mapped computational functions,
not complete R interface or bit-for-bit compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `dip` | `diptest_api` | `dip`, `dip_compute` | substantial |
| `dip.test` | `diptest_api` | `dip_test`, `dip_test_simulate`, `dip_p_value` | substantial |

No exported computational functions are untranslated. Registered printing and
plotting methods are presentation-only. Internal `getCM` is implemented within
`dip_compute` and is not separately counted. See
[TRANSLATION_COVERAGE.md](TRANSLATION_COVERAGE.md); the machine-readable inventory
is in `fpm.toml`.

## Verification and provenance

Deterministic tests include 132 actual-R data/reference cases, full hulls and
predecessor arrays for 65 of those cases, 987 table interpolation cases, a
257-replicate R-uniform Monte Carlo oracle, analytic two-point/grid formulas,
independent envelope geometry, error contracts, and extreme-range invariances.
Fixtures include `statfaculty`, `faithful$eruptions`, `precip`, and the historical
Mechler/Yong Lu regression cases. See [VALIDATION.md](VALIDATION.md).

The original C/R sources and table are retained under `upstream/` for attribution
and inspection only; they are not built or runtime dependencies. Optional R
scripts under `tools/` regenerate the table and reference fixtures.

See [NOTICE.md](NOTICE.md), [UPSTREAM_PROVENANCE.md](UPSTREAM_PROVENANCE.md),
[LICENSE](LICENSE), and [test/data/README.md](test/data/README.md) for authors,
source references, license scope and dataset attribution.
