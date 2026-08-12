# Rmalschains-fortran 0.1.0

Modern Fortran 2018 translation of the computational path exposed by the R package
**Rmalschains 0.2-11**.

Rmalschains implements the MA-LS-Chains family of memetic algorithms for bounded
continuous optimization.  The R interface always minimizes and constructs a
steady-state genetic algorithm (SSGA) hybridized with a persistent local-search
chain.

## Implemented R-facing computation

The Fortran port reproduces the numerical path constructed by
`RmalschainsWrapper.cpp`:

* bounded real-valued population initialization;
* steady-state GA;
* NAM-3 parent selection;
* BLX-alpha crossover;
* BGA mutation with the librealea default mutation probability 0.125;
* replace-worst survivor selection;
* cumulative MA-LS-Chains EA/local-search effort balancing using the original
  `calculateFrec()` equation;
* selection of the best individual not previously marked as a failed local
  search;
* persistent per-individual local-search state (the defining "LS chains"
  behavior);
* the documented local searches `cmaes`, `sw`, `simplex`, and `ssw`;
* the additional native-selector local searches `mts1`/`mtsls1` and `mts2`;
* target-value stopping, relative-improvement thresholding, supplied initial
  populations, and the package result counters/timing diagnostics.

The implementation is self-contained except for BLAS/LAPACK, used by the
CMA-ES covariance eigensolver.

## Public API

```fortran
use rmalschains, only : dp, mals_control, mals_result, &
  malschains_control, malschains_optimize

type(mals_control) :: control
type(mals_result) :: result
real(dp) :: lower(10), upper(10)

lower = -5.0_dp
upper =  5.0_dp

control = malschains_control(popsize=20, ls='cmaes', istep=150, &
  effort=0.5_dp, alpha=0.5_dp, seed=1234_8)

result = malschains_optimize(objective, lower, upper, 8000, control)
```

For `initialpop`, Fortran uses shape `(dimension, number_of_individuals)`.

## Local-search names

* `cmaes`, `cmaesmyrandom`, `cmaesalways`
* `sw`
* `ssw`
* `simplex`
* `mts1`, `mtsls1`
* `mts2`

The R wrapper rewrites the documented `cmaes` choice to `cmaesmyrandom` so that
its seed is respected.  The Fortran implementation always uses its explicit
`control%seed` RNG state, so the three bounded CMA-ES aliases are numerically
identical here.

## Compatibility details

### Population size

`malschains_control(popsize=...)` follows the R helper: population sizes are
rounded to the nearest multiple of 10 and a rounded value of zero becomes 10.

### Initial-population evaluation accounting

The R wrapper evaluates supplied initial-population rows through `Problem::eval`
after initialization, bypassing the `Running` evaluation counter.  Therefore
those evaluations do not appear in `numEvalEA` or `numEvalLS`.

This port preserves those counters and additionally reports:

```fortran
result%actual_nfe
```

which includes all objective callback invocations.

### `lsOnly=TRUE` zero-start quirk

Without an initial population, the original wrapper creates a random population
but starts local search from a zero vector whose stored fitness is also set to
zero **without evaluating the objective**.  This can return an impossible
objective value for a nonnegative function whose optimum is not at zero.

For source compatibility the default is:

```fortran
control%legacy_ls_only_zero_start = .true.
```

For normal numerical use, set:

```fortran
control%legacy_ls_only_zero_start = .false.
```

Then the random setup population is evaluated and its best member is used.
Those setup evaluations are visible in `actual_nfe`.

### SSW initialization

The original C++ `SWN2Dim::getInitOptions()` tests an uninitialized
`option->delta[i]` while trying to clamp the newly computed `delta_init[i]`.
Undefined memory state cannot be translated reproducibly.  This port performs
the evident intended operation: `delta_init` is clamped to `[1e-15, 0.4]`.

### Random-number streams

The original package uses librealea's `SRandom`/CMA-ES RNG machinery.  This
translation uses a deterministic portable Fortran RNG.  Equal integer seeds are
reproducible within this port but do not imply bit-identical trajectories to R.

## Translation boundary

Rmalschains vendors roughly 40,000 lines from the older librealea framework,
including many experimental evolutionary algorithms (`JADE`, `jDE`, `SaDE`,
PSO, CHC, etc.) that are not instantiated by the exported R function
`malschains()`.

Version 0.1.0 translates the **complete computational path reachable from the R
package API**, plus MTS1/MTS2 paths already accepted by the package's native
local-search selector.  Dormant librealea algorithms that cannot be selected by
`malschains()` are retained under `original/` for provenance but are not exposed
as Fortran APIs in this release.

Similarly, the CMA-ES chain follows the same Hansen covariance/path/step-size
adaptation equations and persistent-state semantics, but is expressed as a
modern Fortran implementation rather than a line-for-line conversion of the
2,891-line historical `origcmaes.cc` implementation.

## Building

With FPM and BLAS/LAPACK available:

```text
fpm build
fpm test
fpm run --example sphere
```

The supplied `fpm.toml` links `lapack` and `blas`.

## Validation

The release test suite covers:

1. full MA-LS-Chains with CMA-ES local chains on a sphere objective;
2. all six translated local-search families;
3. a 10-dimensional Rastrigin hybrid optimization;
4. the historical `lsOnly` zero-start behavior and corrected mode;
5. initial-population and evaluation-accounting semantics.

On the release build, the 5-D CMA-ES-chain sphere regression reaches about
`1.64e-12`, with 2500 counted EA evaluations and 2500 counted LS evaluations.
The 10-D simplex/Rastrigin regression reaches about `1.66e-9`.

See `VALIDATION.md` for compiler details.

## Licensing

Rmalschains/librealea is GPL-3.  The original package also contains code with
additional provenance/license notices (RcppDE-derived code, Nikolaus Hansen's
CMA-ES, Robert Davies' newmat, and Richard J. Wagner's ConfigFile class).

The original package tree and its `inst/original_licenses_of_code_parts` file
are retained under `original/`.  See `LICENSES.md` and `COPYING`.
