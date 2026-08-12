# matchingMarkets-fortran

Modern Fortran/FPM translation of the computational matching mechanisms in
R package `matchingMarkets` 1.0-5.

The project is native Fortran. It does not use Rcpp, Armadillo, Java, rJava,
or the package's Java stable-matching launcher.

## Implemented in v0.1.0

### Matching mechanisms

- `hri` -- exact enumeration of stable hospital/resident or stable-marriage
  assignments, including capacities and incomplete numeric preference lists.
  The Fortran implementation uses exhaustive bounded enumeration rather than
  the original Java/constraint encoding, so it is intended for small and
  moderate enumeration problems.
- `hri2` -- exact small/medium-market couples fallback. Couples are assigned
  jointly from their listed college-pair preferences and candidate solutions
  are checked for individual and joint blocking opportunities. The original
  Roth-Peranson/Bacchus C++ matcher remains in `original/`; this reference
  implementation is not intended to reproduce its large-market performance.
- `hri3` -- efficiency-adjusted deferred acceptance / immediate acceptance,
  translated from the package's `eadam.cpp` round-level algorithm.
- `iaa` -- Boston immediate acceptance and Gale-Shapley deferred acceptance,
  with heterogeneous college capacities and incomplete preference lists.
- `sri` -- exact enumeration of stable-roommates solutions for small/moderate
  problems.
- `rsd` -- random serial dictatorship with capacities.
- `ttc` -- top trading cycles with existing tenants/vacant houses.
- `ttc2` -- school-market top trading cycles with capacities.
- `ttcc` -- top trading cycles and chains for kidney exchange.
- `stabchk` -- blocking-pair detection for two-sided assignments.

The already validated `matchingR-fortran` implementation is vendored as an
FPM path dependency and also exposes Gale-Shapley, college admissions, Irving
stable roommates, and standard top-trading-cycles APIs.

### Statistical and optimization utilities

- `plp` -- the transferable-utility stable-roommates partitioning binary LP,
  using the vendored `lpSolve-fortran` dependency.
- `khb` -- matrix-level Karlson-Holm-Breen probit coefficient comparison.
- `probit_fit`, `ols_fit`, and inverse-Mills-ratio helpers.
- `stabsim` and `stabsim2` native matching-data simulation helpers.
- `pair_combinations`, `coalition_partitions`, and `consensus_mc`, translating
  computational helpers used by the structural-model R layer.

### R-name compatibility

The module provides generic names `hri`, `hri2`, `hri3`, `sri`, `ttc`, `ttc2`,
`ttcc`, and `stabchk`, in addition to more descriptive implementation names.
Fortran uses 1-based participant/college IDs and `0` for unmatched or missing
preference entries.

## Deliberately deferred

The large `stabit()` and `stabit2()` structural Gibbs-estimation kernels are
not represented by a generic Heckman approximation in this release. They are
the main remaining computational family for a future version. Their complete
R/Rcpp source is retained under `original/matchingMarkets-master/`.

R formula/data-frame handling, S3 methods, plotting, Java/rJava glue, progress
bars, printing and parallel orchestration are also omitted.

## Build

```sh
fpm build
fpm test
```

BLAS/LAPACK are linked for the matrix-level statistical routines. The project
vendors:

- `vendor/matchingr-fortran` -- GPL-2.0-or-later
- `vendor/lpsolve-fortran` -- LGPL-2.0-only

## Example

```fortran
use matchingmarkets

type(assignment_result_t) :: ans
integer :: sp(2,4), cp(4,2), slots(2)

sp(:,1) = [1,2]
sp(:,2) = [1,2]
sp(:,3) = [2,1]
sp(:,4) = [2,1]
cp(:,1) = [2,1,3,4]
cp(:,2) = [4,3,1,2]
slots = [1,1]

ans = iaa(sp,cp,slots,'deferred')
print *, ans%assignment
```

## Validation

The source is tested with GNU Fortran 14.2.0 using both `-O2` and:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

See `VALIDATION.md` for the exact test coverage.

## License

`matchingMarkets` declares GPL (>= 2), therefore the new Fortran translation
is GPL-2.0-or-later. Original source/copyright material is retained under
`original/`. Vendored dependencies retain their own licenses; see
`LICENSES.md`.
