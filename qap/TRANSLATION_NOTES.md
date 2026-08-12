# Translation notes

## Original package

- Package: `qap`
- Version: 0.1-2
- License: GPL-3
- Main computational source: `src/qapsim.f`
- R-level computational helpers: `R/qap.R`, `R/qapSA.R`, `R/read_qaplib.R`

The original source is retained verbatim under `original/qap-master/`.

## Modernization

The legacy `qaph4` routine was rewritten as free-form Fortran 2018 with:

- `implicit none` everywhere;
- explicit module interfaces;
- `real64`/`int64` kinds from `iso_fortran_env`;
- allocatable result/problem types;
- no R headers or ABI dependency;
- no implicit external procedures;
- a standalone deterministic RNG;
- explicit matrix/permutation validation.

The inner simulated-annealing equations follow `qapsim.f`, including the
symmetric swap-delta formula and its cooling/stopping logic.

## Deliberate safety hardening

The old routine initializes the best objective to a rough mean estimate and
leaves `ope` unwritten until an observed state beats that threshold. On a
pathological problem it could therefore return an uninitialized permutation.
The translation preserves the threshold behavior, but initializes the output
permutation to the starting permutation and falls back to its true objective if
no state ever reaches the threshold.

The R wrapper accepts `ft=0` even though its message describes `(0,1)`. The
Fortran API also permits zero. When temperature reaches zero, worsening moves
are rejected explicitly instead of relying on floating-point division by zero.

## RNG

The R implementation calls `sample(n)` for each restart and R's `unif_rand()`
inside the Fortran routine. The standalone translation uses a Park-Miller RNG
with Schrage arithmetic. Thus identical Fortran seeds are portable and
reproducible, but not bit-identical to R's `set.seed()` stream.

## QAPLIB

The complete `inst/qaplib` data directory is copied to `data/qaplib`. The
reader follows the R implementation and performs no interpretation of the
orientation of solution permutations. This matters for a few historical
QAPLIB `.sln` files whose permutation is expressed in the inverse convention.

## Omitted code

Only R-specific glue is omitted:

- `src/init.c` registration;
- `src/RNG_wrapper.c` wrappers around R's RNG;
- R list/attribute handling.

There is no plotting code in the package.
