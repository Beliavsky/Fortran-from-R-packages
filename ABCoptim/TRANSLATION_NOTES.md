# Translation notes

## Upstream

- Package: ABCoptim
- Version: 0.15.0
- Date: 2017-11-05
- License: MIT + upstream LICENSE file
- Main computational sources:
  - `R/abc_optim.R`
  - `src/abc.cpp`

## Source mapping

The following package-owned numerical operations are translated directly:

| ABCoptim source | Fortran implementation |
| --- | --- |
| `CalculateFitness` | `calculate_fitness` |
| `initial` / C++ initialization | `initialize_foods` |
| `init` | `initialize_source` |
| `SendEmployedBees` | `employed_bees` / `propose_and_select` |
| `CalculateProbabilities` | `calculate_probabilities` |
| `SendOnlookerBees` | separate R and C++ onlooker loops |
| `SendScoutBees` | `scout_bee` |
| `MemorizeBestSource` | `memorize_best` |
| `abc_optim` | public `abc_optim` |
| `abc_cpp_` + R wrapper | public `abc_cpp` |

The objective callback is a Fortran procedure. R `...` arguments are naturally
represented through a module procedure or an internal procedure using host
association.

## Preserved behavioral differences

The pure-R and Rcpp implementations are not numerically identical. The port
therefore does not force them into one behavior. Important source differences
preserved include binary mutation in the R implementation, its onlooker index
progression, persistence comparison (`>` versus `>=`), history placement, and
the different small constants used in probability normalization.

The R code evaluates its initial `GlobalMin` as raw `fn(par)` before defining
and using `fn(par/parscale)/fnscale` for the population. This can mix scales
when `parscale` or `fnscale` is nontrivial. The public R-style entry point keeps
this by default through `legacy_r_scaling`; users may disable it.

The R binary implementation still begins with the same equally spaced
continuous initial food matrix as the source package. Binary 0/1 values enter
through mutations and scout initialization. This is deliberate compatibility,
not a typo in the Fortran port.

## Intentional corrections

The R wrapper around `abc_cpp_` uses `rep(lb, length(par))` and the same for
`ub`, which can expand an already vector-valued bound to length `D*D`. The
Fortran API accepts either a scalar bound or a length-D bound and rejects other
sizes.

At exact `maxCycle` exhaustion, the upstream R/C++ loops have an off-by-one
count/history edge case. The Fortran port retains the intended maximum number
of actual ABC generations (`max_cycle-1` after initialization) but always
returns a valid trimmed history.

Non-finite objective values are rejected consistently rather than reproducing
the R/C++ difference in where NA checking happens.

## Omitted non-computational code

- R S3 `print.abc_answer`.
- R S3 `plot.abc_answer` and graphics calls.
- `.Call` and generated Rcpp registration glue.
- R list/class construction and help-page infrastructure.

## Numerical/runtime choices

The port is self-contained and uses no BLAS/LAPACK. Random numbers come from a
Park-Miller generator with explicit state so runs are reproducible in Fortran.
This intentionally does not emulate R's RNG bit-for-bit.
