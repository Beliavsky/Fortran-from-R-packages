# Translation notes

## Source mapping

| R source | Fortran implementation |
|---|---|
| `R/cmaes.R:cma_es` | `src/cmaes.f90:cma_es` |
| `R/cmaes.R:extract_population` | `src/cmaes.f90:extract_population` |
| `R/functions.R:f_sphere` | `src/cmaes_functions.f90:f_sphere` |
| `R/functions.R:f_rand` | `src/cmaes_functions.f90:f_rand` |
| `R/functions.R:f_rosenbrock` | `src/cmaes_functions.f90:f_rosenbrock` |
| `R/functions.R:f_rastrigin` | `src/cmaes_functions.f90:f_rastrigin` |
| `shift_function` | `shifted_value` |
| `rotate_function` | `rotated_value` |
| `bias_function` | `biased_value` |
| `stats::rnorm` | local Box-Muller generator in `cmaes_rng` |
| `eigen(..., symmetric=TRUE)` | LAPACK `DSYEV` in `cmaes_linalg` |

## Numerical representation

The algorithm remains double precision. `C`, `B`, and `BD` are dense matrices,
matching the original implementation. Eigenvalues returned by LAPACK are
reversed into descending order so the update and positive-definiteness check
have the same ordering used by R's symmetric `eigen()`.

The RNG is intentionally self-contained instead of calling R. This changes the
random trajectory but not the CMA-ES equations.

## Box constraints

The original evaluates the objective at the clamped point `vx`, multiplies its
fitness by a squared-distance penalty, and then performs selection/recombination
with the original sampled point `arx`. The translation retains this distinction.

Because `pen = 1 + squared_error` is evaluated in finite precision, a squared
violation smaller than half an ulp near one can round away. Such a trial then
passes the source's `pen <= 1` validity test. Regression checks therefore allow
roundoff-scale boundary excursions instead of changing the original formula.

## Diagnostics

The R code sorts `arfitness`, keeps `aripop` as indices into the original
population, then records `arfitness[aripop]`. This is a second indexing of the
already-sorted vector. The Fortran diagnostic history intentionally preserves
that observable behavior. Selection itself uses `aripop` correctly on `arx`
and `arz`.

## Omitted R-only behavior

- S3 result class and deprecated `cmaES()` R wrapper
- parameter names/rownames
- `...` argument forwarding (Fortran callers can use internal procedures or
  module state to capture additional objective data)
- R warnings/messages and `NULL` values

When R would leave `best.par=NULL` because no unpenalized point has ever been
accepted, the Fortran result uses a zero-length `par` array.
