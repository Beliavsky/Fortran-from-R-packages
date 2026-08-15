# Porting notes

Upstream: Delaporte R package 8.4.3 (2026-01-08).

## Mapping

| Upstream operation | Standalone Fortran |
|---|---|
| `ddelap` / `ddelap_f_s` | `ddelap`, `ddelap_vec` |
| `pdelap` / `pdelap_f_s` | `pdelap`, `pdelap_vec` |
| `qdelap` / `qdelap_f_s` | `qdelap`, `qdelap_vec` |
| `qdelap(..., exact=FALSE)` | `qdelap_approx` |
| `rdelap(..., exact=TRUE)` | `rdelap(..., exact=.true.)` |
| `rdelap(..., exact=FALSE)` | `rdelap(..., exact=.false.)` |
| `MoMdelap` / `momdelap_f` | `momdelap` |
| R/C registration wrappers | omitted |
| OpenMP get/set wrapper API | omitted |

The PMF, CDF, exact quantile, and method-of-moments formulas are direct
adaptations of the upstream compiled Fortran. The exact random path still uses
uniform-to-quantile inversion, but uses Fortran `random_number` instead of R's
RNG. The mixture random path uses standalone Gamma and Poisson generators.

`qdelap_approx` preserves the upstream Monte-Carlo construction and R type-8
sample-quantile interpolation. Because the underlying RNG is not R's RNG,
identical seeds do not imply identical simulated samples across R and Fortran.

The R wrapper's interactive guard for `pdelap(q >= 32768)` is intentionally not
ported. The standalone routine also guards values beyond the exact `int64`
range before integer conversion.
