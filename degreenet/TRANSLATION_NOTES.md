# Translation notes

## Upstream

- Package: `degreenet`
- Upstream version: 1.3-7
- Date in uploaded source: 2026-05-24
- Author/copyright holder: Mark S. Handcock and statnet contributors
- License: GPL-3 plus the Section 7 statnet attribution requirements in
  `LICENSE`

Every Fortran source includes the statnet attribution marker and the example
prints the required attribution.

## R/C to Fortran mapping

| Upstream area | Fortran implementation |
|---|---|
| `src/zeta.c`, `R/zeta.R` | `degreenet_math:zeta_r` |
| `src/cmp.c`, `R/cmp*.R` | CMP routines in `degreenet_distributions` |
| `src/dpoilog.c`, `R/poissonlognormal.R` | `dpln`, `ldpln`, `simpln_one` |
| Yule/Zipf/Waring/DQE/GHDI/NB/PE densities | `degreenet_distributions` |
| stopped Yule/Waring/Zipf models | `degreenet_compound` |
| repeated `ll*` functions | `degreenet_models:loglik_model` |
| repeated `a*mle`, `*mle` wrappers | `degreenet_fit:fit_degree_model` |
| `groupedmodels.R` binning | `grouped_probability`, `grouped_loglik` |
| `llr*` rounding tables | `rounded_bin`, `rounded_probability`, `rounded_loglik` |
| `mands.R` | `degreenet_diagnostics` |
| `bs.R` / parametric bootstrap families | `bootstrap_degree_model` |
| `simdist.R` | `sample_model` plus distribution-specific RNG helpers |
| `reedmolloy.R`, `ryule.R` | `degreenet_graph` |
| plotting / PDF / graphics | intentionally omitted |

## Deliberate numerical/API changes

1. **One common MLE engine.** Upstream repeats `optim()` calls for each family,
   often trying BFGS and Nelder-Mead and keeping the better result.  The
   translation exposes one bounded Nelder-Mead engine with a finite-difference
   observed Hessian.  This preserves the likelihood targets while making the
   Fortran API much smaller and type-safe.
2. **Poisson-lognormal.** Upstream contains both an R Gauss-Hermite
   approximation and a C adaptive integrator.  The Fortran PMF uses normalized
   Gauss-Hermite quadrature (32 points by default), avoiding an R integration
   dependency.  The order is caller-selectable.
3. **CMP.** The normalizing constant is evaluated with a log-sum-exp recurrence
   instead of the C routine's overflow-prone ordinary sum/asymptotic switch.
   It targets the same CMP normalizer.
4. **Conditioning normalization.** A few upstream `ld*` functions divide a log
   probability by a conditioning probability rather than subtracting its log.
   The Fortran density functions use mathematically correct conditional PMFs.
5. **Graph generation.** Upstream calls `igraph::sample_degseq()` and returns a
   `network` object.  The Fortran `reed_molloy()` returns an explicit edge list
   and uses Havel-Hakimi realization; `yule_graph()` repeatedly samples a Yule
   degree sequence until it is graphical.  This preserves the degree-sequence
   task but not igraph's particular random graph draw.
6. **Bootstrap output.** R-specific quantile/list/file-save wrappers are not
   reproduced.  `bootstrap_degree_model()` returns the bootstrap parameter
   matrix so callers can compute any desired summaries.
7. **Grouped and rounded fits.** The translated package exposes the grouped and
   rounded likelihood primitives explicitly.  The many named upstream wrappers
   can be reproduced by optimizing these primitives; they are not duplicated
   one-for-one as separate Fortran subroutine names.

## Model constants

`MODEL_YULE`, `MODEL_DP`, `MODEL_WARING`, `MODEL_DQE`, `MODEL_GHDI`,
`MODEL_NB`, `MODEL_PE`, `MODEL_CMP`, `MODEL_PLN`, `MODEL_GYULE`,
`MODEL_GEODP`, `MODEL_NBYULE`, `MODEL_NBWAR`, `MODEL_GWAR`, `MODEL_GEOM`,
and `MODEL_POIS`.

## Testing

Tests include fixed independent reference values for zeta, Yule, Zipf, Waring,
DQE, CMP (including the Poisson special case), Poisson-lognormal quadrature,
stopped-process probabilities, likelihood improvement under MLE, grouped and
rounded bins, diagnostics, and graph degree realization.
