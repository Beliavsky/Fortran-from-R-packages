# RLRsim (Fortran translation)

This directory translates the portable computational core of the R package **RLRsim 3.1-9** to modern free-form Fortran. RLRsim implements simulation-based finite-sample likelihood-ratio and restricted likelihood-ratio tests for a single variance component using the distribution derived by Crainiceanu and Ruppert (2004).

The Fortran API works directly with numeric design matrices rather than R `lm`, `lme`, or `lmerMod` objects. It does not include plotting, interactive code, S3 classes, model-frame/formula processing, or parallel R backends.

## Build

```text
fpm build
fpm test
fpm run --example basic_simulation
```

The implementation has no external BLAS/LAPACK requirement. It uses a package-local real kind (`dp = real64`), modified Gram-Schmidt residualization, a symmetric Jacobi eigensolver for squared singular values, and Fortran intrinsic random numbers with gamma/chi-square transforms.

## Public API

- `lrt_sim`: simulation underlying R `LRTSim`.
- `rlrt_sim`: simulation underlying R `RLRTSim`.
- `exact_lrt_from_design`: model-independent equivalent of the computational portion of `exactLRT` when the observed LRT statistic is already available.
- `exact_rlrt_from_design`: model-independent equivalent of the computational portion of `exactRLRT` when the observed RLRT statistic is already available.
- `rlrsim_result`: simulated statistics plus maximizing lambda values.
- `exact_test_result`: observed statistic, simulation p-value, simulated null sample, and maximizing lambda values.

Seeded runs are deterministic within this implementation, but the stream is not intended to reproduce R's RNG bit-for-bit. Tests cover the general simulation, both specialized `RLRTSim` approximation branches, the general eigenvalue-truncation branch, positive `lambda0`, and rank-deficient fixed-effect designs.

## Translation coverage

Package status: **substantial**. Coverage is **4 of 5 (80%)** exported computational R functions. This fraction measures mapped computational functions, not complete R interface or behavioral compatibility.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `LRTSim` | `rlrsim_api` | `lrt_sim` | substantial |
| `RLRTSim` | `rlrsim_api` | `rlrt_sim` | substantial |
| `exactLRT` | `rlrsim_api` | `exact_lrt_from_design` | substantial |
| `exactRLRT` | `rlrsim_api` | `exact_rlrt_from_design` | substantial |
| `extract.lmeDesign` | - | - | untranslated |

`extract.lmeDesign` remains untranslated because it is tightly coupled to R's `nlme` fitted-object, formula, grouping-factor, and model-matrix machinery. The unexported `extract.lmerModDesign` is outside the stated coverage denominator. The translated exact-test helpers accept the extracted matrices and observed likelihood-ratio statistic explicitly.

For `RLRTSim`, all three upstream `use.approx` paths are included: the balanced-ANOVA closed-form shortcut, the one-dominating-eigenvalue shortcut, and the general leading-eigenvalue truncation path. The same six-decimal balanced-spectrum test and dominating-eigenvalue ratio test are used as in the R implementation.

## Numerical compatibility notes

The upstream C++ kernel initializes each simulated likelihood-ratio maximum at zero, searches an increasing lambda grid, and stops after the objective first decreases. The Fortran kernel preserves that behavior. For the common `lambda0 = 0` case, the standard grid is `0` followed by log-equidistant positive values. For positive `lambda0`, the translated asymmetric grid now follows the R construction and contains exactly the requested `gridlength` points, including zero and `lambda0`.

The upstream `RLRTSim.R` contains an expression `log(-10 * lambda0)` in the branch that expands a too-high lower grid limit for positive `lambda0`; that expression is not real-valued. This translation uses `log(lambda0 / 10)`, matching the apparent intent to extend the lower endpoint to one decade below `lambda0`.

## References

Crainiceanu, C. M. and Ruppert, D. (2004). Likelihood ratio tests in linear mixed models with one variance component. *Journal of the Royal Statistical Society: Series B*, 66, 165-185.

See `NOTICE.md`, `upstream/DESCRIPTION`, and `upstream/CITATION` for upstream authorship and provenance.
