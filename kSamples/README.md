# kSamples - modern Fortran translation

This directory translates the computational core of the R package **kSamples 1.2-12** to modern free-form Fortran with FPM. kSamples provides k-sample rank tests and exact/randomization distributions, including Anderson-Darling, QN/Kruskal-Wallis score tests, Steel multiple-comparison tests, Jonckheere-Terpstra tests, and 2 x t contingency-table combinations.

The public Fortran facade is module `ksamples`. The API uses explicit numeric arrays and derived result types rather than R formulas, lists, S3 objects, or printing methods.

## Build

Place this directory at the top level of `Fortran-from-R-packages` beside `rfortran-core` and `SuppDists`, then run:

```text
fpm build
fpm test
fpm run --example rank-tests
```

No separately installed BLAS or LAPACK library is required. The current translation uses `rfortran-core` for shared R-like numerical helpers and the existing `SuppDists` translation for normal-order scores.

## Public API

The main procedures are:

- `ad_pval`, `ad_test`, `ad_test_combined`
- `qn_test`, `qn_test_combined`
- `steel_test`, `steel_confint`
- `jt_test`, `djt`, `pjt`, `qjt`
- `contingency2xt`, `contingency2xt_comb`
- `conv`

Result types (`ad_result`, `qn_result`, `jt_result`, `steel_result`, `contingency_result`, `combined_result`, and `steel_confint_result`) are also re-exported from `ksamples`.

## Translation coverage

Package status: **substantial**.

**14 of 14 (100%)** exported computational R functions have a meaningful Fortran mapping. This percentage measures mapped computational functions; it does **not** mean complete R interface or behavioral compatibility. R-only formula/list dispatch, NA-removal warnings, S3 construction/printing, graphics, and bit-identical R RNG behavior are intentionally omitted. `ad.pval` and `SteelConfInt` retain material differences documented below.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `ad.pval` | `ksamples_ad` | `ad_pval` | partial |
| `ad.test` | `ksamples_ad` | `ad_test` | substantial |
| `ad.test.combined` | `ksamples_ad` | `ad_test_combined` | substantial |
| `contingency2xt` | `ksamples_contingency` | `contingency2xt` | substantial |
| `contingency2xt.comb` | `ksamples_contingency` | `contingency2xt_comb` | substantial |
| `conv` | `ksamples` | `conv` | complete |
| `qn.test` | `ksamples_qn` | `qn_test` | substantial |
| `qn.test.combined` | `ksamples_qn` | `qn_test_combined` | substantial |
| `Steel.test` | `ksamples_steel` | `steel_test` | substantial |
| `SteelConfInt` | `ksamples_steel` | `steel_confint` | partial |
| `djt` | `ksamples_jt` | `djt` | complete |
| `pjt` | `ksamples_jt` | `pjt` | complete |
| `qjt` | `ksamples_jt` | `qjt` | substantial |
| `jt.test` | `ksamples_jt` | `jt_test` | substantial |

The `ad_pval` port preserves the upstream tabulated calibration but replaces R `smooth.spline` interpolation with deterministic linear interpolation. `steel_confint` implements the asymptotic simultaneous confidence-bound inversion; exact/simulated confidence-level refinement remains untranslated. See `docs/API_COVERAGE.md` and the per-function entries in `fpm.toml`.

## Numerical and interface notes

Exact tests enumerate reallocations conditional on the pooled data/tie pattern. When the exact allocation count exceeds `nsim`, the test routines fall back to simulation, mirroring the upstream package's computation guard. Simulation uses Fortran `random_number`; deterministic tests do not depend on RNG output.

The Fortran API expects already-clean numeric arrays. It does not silently remove NaNs. Callers should perform missing-value policy before invoking the tests.

`dp` is imported from `r_kinds` in `rfortran-core` and used consistently for real arithmetic. The source does not use `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent literals.

## Provenance and license

The upstream computational source snapshot used for the translation is under `upstream/`; it is retained for provenance and is not part of the FPM build. kSamples is GPL (>=2), and this translation is GPL-2.0-or-later. See `NOTICE.md` and `LICENSE`.
