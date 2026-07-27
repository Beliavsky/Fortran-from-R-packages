# Validation record

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Python 3.13 used for release audits
- No external BLAS, LAPACK, LP, statistics, or R runtime dependency

FPM was not installed in the validation environment. The FPM manifest was parsed
as TOML and all targets were also built directly with GNU Fortran.

## Checked build flags

```text
-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace
```

## Optimized build flags

```text
-std=f2018 -O2 -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror
```

## Test programs

| Test | Coverage |
|---|---|
| `test_simplex` | bounded LP optimum and infeasibility detection |
| `test_optimizer` | long/short targets, factor/category neutrality, investability, turnover |
| `test_simulation` | internal transfers, market fills, volume limits, transaction costs |
| `test_lifecycle` | multi-day holdings and delisting return/liquidation |
| `test_stats` | rank normalization, grouped normalization, neutralization, drawdown |
| `test_exposures_performance` | factor/category exposure and performance statistics |
| `test_data` | cross-section carry-forward, replacements, statistics, split adjustment |

The demo and both examples are compiled and executed by the validation script.

The release contains 2,607 lines across 19 Fortran source, test, demo, and example files.

## Reference checks

- The simple LP optimum is independently identifiable at `(2, 2)` with objective
  `-10` and is checked to `1e-10`.
- Portfolio optimizer tests have exact integer-share solutions and verify factor
  and category equality constraints after rounding.
- Rank-normal and repeated-neutralization tests use fixed deterministic values.
- Performance annualization and volatility use independently calculated numeric
  references.
- Delisting and P&L tests use direct share-times-return identities.

## Release audits

The release validation checks:

- valid `fpm.toml` syntax;
- ASCII-only translated text;
- maximum Fortran line length of 132 characters;
- `implicit none` in every Fortran program/module;
- `GPL-3.0-only` SPDX headers in every Fortran file;
- no compiled objects or module files in the release tree;
- original/archive/translated SHA-256 manifests;
- a clean rebuild and rerun after extracting the exact final ZIP.
