# Computational coverage

Upstream package: `PortfolioOptim` 1.1.1.

| Upstream routine | Fortran routine | Status |
|---|---|---|
| `BDportfolio_optim` | `bdportfolio_optim` | Implemented |
| `PortfolioOptimProjection` | `portfolio_optim_projection` | Implemented |
| `.RISK_post` | `risk_post` | Implemented |
| `.make_diag` | `diagonal_matrix` | Implemented |
| `.F_func` | `f_func` | Implemented |
| `.ZI_projection` | `zi_projection` | Implemented |
| `Rsymphony_solve_LP` dependency | `solve_lp` | Replaced by native two-phase simplex |

## Additional Fortran interfaces

- `quadratic_lp_projection`: stable least-distance projection onto an LP optimal
  set using a simplex first stage and ADMM quadratic second stage.
- `risk_measure`: direct evaluation of the four supported risk measures.
- `risk_code`: case-insensitive conversion from package risk names.
- typed `risk_result`, `lp_result`, `projection_result`, and `portfolio_result`
  structures.

## Included risk measures

- Conditional Value-at-Risk (`CVAR`)
- Deviation Conditional Value-at-Risk (`DCVAR`)
- Lower Semi Absolute Deviation (`LSAD`)
- Mean Absolute Deviation (`MAD`)

## Deliberately excluded

The upstream package has no plotting or compiled-code layer. The following R
infrastructure is not translated because it is not a numerical algorithm:

- R list, matrix-label, and printing behavior;
- `match.arg`, `all.equal`, and R input-dispatch conventions;
- CRAN package metadata and roxygen-generated help rendering; and
- optional example-only dependencies such as `mvtnorm` and `Rglpk`.

The unmodified R source and documentation remain under `original/`.
