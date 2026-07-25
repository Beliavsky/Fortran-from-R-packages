# API Map

## R to Fortran mapping

| Original R routine | Modern Fortran routine | Status |
|---|---|---|
| `NelsonSiegel(rate, maturity, doplot)` | `fit_nelson_siegel(rate, maturity, fit, ...)` | Computational path translated and tested |
| Internal Nelson-Siegel `fx` | `nelson_siegel_curve(maturity, parameters)` | Translated and tested |
| Internal Nelson-Siegel `func` | Internal SSE objective | Translated and tested |
| Nelson-Siegel grid search | Internal log-spaced grid plus SVD least squares | Translated with documented basis correction |
| `Svensson(rate, maturity, doplot)` | `fit_svensson(rate, maturity, fit, ...)` | Computational path translated and tested |
| Internal Svensson `fx` | `svensson_curve(maturity, parameters)` | Translated and tested |
| Effective Svensson `func` | L1 objective, default | Translated and tested |
| Earlier SSE expression in Svensson `func` | `objective="sse"` | Exposed and tested as an optional objective |
| Svensson two-dimensional grid search | Internal positive-decay grid plus SVD least squares | Translated and tested |
| `nlminb` | Bounded Nelder-Mead with log decay parameters | Tested numerical analogue |
| Plotting calls | None | Excluded plotting infrastructure |

## Result mapping

The R routines return the optimizer list. `term_structure_fit` exposes the
numerically useful equivalent fields:

- `parameters`
- `fitted`
- `residuals`
- `objective_value`
- `sse`
- `mae`
- `rmse`
- `converged`
- `iterations`
- `evaluations`
- `status`

## Complete computational scope

The original package has no additional numerical bond-pricing, duration,
convexity, cash-flow, or fixed-income routines. Its complete computational
surface is the two term-structure fits above.
