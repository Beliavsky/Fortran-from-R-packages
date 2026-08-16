# CCd R to Fortran API map

| R API | Fortran API | Notes |
|---|---|---|
| `dcc` | `dcc` | Elemental PMF/log-PMF; `logged` maps directly. |
| `pcc` | `pcc` | Scalar integer CDF formula matching upstream implementation. |
| `qcc` | `qcc` | Elemental left quantile; uses bracketing/binary search rather than `optimize`. |
| `cc.mle` | `cc_mle` | Returns `cc_fit_result`; Nelder-Mead optimization. |
| `cc.mle0` | `cc_mle0` | Golden-section one-dimensional optimization. |
| `cc.reg` | `cc_reg` | Adds an intercept internally, matching `model.matrix(y ~ .)`. |
| `loc0.test` | `loc0_test` | Returns statistic and chi-square(1) upper-tail p-value. |

There are no plotting exports in the supplied package.
