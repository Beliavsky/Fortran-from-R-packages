# Computational API map

| Upstream R/C routine | Fortran API | Notes |
|---|---|---|
| `int(..., type="Romberg")` | `integrate_romberg` | Scalar callback; callers loop for vector-valued R semantics. |
| `int(..., type="TOMS614")` | `toms614_integrate` | Direct modern port of INTHP/TOMS 614. |
| `int2` | `integrate_2d` | Mapped Gauss-Legendre product quadrature. |
| `runge.kutta` | `runge_kutta` | Scalar RK4. |
| `mexp` | `matrix_exp` | Scaling/squaring Taylor implementation. |
| `lin.diff.eqn` | `lin_diff_eqn` | Uses `matrix_exp`. |
| `gauss.hermite` | `gauss_hermite` | Nodes and normalized weights. |
| `gettvc` / `gettvc_f` | `gettvc` | Numerical covariate alignment only. |
| `contr.mean` | `contrast_mean` | Numeric contrast matrix. |
| `capply(..., sum)` | `group_sum` | Integer group labels. |
| `capply(..., mean)` | `group_mean` | Integer group labels. |
| `mu1.*`, `mu2.*` | `mu1_*`, `mu2_*` | All 13 PK/PD mean functions. |

## Distribution mapping

Every exported upstream distribution has `d`, `p`, `q`, and `r` entry points
using the same base name: `invgauss`, `laplace`, `levy`, `pareto`, `simplex`,
`twosidedpower`, `boxcox`, `burr`, `gextval`, `ggamma`, `ginvgauss`, `glogis`,
`gweibull`, `hjorth`, `powexp`, `skewlaplace`, `betabinom`, `doublebinom`,
`multbinom`, `doublepois`, `multpois`, `pvfpois`, `gammacount`, and `consul`.

Example: R `qgammacount(p,m,s)` maps to Fortran
`qgammacount(p,m,s)`.
