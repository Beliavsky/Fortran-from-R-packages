# API coverage

This file maps the computational portions of mc2d 0.2.2 to the Fortran port.
R vectorization is generally represented by scalar numerical functions plus
Fortran arrays/loops, while R's dynamic dispatch is represented by generic
interfaces, derived types, and procedure callbacks.

## Distributions

| R family | Fortran implementation |
|---|---|
| `d/p/q/rbern` | `dbern`, `pbern`, `qbern`, `rbern` |
| `d/p/q/rbetagen` | `dbetagen`, `pbetagen`, `qbetagen`, `rbetagen` |
| `d/p/q/rbetasubj` | `dbetasubj`, `pbetasubj`, `qbetasubj`, `rbetasubj` |
| `d/p/q/rlnormb` | `dlnormb`, `plnormb`, `qlnormb`, `rlnormb`; upstream defaults preserved |
| `d/p/q/rpert` | `dpert`, `ppert`, `qpert`, `rpert`; mean form via `*_mean` |
| `d/p/q/rtriang` | `dtriang`, `ptriang`, `qtriang`, `rtriang`; mean form via `*_mean` |
| `d/p/q/rmqi` | `dmqi`, `pmqi`, `qmqi`, `rmqi` |
| `d/p/q/rempiricalD` | `dempiricald`, `pempiricald`, `qempiricald`, `rempiricald` |
| `d/p/q/rempiricalC` | `dempiricalc`, `pempiricalc`, `qempiricalc`, `rempiricalc` |
| `d/rdirichlet` | `ddirichlet`, `rdirichlet` |
| `d/rmultinomial` | `dmultinomial`, `rmultinomial` |
| `d/rmultinormal` | `dmultinormal`, `rmultinormal`, plus varying-parameter forms |
| `rtrunc` | `rtrunc` using CDF/quantile procedure callbacks |
| `lhs` | `lhs` using a quantile procedure callback |

## Two-dimensional Monte Carlo framework

| R function/concept | Fortran implementation |
|---|---|
| `mcnode` / `mcdata` | derived type `mcnode`, generic `mcdata` |
| `ndvar`, `ndunc` | `ndvar`, `ndunc` |
| `typemcnode`, dimensions | `typemcnode`, `dimmcnode`, `dimmc` |
| `mc` | derived type `mc`, constructor `make_mc` |
| `Ops.mcnode` | overloaded `+`, `-`, `*`, `/`, `**` with broadcasting |
| `pmin`, `pmax` | `pmin_node`, `pmax_node` |
| `extractvar`, `addvar` | same names |
| `unmc`, `outm` | `unmc`, `outm_set` |
| `mcstoc` | callback-driven `mcstoc` |
| `mcprobtree` | `mcprobtree_weights`, `mcprobtree_switch` |
| `cornode` | generic `cornode`, matrix and node forms |
| `mcapply` | `mcapply_reduce`, `mcapply_elemental` |
| `summary.mcnode` | `node_summary`; `node_summary_each` for multivariate `outm="each"` |
| `quantile.mcnode` | `node_quantiles`; `node_quantiles_each` for multivariate `outm="each"` |
| `mcratio` | `mcratio` |
| `tornado` | `tornado`, including multivariate `outm="each"` outputs |
| `tornadounc` | `tornadounc`, including mean/SD/quantile statistics for VU nodes |
| `converg` computation | `running_convergence` (plotting excluded) |
| `mcmodel`, `evalmcmod` | type `mcmodel`, `evalmcmod` with callback |
| `mcmodelcut`, `evalmccut` | type `mcmodelcut`; full callback via `evalmccut`, low-memory column/reducer loop via `evalmccut_reduce` |

## Intentional adaptations / omissions

- Plotting and graphics-only files (`plot.*`, `hist.*`, `gghist`, `ggplotmc`,
  `ggspaghetti`, `ggtornado`, `spaghetti`) are excluded as requested.
- R S3 printing and dynamic method lookup are not reproduced; Fortran callers
  access result derived types directly.
- R expression capture (`substitute`, `eval`, environment lookup, dynamic
  function names) is replaced by explicit procedure callbacks. This preserves
  the computational capability without embedding an R interpreter.
- R logical-node predicates such as `is.na.mcnode` map naturally to intrinsic
  Fortran/`ieee_arithmetic` operations on `node%value`; no artificial numeric
  mcnode wrapper is introduced for logical results.
- R's automatic vector recycling is limited to the explicit mcnode broadcasting
  rules and varying-parameter routines; scalar distribution functions are
  intentionally idiomatic scalar Fortran functions.

- R permits `outm` to contain a vector of arbitrary function names. The Fortran
  derived type stores one `outm` string. `"each"` and `"none"` are represented
  directly; arbitrary R function-name reducers are not dynamically looked up.
  For multivariate tornado calculations, use `outm="each"` or reduce the
  variates explicitly in Fortran before calling the statistic.
- The low-memory `evalmccut` expression/list machinery is represented by
  `evalmccut_reduce`: a setup callback, a per-uncertainty-column model callback,
  and a typed reducer returning a fixed-size real statistic vector. This retains
  the looped computation without attempting to serialize arbitrary nested R
  lists or captured expressions.
