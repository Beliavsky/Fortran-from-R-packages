# Translation coverage

## Directly translated / implemented

| Upstream functionality | Fortran status |
|---|---|
| `rq.fit.fnb` / `rqfnb.f` | Direct modern translation |
| `rq.fit.fnc` / `rqfnc.f` | Direct modern translation |
| `rq.fit.qfnb` / `qfnb.f` | Implemented using the translated FNB core |
| `rq.wfit` with FNB | Implemented |
| `rq.fit.lasso` | Implemented |
| `rq.fit.scad` | Implemented |
| `rq.fit.pfn` | Implemented at numeric-array level |
| `lprq` | Implemented |
| nonlinear QR | Iterative-linearization Fortran API |
| `kuantile`, `q489`, `qselect` | Implemented |
| XY bootstrap | Implemented |
| `rls` | Implemented |
| `combos` | Implemented |
| exponential random helper | Implemented |

The Frisch-Newton iteration keeps the upstream Mehrotra predictor-corrector
structure, primal/dual step restrictions and centering formula. BLAS/LAPACK
operations are expressed with Fortran array operations plus a self-contained
SPD Cholesky solve.

## Deliberately not represented as if translated

The following are important parts of the full R package but are **not yet**
part of the compiled v0.1.0 Fortran API:

- Barrodale-Roberts `rq.fit.br`, including complete tau solution paths and
  rank-inversion confidence intervals (`rqbr.f`).
- Sparse Frisch-Newton (`rq.fit.sfn`, `rq.fit.sfnc`) and the associated
  SparseM/sparse-Cholesky machinery (`srqfn*`, `cholesky.f`, `sparskit2.f`).
- `rqss` total-variation smoothing and its sparse penalty machinery.
- Censored quantile regression (`crq`) Powell/Portnoy/Peng-Huang kernels.
- Specialized bootstrap kernels `pwxy`, `pwy`, `wxy`, `xys` beyond the
  portable pairs bootstrap included here.
- Full inference/summary machinery, rank tests, Khmaladze/Pareto tests and
  extreme-value inference.
- `dynrq`, formula/model-frame processing and R time-series integration.
- `rq.fit.conquer`, whose algorithm belongs to the external `conquer` package.
- Portfolio/front-end routines that depend on R object ecosystems.

Those original sources are retained unchanged under `original/` for
provenance and for future translation work. No unrelated replacement algorithm
is exposed under an upstream method name.

## Representation differences

- R lists/classes are replaced by Fortran derived result types.
- Dense matrices are conventional Fortran `(n,p)` arrays; callers do not need
  to manually transpose them into the old `.Fortran` ABI layout.
- `kuantiles` sorts a working copy for multiple requested quantiles rather than
  using the original repeated selection kernel. Numerical definitions are the
  same Hyndman-Fan types.
- `nlrq_fit` exposes the numerical model/Jacobian callback directly instead of
  parsing an R formula or `numericDeriv` attributes.
