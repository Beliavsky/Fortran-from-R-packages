# Porting notes

## Scope

All five exported computational functions are represented:

| R function | Fortran |
|---|---|
| `CLA` | `critical_line`; compatibility generic `CLA` |
| `MS` | `mean_sigma`; compatibility `MS` |
| `findSig` | `find_sigma`; compatibility `findSig` |
| `findMu` | `find_mu`; compatibility `findMu` |
| `muSigmaGarch` | `mu_sigma_garch`; compatibility `muSigmaGarch` |

The internal kernels `initAlgo`, `getMatrices`, `computeInv`, `computeW`, and
`computeLambda` are translated into private typed procedures in `cla_core`.
The historical `Ver8`, `Ver9`, and environment copies are retained in
`original/CLA` but are obsolete alternative implementations rather than
separate exported algorithms.

## Numerical substitutions

- R's `solve` is replaced by LAPACK `DGESV`.
- Positive-definiteness checking uses LAPACK `DPOTRF`.
- R's `uniroot` in `findMu` is replaced by bracketed bisection.
- R's `fGarch::garchFit` dependency is replaced by a self-contained
  constant-mean GARCH(p,q) maximum-likelihood estimator using Nelder-Mead.
- Normal and standardized Student-t innovations are supported. The original
  default `cond.dist="std"` maps to standardized Student-t.

The original `formula` argument could delegate arbitrary formulas to
`fGarch`. The typed Fortran interface instead accepts explicit ARCH and GARCH
orders. Mean ARMA terms and other `fGarch` distribution names are not part of
CLA's own algorithms and are not emulated.

## Stabilized turning-point filtering

The R source repeatedly halves `tol.lambda` almost to floating-point underflow
when a candidate lambda is numerically equal to the current one. That behavior
can generate platform-dependent duplicate turning points. The Fortran port
stops tolerance reduction at `100*epsilon`, rejecting numerically identical
cycles while retaining genuinely smaller lambdas.

## Presentation differences

S3 classes, row names, printing, plotting, R expressions, and data-frame
presentation are omitted. `free_indices` is represented by a logical
`free_mask` matrix, which is simpler and safer for compiled use.
