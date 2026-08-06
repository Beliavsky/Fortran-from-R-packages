# Porting notes

## Scope

`mgcv` 1.9-4 contains roughly two megabytes of R and C numerical code plus a
large R modelling framework. Reproducing its complete formula language, object
system, sparse model representation, specialized families, plotting, and
third-party integrations would require an R-compatible runtime. Version 0.1.0
therefore provides a native, array-oriented numerical core rather than a false
claim of drop-in parity.

## Data and model representation

R formulas and S3 smooth objects are replaced by:

- explicit design matrices,
- explicit three-dimensional penalty arrays `S(p,p,m)`,
- `smooth_spec_t` metadata for basis prediction, and
- `gam_model_t` fitted-model results.

The caller controls column placement. `embed_penalty` places a term penalty in
the corresponding coefficient block.

## Smoothing-parameter optimization

R `mgcv` has sophisticated Newton and quasi-Newton algorithms using first and
second derivatives of REML, ML, GCV, UBRE, and NCV criteria. This port uses a
bounded deterministic coordinate search in log smoothing parameters. GCV and
UBRE use standard dense effective-df calculations. `method_reml` is explicitly
reported as **REML-like** because it uses a dense log-determinant criterion and
omits some constants and derivative machinery of the R implementation.

## Smooth bases

- Cubic regression splines use the translated `splines` natural cubic basis.
  Their penalty is obtained by numerical integration of squared second
  derivatives on a fine grid.
- P-splines use a B-spline basis and a configurable finite-difference penalty.
- Cyclic smooths use a Fourier basis. This is periodic and well-conditioned,
  but it is not byte-for-byte the cyclic cubic regression basis in R `mgcv`.
- The 1D and 2D thin-plate-style constructors use low-rank radial functions
  with an identity penalty on radial coefficients. They preserve the intended
  smooth/null-space split but do not reproduce `tprs.c`'s exact truncation and
  reparameterization.
- Tensor products currently accept two margins. Their design is the row-wise
  Kronecker product and their two penalties are the standard Kronecker margin
  penalties.

## Penalized IRLS

The solver forms dense `X'WX + sum(lambda_j S_j)` systems and applies a small
scale-aware diagonal stabilization before Cholesky solution. Step halving is
used when the penalized deviance increases. This is intended for small and
medium dense models, not `bam`-scale data.

## Constraints

`pcls_fit` starts from the unconstrained penalized solution and alternates a
gradient step with projections onto equality and inequality constraints. This
is practical for moderate problems but is not the original active-set C
implementation in `qp.c`.

## Probability routines

`weighted_chisq_cdf` uses numerical inversion of the characteristic function.
It is a portable approximation to `psum.chisq`/Davies-style calculations and
allows the integration grid to be increased for difficult tails.

`rtweedie` implements the compound Poisson-gamma representation for
`1 < p < 2`. It does not cover all limiting cases handled by R `mgcv`.

## Numerical dependencies

No BLAS or LAPACK library is required. Dense Cholesky, Gaussian elimination,
and symmetric Jacobi eigensystems are implemented in Fortran. This maximizes
portability but will be slower than tuned vendor libraries on large matrices.
