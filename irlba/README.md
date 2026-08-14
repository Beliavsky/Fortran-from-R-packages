# irlba-fortran v0.1.0

Modern Fortran/FPM translation of the computational code in the R package
`irlba` 2.3.7.

The library implements fast truncated singular value decomposition using the
Baglama-Reichel implicitly restarted Lanczos bidiagonalization algorithm, plus
the other exported numerical algorithms in the package.

## Implemented functionality

- `irlba_svd`: real dense and native CSC sparse truncated SVD.
- Thick restart from a previous `irlba_result`.
- Implicit column centering and scaling without forming a centered dense matrix.
- Square-matrix scalar shifts.
- Full reorthogonalization and the upstream residual / singular-value-change
  convergence tests.
- Native CSC sparse matrix-vector and block matrix products; no R Matrix or
  CHOLMOD dependency is required.
- `partial_eigen`: partial eigenvalue decomposition for real symmetric dense or
  CSC matrices using the same shift strategy as upstream.
- `prcomp_irlba`: truncated PCA with optional centering/scaling, scores,
  standard deviations, and explained-variance summaries.
- `ssvd`: Shen-Huang style L1-thresholded sparse SVD/PCA iteration.
- `svdr`: randomized block truncated SVD for dense and CSC sparse matrices.
- `irlba_complex`: complex dense SVD compatibility path.
- Public `linear_operator` abstraction and `irlb_operator` entry point for
  matrix-free/custom matrix-vector products.

## Deliberate implementation boundaries

The upstream package has two implementations of `irlba`: a fast C path for
real dense/`dgCMatrix` matrices and a slower R reference path that also handles
complex matrices, `right_only`, smallest singular values, and arbitrary R
matrix classes. This port translates the real fast path natively and exposes a
Fortran matrix-operator interface for custom products.

For correctness in v0.1.0:

- `smallest=.true.` uses a LAPACK full-SVD fallback. For CSC input this fallback
  materializes the matrix, so it is not suitable for huge sparse problems.
- complex-valued input uses LAPACK `ZGESSD`/`ZGESDD`-style full SVD rather than
  an iterative complex IRLB path.
- the R `right_only` memory-saving mode is represented functionally by choosing
  the desired returned rank; the internal real IRLB basis remains two-sided.
- R S3/S4 classes, formula/data-frame handling, messages, plotting, and CHOLMOD
  object adapters are intentionally omitted.

## FPM

The package requires BLAS and LAPACK:

```toml
[build]
link = ["lapack", "blas"]

[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
```

Typical use:

```fortran
use irlba

type(irlba_result) :: s
type(irlba_control) :: ctl
real(dp) :: a(100, 30)

ctl%tol = 1.0e-8_dp
s = irlba_svd(a, 5, control=ctl)
```

## License

GPL-3.0-or-later, following the supplied upstream package. See `LICENSE` and
`UPSTREAM.md`.
