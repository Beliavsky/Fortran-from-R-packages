# Porting notes

## FBM representation

Upstream `bigstatsr` uses R reference classes plus memory mapping supplied by
its C++/rmio infrastructure. Standard Fortran has no portable memory-mapping
facility. The Fortran port therefore stores matrices as raw column-major
stream files and performs explicit column/block I/O.

This has two advantages: no R runtime and no POSIX-only system dependency. It
also means highly random element access is slower than an mmap-backed FBM.
The numerical routines are written to use column/block access rather than
individual element reads.

Version 0.1.0 supports real64 FBMs and the package's important 0..255 coded
matrix case. Upstream's float, signed integer, and unsigned-short storage modes
remain future storage-layer work; they do not change the statistical formulas.

## Indices

The Fortran API uses normal Fortran 1-based integer indices. No conversion from
R indices is necessary.

## `big_univLinReg`

The R wrapper computes an SVD of the intercept/covariate matrix, retains left
singular vectors whose singular values satisfy

`d / (sqrt(n) + sqrt(p) - 1) > thr.eigval`,

and applies the optimized Frisch-Waugh-Lovell formulas. The Fortran port uses
the same rank criterion and projection formulas.

## `big_univLogReg`

The Rcpp kernel uses IRLS and falls back to R's `glm()` for columns that do not
converge. The Fortran implementation uses a standalone IRLS/Newton solver for
every column. It reports convergence explicitly; there is no R `glm()` fallback.

## Sparse regression

The upstream C++ code is a modified version of `biglasso` coordinate descent,
with strong-set screening and an outer R CMSA/grid-search workflow. The port
implements the Gaussian and logistic elastic-net objectives directly with
warm-started coordinate-descent paths. It does not reproduce the R-only CMSA
or cluster-parallel model-selection layer in v0.1.0.

## `big_randomSVD`

Upstream calls `RSpectra::svds()` on function-based matrix products. The Fortran
port does the same computationally through the vendored Fortran RSpectra layer:
a derived `linear_operator` reads the FBM by columns and ARPACK performs the
partial eigensolver iterations.

## Parallelism

`bigstatsr` coordinates R workers and also uses OpenMP in selected C++ kernels.
The v0.1.0 algorithms are thread-safe at the numerical level but do not create
R-style worker clusters. BLAS/LAPACK/ARPACK may use whatever threading is
provided by the linked system libraries.

## AUC bootstrap

The sorted AUC and tabulated-bootstrap formulas follow upstream `src/AUC.cpp`,
including its tie convention. Random bootstrap indices use Fortran's
`random_number`.

## Plotting and R object infrastructure

Plotting methods, themes, plotly text construction, S3/S4 display methods,
RDS attach/save metadata, and data-frame factor expansion are intentionally
not translated because they are presentation/runtime infrastructure rather
than numerical kernels.
