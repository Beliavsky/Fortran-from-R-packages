# Porting notes

## Numerical representation

R S3/S4 classes, environments, attributes and caches are represented by typed
Fortran derived types and arrays. Status integers replace R conditions where a
recoverable numerical failure needs to be reported.

Rssa's FFTW-backed trajectory matrix layer is replaced with direct dense Hankel,
Toeplitz and block-Hankel construction/multiplication. The result is portable and
deterministic but can use substantially more memory and time for large problems.

## Decomposition

Real and complex decompositions reuse the sibling `svd` package rather than
copying Rssa's dependency code. The sibling complex `ztrlan_svd_result` returns
only singular values and left singular vectors, so CSSA reconstructs each right
singular vector as `A^H u_j / sigma_j`; numerically zero singular values map to
zero right vectors. Toeplitz SSA forms the lag-covariance matrix explicitly.
MSSA concatenates channel trajectory matrices. Rectangular 2-D SSA forms the
block-Hankel trajectory matrix directly.

General shaped/circular multidimensional window masks are not implemented.
Projection SSA accepts explicit orthonormal row and column projector bases in
place of R formula/list construction.

## Forecasting and parameter estimation

LRR, recurrent forecasting, vector forecasting and LS-ESPRIT are implemented
with typed group-index inputs. MSSA vector/recurrent forecasting currently uses
the common column-direction formulation. Projection-SSA recurrent/vector
corrections are translated with dense least-squares and direct Hankelization.
Bootstrap forecast intervals use the Fortran intrinsic RNG and empirical order
statistics, so they are not R-seed-identical.

## Missing data and grouping

Sequential and iterative gap filling support ordinary SSA, complex SSA and MSSA;
iterative filling also supports rectangular 2-D SSA. Iterations are explicitly
bounded in the Fortran API.

Automatic w-correlation grouping uses deterministic complete-link clustering.
Spectral grouping uses a direct DFT rather than R's FFT/periodogram stack.
Consequently tie handling and floating-point roundoff can differ while the
underlying grouping criterion is preserved.

## Weighted and oblique SSA

Diagonal weighted-oblique decomposition is available as `decompose_wossa`; it
uses the diagonal metric case implemented upstream and does not claim the
unimplemented non-diagonal metric branch. `owcor_ssa` and `wcor_ossa` provide
typed oblique and ordinary correlations for selected reconstructed groups.

`fossa_ssa` implements the real 1-D Filter-adjusted O-SSA calculation for a
single reversed FIR vector, including finite-gamma and filter-only modes.
MSSA/nD filter arrays, shaped masks, circular filtering, and simultaneous
multiple filters remain outside this partial mapping. `decompose_ossa` preserves
upstream's explicit refusal to continue an Oblique SSA decomposition by
returning `rssa_not_supported`.

`iossa_ssa` implements real 1-D iterative O-SSA. R's nested list of component
indices is represented by a selected-index vector plus one integer group label
per selected component. The dense reference implementation retains the upstream
kappa separation, left/right balance, low-rank oblique transform, Hankelization,
and RMS convergence logic. R cache mutation, printed diagnostics, and non-1-D
object families are omitted.

`eossa_ssa` implements the default real 1-D column-space LS EOSSA route: ESPRIT
shift solve, general complex eigenvectors, complete-link root clustering, SVD
realification of clustered eigenspaces, and oblique basis rotation. The root
clusters are labeled deterministically in frequency order; exact R `hclust`
cluster labels are not promised. Row-subspace, TLS, multidimensional/circular
mask, and beta-weighted multi-axis EOSSA paths remain omitted.

All 122 functions in the declared computational coverage basis now have a
meaningful mapping, but the package status remains `substantial` because the
partial mappings above are not complete option-level R compatibility.
