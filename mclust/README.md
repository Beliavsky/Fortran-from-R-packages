# mclust-fortran v0.1.2

Modern Fortran/FPM computational port of the R package **mclust 6.1.3**
(Gaussian finite-mixture modelling for clustering, classification, and
density estimation).

## What is included

The port keeps the mature upstream MCLUST covariance-estimation kernels and
places them behind a modern free-form Fortran API with explicit interfaces,
derived result/control types, allocatable arrays, and FPM packaging.  The
high-level R numerical workflows needed to use those kernels without R have
been reimplemented in Fortran.

Implemented functionality includes:

- one-dimensional `E` and `V` Gaussian mixtures;
- all 14 multivariate covariance models: `EII`, `VII`, `EEI`, `VEI`, `EVI`,
  `VVI`, `EEE`, `EVE`, `VEE`, `VVE`, `EEV`, `VEV`, `EVV`, and `VVV`;
- EM fitting with a common modern E-step and the original constrained M-step
  algorithms;
- hierarchical agglomerative initialization (`E`, `V`, `EII`, `VII`, `EEE`,
  `VVV`) and the upstream `VARS`, `STD`, `PCS`, `PCR`, `SVD`, and `SPH`
  transformations;
- BIC/ICL parameter counting, model-grid fitting, and model selection;
- mixture/component density evaluation, posterior probabilities, MAP classes,
  uncertainty, 1-D CDF/quantiles, and HDR levels;
- simulation from a fitted mixture;
- model-based discriminant analysis and semi-supervised classification;
- MclustDR dimension reduction and projection;
- conditional-Gaussian missing-value imputation;
- entropy-based cluster combination;
- weighted EM fitting;
- parametric bootstrap likelihood-ratio and parameter bootstrap helpers;
- covariance-weighted statistics, adjusted Rand index, Brier score,
  `logsumexp`, `softmax`, counting, random orthogonal matrices;
- `unmap`, majority-vote and optimal cluster-label matching utilities;
- canonical discriminant coordinates (`crimcoords`).

Plotting, printing, R formula/data-frame/S3 plumbing, and R graphics helpers are
not ported.

## Important v0.1.2 scope boundary

This is a computational-core port, not yet a bit-for-bit translation of every
optional high-level mclust branch.  The following computational features are
not implemented in v0.1.2:

- conjugate-prior/Bayesian regularization (`priorControl`, `defaultPrior`, and
  the legacy `*p` fitting paths);
- the optional uniform noise/background component used by some `mclustBIC`
  workflows;
- `gmmhd` high-density connected-component clustering;
- `MclustDRsubsel` variable/direction subset selection;
- the full `cvMclustDA` cross-validation wrapper;
- `imputePairs`, `hcRandomPairs`, and a few small R partition/object helpers.

The underlying upstream source routines are retained where applicable, so
these are natural targets for a later compatibility release.

## FPM legacy-language compatibility

The translated legacy numerical kernels are free-form `.f90` sources, but they intentionally retain valid historical Fortran features such as implicit external procedure interfaces, implicit typing in some original routines, labeled `DO`, arithmetic `IF`, and computed `GO TO`.  Modern fpm releases disable implicit typing and implicit external interfaces by default, so `fpm.toml` explicitly enables those two language features while keeping `source-form = "free"`.


## Build with FPM

BLAS and LAPACK are required.  The manifest links `lapack` and `blas`.

```text
fpm build
fpm test
fpm run --example basic_mclust
```

The three upstream numerical kernel files are now ordinary free-form `.f90` sources. This is a source-form conversion only: valid legacy constructs such as labeled `DO`, `COMMON`, computed `GO TO`, and arithmetic `IF` are not rewritten merely for style.

The public umbrella module is `mclust`:

```fortran
program demo
  use mclust
  implicit none
  real(dp) :: x(100,2)
  type(mclust_fit) :: fit
  integer :: status

  ! Fill x here.
  call mclust_select(x,fit,status=status)
  if(status /= 0) error stop 'mclust fit failed'

  print *, fit%model_name, fit%g, fit%bic
end program demo
```

## Source organization

`src/*.f90` contains the modern Fortran API and translated high-level
algorithms. `src/legacy/*.f90` contains the numerically mature upstream MCLUST
kernels converted from fixed source form to free source form, retained to avoid replacing the specialized
Celeux-Govaert covariance updates with approximations.  These routines are
called only through explicit interfaces in `mclust_legacy_interfaces.f90`.

This hybrid approach modernizes normal application-facing use while preserving
upstream numerical algorithms and provenance.

See `API_MAPPING.md`, `ALGORITHM_NOTES.md`, `VALIDATION.md`, and `UPSTREAM.md`.
