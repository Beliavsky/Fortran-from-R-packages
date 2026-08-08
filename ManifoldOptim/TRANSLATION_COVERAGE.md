# Translation coverage

## Upstream package

The supplied source is ManifoldOptim 1.0.2. Its DESCRIPTION identifies it as
an R interface to ROPTLIB 0.3 and declares `GPL (>= 2)`. The complete supplied
source tree is retained under `original/ManifoldOptim-master`.

## Computational surface

The Fortran package covers the computational functionality exposed by the R
interface:

- Euclidean, Sphere, Stiefel, Grassmann, SPD, LowRank, and OrthGroup domains
- product domains and `numofmani` multiplicity
- tangent projection, retraction, differentiated retraction, vector transport,
  inverse transport/cotangent operations, metrics, feasible random points, and
  Euclidean-to-Riemannian derivative conversion
- objective-only numerical gradients and numerical Riemannian Hessian actions
- all eleven solver names exposed by `manifold.optim`
- `orthonorm` functionality
- solver histories and evaluation counters

## Architectural gaps closed in 0.2.0

Compared with 0.1.0, this release removes the main computational shortcuts:

1. **Stiefel ParamSet=2**
   - Implements the constructed retraction using the skew block representation
     and matrix exponential rather than silently using QF.
   - ParamSet=1 continues to use QF.

2. **Sphere parameter sets**
   - ParamSet=2 uses the exponential retraction.
   - ParamSets 3 and 4 use QF with the sphere parallel-translation transport.
   - ParamSet=4 computes the differentiated-retraction locking `Beta` scaling.
   - ParamSet=1 uses intrinsic-coordinate parallelization semantics.

3. **Grassmann / Stiefel intrinsic parallelization**
   - Flat extrinsic vectors are converted through deterministic orthogonal
     complements so that transport copies intrinsic coordinates rather than
     merely projecting the old extrinsic vector at the new point.

4. **LowRank quotient geometry**
   - `U,D,V` storage is retained at the public boundary.
   - Vertical Stiefel gauge components are absorbed into the `D` component as
     in upstream `LowRank::ObtainIntr`/`ObtainExtr`.
   - The metric includes the upstream `D` and `D^T` scaling of the U/V tangent
     blocks.
   - Vector transport preserves the scaled intrinsic U/V coordinates and
     rescales them using the new `D` factor.

5. **LRTRSR1**
   - Uses limited-memory S/Y pairs and the compact SR1 action
     `(Y-gamma*S) (SY-gamma*SS)^(-1) (Y-gamma*S)^T` plus `gamma I`.
   - It no longer allocates or updates a dense n-by-n SR1 matrix.

6. **Line search**
   - Armijo, weak Wolfe, strong Wolfe, and exact scalar search modes are
     available through `solver_options%line_search`.
   - `LINESEARCH_INPUTFUN` and `solver_options%line_search_proc` provide the
     upstream-style user step-selection hook.

7. **Quasi-Newton details**
   - RBFGS/LRBFGS use manifold `Beta` scaling in the secant gradient vector.
   - ROPTLIB-style `nu`/`mu` curvature safeguards are exposed as `qn_nu` and
     `qn_mu`.
   - `RBroydenFamily` implements the Broyden-family rank-one term. The upstream
     default Phi=1 remains the default; `broyden_phi` additionally permits
     other family members.
   - RWRBFGS performs its update at the old point and then transports the
     inverse-Hessian approximation.

8. **RCG variants**
   - Fletcher-Reeves, Polak-Ribiere, PR+, FR-PR, Hestenes-Stiefel,
     Dai-Yuan, and Hager-Zhang beta formulas are supported.
   - Hestenes-Stiefel (`"HS"`) is the default, matching upstream RCG.

9. **Euclidean Hessian conversion**
   - Stiefel/OrthGroup and Grassmann analytical Hessian-vector callbacks now
     receive the connection correction involving the Euclidean gradient before
     tangent projection, following upstream `EucHvToHv` logic.

10. **Cotangent actions used by RWRBFGS**
    - QF Stiefel/OrthGroup, QF Grassmann, and exponential/QF Sphere use explicit
      cotangent formulas rather than treating cotangent transport as ordinary
      projection.

## Remaining implementation differences

The following differences are intentional or are not fully eliminated:

- The public API remains flat Fortran arrays rather than reproducing ROPTLIB's
  C++ `Variable`/`Vector`/`SharedSpace` class hierarchy. This is an API/design
  difference, not a missing numerical method.
- Intrinsic Stiefel/Grassmann parallelization uses a deterministic orthogonal
  complement generated in Fortran rather than ROPTLIB's cached LAPACK
  Householder reflector representation. It represents the same tangent-space
  decomposition but will not be bit-for-bit identical to ROPTLIB coordinates.
- ROPTLIB's specialized LowRank `coTangentVector` routine is not copied in full.
  For RWRBFGS on LowRank, the Fortran implementation uses the inverse
  intrinsic-coordinate transport. All other translated LowRank geometry uses
  the quotient/scaling representation described above.
- ROPTLIB contains more elaborate interpolation and initial-step heuristics
  inside each line-search algorithm. The Fortran implementation supplies the
  same four algorithm classes and conditions but uses a smaller standalone
  bracketing/interpolation implementation.
- Temporary-data caches and LAPACK/BLAS acceleration are replaced by allocatable
  Fortran work arrays and self-contained dense linear algebra. This can affect
  performance and last-bit floating-point results, not the exposed algorithms.
- Rcpp/R registration, printing, timers, debugging infrastructure, MATLAB/mex
  compatibility, plotting/presentation code, and R documentation plumbing are
  intentionally not translated.

These remaining differences are documented so this package does not claim
bit-for-bit identity with ROPTLIB. Computationally, the major 0.1.0 shortcuts
(Stiefel ParamSet=2, LowRank quotient scaling, dense LRTRSR1, Armijo-only line
search, and missing locking/Beta behavior) are removed in 0.2.0.
