# Porting notes

Upstream contains roughly 10,500 lines across more than 200 R files. Version 0.3.0 extends the reusable numerical nucleus substantially.

## v0.2.0 additions

- Circular MLEs: cardioid and circular exponential.
- vMF MLE plus vMF mixture EM.
- Spherical Cauchy and PKBD estimators.
- Directional k-NN and vMF discriminant analysis.
- Rayleigh and Kuiper tests, including simulation-based p-values.
- Fisher-Bingham saddlepoint normalization and Kent density.
- Bingham, Kent, and matrix-Fisher simulation; matrix-Fisher sample-mean SVD.
- Internal Gaussian elimination, Jacobi symmetric eigensolver, and 3x3 SVD support.

The Kent random generator is distribution-compatible but uses a simpler uniform-sphere rejection envelope instead of the upstream Fisher-Bingham proposal. The PKBD optimizer maximizes the same mesos-parameter likelihood but uses a safeguarded numerical-gradient search rather than the upstream analytic Hessian Newton step.

Still not full parity: arbitrary-dimensional ESAG, ESAG/IAG/Kent/GCPC/SPML/SIPC/CIPC regression families, Kent MLE, mixture PKBD/spherical-Cauchy models, the large specialized ANOVA/bootstrap/permutation family, distance-correlation variants, KNN tuning/regression, and plotting/map helpers. Plotting/geographic presentation code remains intentionally out of scope.


## v0.3.0 additions

- General d-dimensional ESAG parameter mapping, density, and simulation (3-D remains available as `desag3`).
- PKBD rejection simulation and complete spherical-Cauchy/PKBD mixture density, simulation, and EM fitting workflows.
- Kent MLE with the upstream orientation construction and Fisher-Bingham saddlepoint normalizer; optimization uses a constrained pattern search enforcing `0 <= 2 beta < kappa`.
- Generic distance correlation plus circular and spherical wrappers.
- Multivariate KNN regression and deterministic K-fold tuning, including spherical response normalization/loss.

The remaining gaps are now dominated by large specialized regression families (ESAG/IAG/GCPC/SPML/SIPC/CIPC), specialized two-sample/ANOVA bootstrap-permutation procedures, and plotting/geographic presentation helpers. Those regression/test families require a broader optimization/design-matrix layer and are intentionally left for a later parity release rather than adding shallow wrappers.


## v0.4.0 additions

- General-dimensional isotropic angular Gaussian MLE.
- SIPC, CIPC, and GCPC MLEs.
- Embedding, high-concentration, and heterogeneous directional two-sample permutation procedures.
- New `directional_advanced` and `directional_permutation` modules, with `test_parity_v04`.

The optimization strategy is intentionally self-contained: where upstream R delegates to `optim`, `nlm`, or Rfast, the Fortran port uses deterministic Newton or pattern-search iterations. The largest remaining worthwhile gaps are the regression families (`esag.reg`, `iag.reg`, `spml.reg`, `sipc.reg`, `cipc.reg`, `gcpc.reg` and related model-selection helpers), plus some specialized LR/bootstrap ANOVA variants. Plotting/geographic presentation remains out of scope.

## v0.5.0 parity work

The principal standalone numerical gaps identified after v0.4.0 were the regression families and general-group directional ANOVA. v0.5.0 adds those families using the upstream Directional likelihood formulas while replacing R's `optim`, `nlm`, and the still-untranslated Rfast `spml.reg` backend with deterministic, self-contained Fortran optimization.

Regression APIs take an explicit numeric design matrix. Unlike R's `model.matrix`, the Fortran routines do not silently add an intercept; include a column of ones when an intercept is wanted. Responses are unit-vector matrices (including two-column cosine/sine responses for circular models).

`esag_reg` and `sespc_reg` implement the 3-D Directional package models. Their upstream orientation-grid wrapper is not reproduced: the joint likelihood is optimized directly from an OLS directional start. This avoids a large coarse grid search and keeps the Fortran API deterministic.

The ANOVA routines support arbitrary integer group labels from 1 through `maxval(group)`. Their optional `mc_reps` path returns a Monte Carlo permutation p-value and avoids requiring a separate special-function implementation of the F distribution. Two-sample centered bootstrap forms are provided for the embedding, high-concentration, and heterogeneous tests.

The v0.3 runtime-check warnings caused by noncontiguous row slices passed to mixture update routines and quaternion conversion were removed by explicit contiguous temporaries. All v0.1-v0.5 tests now run under `-fcheck=all` without those warnings.

Remaining upstream-only code is primarily plotting/contours/maps, datasets and I/O wrappers, formula/data-frame conveniences, model-selection wrappers around already translated kernels, and narrowly specialized duplicate circular ANOVA wrappers. These are not considered worthwhile numerical parity gaps for the standalone Fortran library.
