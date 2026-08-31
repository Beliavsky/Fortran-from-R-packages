# API coverage

This translation targets `ordinal`'s reusable numerical API rather than R's
formula, S3, plotting, and model-frame layers. Callers supply already-built
numeric design matrices and integer ordered-response codes.

## Implemented

### Link and threshold kernels

- All seven cumulative-link families used by current `ordinal`: logit, probit,
  cloglog, loglog, cauchit, Aranda-Ordaz, and log-gamma.
- CDF, density, and density-gradient kernels corresponding to the upstream
  generalized-normal/logistic/Cauchy, Gumbel, Aranda-Ordaz, and log-gamma code.
- Flexible, symmetric, symmetric2, and equidistant threshold structures,
  including threshold Jacobians and upstream-style starting values.

### Cumulative-link models (`clm`)

- Weighted cumulative-link likelihood for supplied numeric location matrices.
- Location offsets, log-scale submodels and scale offsets.
- Nominal-effects submodels using observation-specific threshold shifts with
  the same threshold-basis construction as upstream `clm`.
- Fixed or estimated Aranda-Ordaz/log-gamma shape parameters.
- Exact fixed-shape CLM gradient and Hessian translated from the upstream
  Newton fitting formulas.
- Analytic Newton fitting with step halving for fixed-shape links; numerical
  BFGS is retained for estimated flexible-link shape parameters.
- Covariance matrix and Hessian diagnostics: gradient norm, eigenvalue range,
  estimated rank, condition measure, and positive-definiteness flag.
- Observed-category fitted probabilities and full category-probability
  prediction, including nominal effects.
- Profile likelihood for fixed-shape CLM coefficients using analytic
  constrained Newton nuisance refits.
- Profile-likelihood and normal-theory Wald confidence intervals.
- Rank detection and rank-deficient-column reduction utilities for numeric
  design matrices.

### Cumulative-link mixed models (`clmm`)

For a single scalar Gaussian random-intercept term:

- Laplace likelihood for `nAGQ = 0` or `1`.
- Adaptive Gauss-Hermite quadrature for `nAGQ > 1`.
- Non-adaptive Gauss-Hermite quadrature for `nAGQ < 0`.
- Generated Gauss-Hermite rules rather than a single hard-coded rule.
- Conditional random-effect modes and conditional variances.
- Random-effect standard-deviation estimation and outer covariance/Hessian
  diagnostics.

For general Gaussian random effects under Laplace approximation:

- Multiple random-effect terms, including crossed or nested grouping layouts
  supplied by group-index columns.
- Random slopes / multiple random coefficients per group.
- Full correlated covariance matrices within each random-effect term,
  parameterized by a packed lower Cholesky factor with log diagonal entries.
- Joint conditional modes and conditional covariance matrix.
- Marginal Laplace likelihood, fitting, numerical outer Hessian/covariance,
  and convergence diagnostics.

The division above mirrors an important upstream limitation: adaptive
Gauss-Hermite quadrature is useful for a single scalar random-effect term,
whereas general multi-term/random-slope models use the Laplace path.

## Intentionally not translated because it is R-specific

- Formula/model-frame construction, factor contrasts, environment updates,
  and automatic S3 dispatch.
- Printing, summaries, plotting, profile/slice visualization, and plotting
  interpolation helpers.
- `Matrix` sparse-object classes and compiled R registration glue.
- R serialization and interactive utilities.

## Remaining differences

The main remaining differences are interface/performance details rather than
missing major statistical model classes:

- General multi-term/random-slope Laplace calculations currently use dense
  matrices. Upstream R code can exploit sparse `Matrix` structures for large
  problems.
- Automatic alias handling is not embedded in every fitter. The public
  rank/drop utilities let callers reduce numeric design matrices explicitly,
  which is the natural equivalent when no R formula/model-matrix layer exists.
- Flexible-link shape fits use deterministic numerical BFGS derivatives for
  the outer shape parameter, matching the need for numerical treatment but
  not every internal control/detail of R's optimizer stack.
- Profile calculations return numerical likelihood/interval results directly;
  R-specific spline/plot objects and profile graphics are intentionally absent.
- Legacy `clm2`/`clmm2` compatibility interfaces are not reproduced because
  their useful numerical functionality is represented by the current CLM/CLMM
  kernels.

These differences are not replaced by R dependencies, C wrappers, or vendored
translated packages.
