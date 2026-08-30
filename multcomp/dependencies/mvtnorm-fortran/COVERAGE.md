# Computational coverage

## Translated numerical areas

| Upstream area | Fortran coverage |
|---|---|
| `dmvnorm`, `rmvnorm` | Multivariate normal log densities, densities, and seeded simulation |
| `dmvt`, `rmvt` | Shifted multivariate-t density and shifted/Kshirsagar simulation |
| `pmvnorm`, `pmvt` | Normal, shifted-t, and Kshirsagar-t rectangular probabilities |
| `qmvnorm`, `qmvt` | Lower-, upper-, and two-sided simultaneous quantiles |
| `GenzBretz` | Randomized conditional-integration control with antithetic shifted Halton points |
| `TVPACK` | Exact univariate and deterministic bivariate normal route; high-accuracy common route otherwise |
| `Miwa` | Compatible control and deterministic common probability route |
| `ltMatrices`, `syMatrices` numerical methods | Packing, unpacking, multiplication, solution, cross-products, conversions, standardization, permutation, correlations and partial correlations |
| `chol2cov`, `cov2chol`, `invchol2chol`, `chol2invchol` | Included |
| `invchol2cov`, `cov2invchol`, `invchol2pre`, `chol2pre` | Included |
| `chol2cor`, `invchol2cor`, `chol2pc`, `invchol2pc` | Included |
| `Dchol`, `invcholD`, `vectrick` | Included as `dchol`, `invchold`, and `vectrick` |
| `marg_mvnorm`, `cond_mvnorm` | Arbitrary-index marginal and conditional Gaussian laws |
| `mvnorm` numerical behavior | Typed model construction, simulation, permutation, means and covariance through `mvnormal_model` |
| `ldmvnorm`, `lpmvnorm`, `ldpmvnorm` | Exact, interval-censored, and mixed-data log likelihoods |
| `sldmvnorm`, `slpmvnorm`, `sldpmvnorm` | Deterministic numerical scores for mean and Cholesky parameters |
| `deperma`, `destandardize` | Numerical Jacobian-based score transformations |
| `lpRR`, `slpRR` | Discrete random-effects probability and numerical score |
| Distribution support | Normal, Student-t, chi-square, incomplete beta and incomplete gamma routines |

## Intentional API redesign

- R S3 classes are replaced by typed Fortran structures and ordinary arrays.
- Observations are rows, following common Fortran numerical-library practice.
- Multiple lower-triangular matrices use rank-three arrays or caller loops
  rather than metadata-bearing packed R matrices.
- Character names, dimnames, formula dispatch, and R recycling rules are not
  reproduced.
- Numerical scores use finite differences with common deterministic seeds.

## Not compiled

- R printing, subsetting, dimension-name, and data-frame methods.
- R's C registration layer and external C API headers.
- Documentation-building and vignette infrastructure.
- Exact line-by-line reproductions of the historical Miwa and trivariate
  TVPACK implementations. Their method selectors use the common modern
  integration engine, with deterministic bivariate normal integration where
  applicable.
- The upstream specialized minimax-tilted C implementation used internally by
  `lpmvnorm`; the Fortran likelihood calls the public rectangle-probability
  engine instead.

The original implementations remain under `original/mvtnorm-1.4-2` for
comparison and provenance.
