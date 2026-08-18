# Computational coverage

## Translated

### Core families

| Upstream area | Fortran coverage |
|---|---|
| `indepCopula` | Arbitrary-dimensional CDF, density, simulation |
| `normalCopula` | Multivariate density/simulation; bivariate deterministic CDF; higher-dimensional conditional integration |
| `tCopula` | Multivariate density/simulation; bivariate deterministic CDF; higher-dimensional scale-mixture integration |
| `claytonCopula` | Arbitrary-dimensional CDF and frailty simulation; bivariate density/fitting/transforms |
| `gumbelCopula` | Arbitrary-dimensional CDF and positive-stable frailty simulation; bivariate density/fitting/transforms |
| `frankCopula` | Arbitrary-dimensional positive-parameter CDF/simulation and bivariate negative-parameter support |
| `amhCopula`, `joeCopula` | Bivariate CDF, density, simulation, fitting |
| `fgmCopula`, `plackettCopula` | Bivariate CDF, density, simulation, fitting |
| `moCopula` | Bivariate CDF, tau, tail dependence; ordinary density marked singular |
| `lowfhCopula`, `upfhCopula` | Bivariate bounds, dependence measures, singular-density convention |
| Galambos, Husler-Reiss, Tawn | Pickands function, CDF, numerical density, simulation and selected fitting/dependence operations |

### Transformations and constructions

- 90-, 180-, and 270-degree rotations.
- Weighted mixtures of typed copula models.
- Khoudraji asymmetrization for bivariate base copulas.
- Two-level nested-Clayton CDF with parameter-order validation.
- Bivariate conditional CDF, Rosenblatt transform, and inverse transform.

### Dependence and empirical methods

- Model Kendall tau, Spearman rho, and lower/upper tail dependence.
- Parameter inversion from Kendall tau and Spearman rho.
- Sample Kendall tau-b and Spearman rho.
- Pseudo-observations with averaged ties.
- Empirical copula at arbitrary points.
- Permutation tests for multivariate independence, bivariate exchangeability,
  and radial symmetry.

### Estimation

- Inversion-of-Kendall-tau estimation.
- One-parameter maximum pseudo-likelihood using bounded golden-section search.
- Numerical observed-information standard errors.
- Gaussian and Student correlation fitting, with Student degrees of freedom
  treated as fixed.

### Special numerical tools

- Normal, Student, incomplete-beta, and incomplete-gamma functions.
- Cholesky decomposition, matrix inversion, nearest correlation matrices, and
  sample correlations.
- Stirling numbers of both kinds, Eulerian numbers, Sibuya laws, logarithmic
  series laws, positive-stable simulation, and seeded xorshift random numbers.

## Not translated

The following upstream areas remain outside this release:

- The complete S4 class hierarchy, method dispatch, fixed-parameter objects,
  data frames, formula handling, and R object serialization.
- Plotting, lattice graphics, pair plots, contours, perspective plots, and
  color/annotation infrastructure.
- General `mvdc` models with arbitrary R marginal-distribution callbacks and
  their joint fitting infrastructure.
- Arbitrary nested Archimedean trees, all high-order generator derivatives,
  and every specialized nested-density optimization.
- The complete empirical-copula smoothing catalogue and all Pickands
  nonparametric estimators.
- The upstream serial-independence, multivariate-independence, extreme-value,
  exchangeability, radial-symmetry, goodness-of-fit, multiplier-bootstrap,
  cross-validation, and model-selection test families. The port provides
  smaller dependency-free permutation tests rather than claiming equivalence.
- Full parameter-score arrays, automatic differentiation, every derivative
  with respect to copula parameters, and all sandwich covariance variants.
- Specialized stable-distribution interfaces delegated upstream to `stabledist`,
  GSL, `mvtnorm`, `numDeriv`, `pcaPP`, `pspline`, and other R packages.
- High-dimensional Miwa/TVPACK-style exact probability algorithms.
- Vignettes, demos, data sets, and R package installation infrastructure as
  compiled inputs. They remain retained for provenance.

## Mapping philosophy

The public Fortran module uses familiar names such as `pCopula`, `dCopula`,
`rCopula`, `cCopula`, `tau`, `rho`, `lambda`, `iTau`, `iRho`, and `pobs`, but
uses typed `copula_model` objects instead of R S4 instances.
