# mvtnorm-fortran

A dependency-free modern Fortran/FPM numerical port of `mvtnorm` 1.4-2.

The library provides multivariate normal and Student-t densities, simulation,
rectangular probabilities, simultaneous quantiles, lower-triangular matrix
parameterizations, Gaussian marginal and conditional distributions, and
likelihoods for exact, interval-censored, and mixed observations.

## Build

```text
fpm build
fpm test
fpm run mvtnorm_demo
fpm run --example distributions_and_probabilities
fpm run --example matrices_and_likelihoods
```

The source uses Fortran 2018, `implicit none`, double precision defined by
`dp = kind(1.0d0)`, and no external BLAS, LAPACK, probability, or optimization
library.

## Main interface

```fortran
use mvtnorm
```

### Densities and simulation

```fortran
log_density = dmvnorm_one(x, mean, covariance, .true.)
log_t_density = dmvt_one(x, delta, scale, df, .true.)

normal_draws = rmvnorm(n, mean, covariance, seed)
t_draws = rmvt_shifted(n, scale, df, delta, seed)
kshirsagar_draws = rmvt_kshirsagar(n, scale, df, delta, seed)
```

Matrix-valued `dmvnorm` and `dmvt` evaluate one observation per row.

### Rectangular probabilities

```fortran
type(probability_control) :: control
type(probability_result) :: result

control = genz_bretz(maxpts=250000, abseps=1.0e-5_dp, seed=12345)
result = pmvnorm(lower, upper, mean, covariance, control)
result = pmvt(lower, upper, delta, scale, df, control)
```

`probability_result` contains the probability, a batch-based error estimate,
termination code, evaluation count, and message.

The available control constructors are:

```fortran
control = genz_bretz(maxpts, abseps, releps, seed, batches)
control = tvpack(abseps)
control = miwa(steps, abseps)
```

The modern high-dimensional engine uses randomized shifted Halton points,
antithetic pairing, conditional normal integration, variable reordering, and
batch error estimation. Student-t probabilities add inverse-chi-square
mixture integration. The shifted and Kshirsagar noncentral-t conventions are
both available through `pmvt` and `pmvt_kshirsagar`.

The deterministic two-dimensional normal route evaluates the Plackett
correlation integral adaptively. The `TVPACK` and `Miwa` controls preserve the
upstream method-selection API, but their higher-dimensional implementations
use the common modern conditional-integration engine rather than literal
line-by-line translations of the historical TVPACK and Miwa source. See
`PORTING_NOTES.md`.

### Simultaneous quantiles

```fortran
q = qmvnorm(p, mean, covariance, "lower", control)
q = qmvt(p, delta, scale, df, "both", control)
```

Lower, upper, and two-sided simultaneous quantiles are supported. A safeguarded
bracketed root solver reuses deterministic integration seeds across evaluations.

### Cholesky and precision parameterizations

The array-based API covers the numerical operations behind upstream
`ltMatrices` and `syMatrices`:

```fortran
chol = cov2chol(covariance, ok, message)
invchol = chol2invchol(chol, ok)
covariance = invchol2cov(invchol)
precision = invchol2pre(invchol)
correlation = chol2cor(chol)
partial_correlation = invchol2pc(invchol)

packed = pack_lower(chol, .true., .false.)
chol = unpack_lower(packed, dimension, .true., .false.)
```

Also included are lower-triangular multiplication and solution, cross-products,
log determinants, row/column scaling, standardization, covariance permutation,
`vectrick`, score depermutation, and score destandardization.

Fortran arrays replace R's classes and dimension-name metadata. A batch of
triangular matrices is represented by a rank-three array or processed in a
caller loop.

### Marginal and conditional Gaussian laws

```fortran
marginal = marginal_mvnormal(mean, covariance, indices)
conditional = conditional_mvnormal(mean, covariance, given_indices, given)
```

The conditional implementation uses the standard Schur-complement formulas and
supports arbitrary conditioning-index order.

### Exact and interval-censored likelihoods

```fortran
fit = ldmvnorm(observations, mean, covariance)
fit = lpmvnorm(lower, upper, mean, covariance, control)
fit = ldpmvnorm(exact, lower, upper, mean, covariance, control)
```

Observations are rows. `ldpmvnorm` treats the leading columns as exact and the
remaining columns as interval-censored, matching the upstream block convention.

Score interfaces are also provided:

```fortran
score = sldmvnorm(observations, mean, covariance)
score = slpmvnorm(lower, upper, mean, covariance, control)
score = sldpmvnorm(exact, lower, upper, mean, covariance, control)
```

Scores are returned for the mean followed by column-major packed lower-Cholesky
parameters. Exact-density scores and probability scores use deterministic
central differences. This is slower than the specialized upstream C score
engine but offers one consistent typed interface.

### Random-effects probability

`lpRR` and `slpRR` implement the discrete random-effects mixture probability
and its numerical score.

## Numerical conventions

- Covariance matrices must have positive diagonal entries.
- Densities require a positive-definite covariance matrix.
- Probability integration currently requires a positive-definite correlation
  matrix; semidefinite singular cases return `inform = 3` rather than silently
  adding a material ridge.
- Infinite limits use the Fortran `huge()` convention.
- Random simulation and randomized integration can be made deterministic by
  supplying seeds.
- Probability error estimates measure variation across randomized batches;
  they are statistical estimates, not rigorous deterministic bounds.

## License and provenance

The upstream package declares GPL version 2. This derivative project preserves
`GPL-2.0-only`, includes the complete license text, retains the unmodified
upstream source, and records SHA-256 manifests for original and translated
files.

The principal upstream reference is:

Alan Genz and Frank Bretz, *Computation of Multivariate Normal and t
Probabilities*, Lecture Notes in Statistics, Springer-Verlag, 2009.
