# Porting notes

## Direct translations

- Gaussian complete-data AR(1) maximum likelihood formulas.
- Gaussian missing-data EM sufficient statistics.
- Conditional Gaussian means and covariance blocks for missing intervals.
- Student-t AR(1) latent-scale IRLS updates.
- Student-t AR(1) Gibbs sampling of scales and missing blocks.
- OHLC and volume log-scale wrappers.
- Outlier detection based on one-sided Gaussian or Student-t tail probability.

## Adaptations

### R objects

R vectors, matrices, attributes, lists, `zoo`, and `xts` objects are replaced by
Fortran arrays and derived result types. Missing/outlier locations are stored as
one-based integer arrays in the result.

### Randomness

A deterministic xorshift/Box-Muller/Marsaglia-Tsang generator is included.
Results therefore do not reproduce R's RNG stream, but a fixed seed is portable
across supported Fortran compilers.

### Student-t VAR missing-data estimator

The R package uses stochastic EM with multiple parallel Gibbs chains and can
partition missingness patterns. This translation preserves multivariate-t
weights, weighted VAR updates, scatter estimation, degrees-of-freedom fitting,
and the omit-missing mode. Its non-omit path uses sequential conditional means
for missing contemporaneous components rather than parallel MCMC. This is a
numerical adaptation and can differ from the R estimates on heavily incomplete
data.

### Omitted infrastructure

`plot_imputed`, ggplot/graphics integration, R S3/attribute behavior, parallel
socket clusters, package datasets in RData format, vignettes, and website
infrastructure are not compiled. Original source files are retained in
`original/`.
