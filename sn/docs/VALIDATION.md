# Validation

The test suite is deterministic and contains six programs.

1. `test_univariate`
   - SN/normal and ST/Student-t reductions
   - analytic SN CDF at zero
   - SN and SC CDF/quantile inversion
   - direct/centered SN parameter round trip
2. `test_rng_moments`
   - seeded SN sample mean and variance against analytic cumulants
   - finite ST random generation
3. `test_multivariate`
   - one-dimensional multivariate APIs against univariate APIs
   - seeded multivariate SN mean
4. `test_sun`
   - SN-to-SUN density and CDF equivalence
   - marginal and affine transformations
   - seeded SUN generation
5. `test_fit`
   - SN regression slope and scale recovery
   - fitted/residual identities and prediction
6. `test_utilities`
   - half-vectorization and duplication dimensions
   - block-diagonal trace
   - normal and Student-t product sign probabilities
   - product quantile inversion
   - grouped normal likelihood fit

Run `make check` for a bounds-checked build and `make optimized` for the
optimized suite. `scripts/test.sh` runs the manifest check, both suites, and the
demonstration.
