# Validation

The test suite contains five programs:

1. `test_continuous`: CDF/quantile round trips and density fixtures across all
   continuous families.
2. `test_discrete`: normalization, CDF identities, quantile conventions, and
   tail/log behavior across all discrete families.
3. `test_multivariate`: bivariate, mixture, Dirichlet, multinomial,
   Dirichlet-multinomial, and multivariate-hypergeometric identities.
4. `test_random`: deterministic seeding, support checks, and selected empirical
   moments.
5. `test_api`: log probabilities, upper tails, invalid inputs, zero-length
   draws, and source-compatible/corrected edge modes.

Both checked and optimized builds are required to pass. The checked build uses
GNU Fortran runtime checks; both builds use strict Fortran 2018 mode and treat
warnings as errors.

The upstream R test suite was used as a coverage guide. R was not available in
the build environment, so bit-for-bit comparisons against a live R installation
were not performed. Formula identities, normalization, inverse consistency, and
simulation moments provide independent numerical validation.
