# Validation report

## Environment

The translated release contains 3,747 lines across 22 Fortran source, test, demo, and example files.

The release was validated with GNU Fortran 14.2.0 in Fortran 2018 mode.

Checked compilation uses:

```text
-std=f2018 -O0 -g -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Wno-compare-reals -fcheck=all -fbacktrace -Werror
```

An optimized `-O2` build is also executed.

FPM was not available in the construction environment, so the included direct
validation scripts were used. The FPM manifest is parsed separately and follows
standard automatic source, application, example, and test discovery.

## Test programs

### `test_continuous`

Checks:

- density/CDF/quantile identities for major continuous families;
- numerical moments and limited moments;
- distribution boundaries;
- seeded random generation smoke tests;
- inverse-Gaussian and Gumbel calculations.

### `test_supplements`

Checks moments, limited moments, and MGFs for common base distributions.

### `test_discrete`

Checks:

- logarithmic PMF/CDF/quantile identities;
- zero-truncated normalization;
- zero-modified point masses and CDFs;
- binomial and negative-binomial boundaries;
- Poisson-inverse-Gaussian probabilities and seeded RNG behavior.

### `test_aggregate`

Checks:

- direct convolution;
- exact compound distributions;
- Poisson, binomial, and negative-binomial Panjer recursion;
- aggregate mean and variance;
- VaR and CTE;
- discretization and normal approximations.

### `test_phase_credibility`

Checks:

- phase-type density, CDF, moments, MGF, and random generation;
- matrix exponential behavior;
- Buhlmann-Straub credibility;
- Poisson-gamma, Bernoulli-beta, and normal-normal credibility.

### `test_grouped_risk`

Checks:

- grouped mean, variance, quantile, and ogive interpolation;
- deductible, limit, and coinsurance transformations;
- adjustment-coefficient root solving;
- explicit exponential-claim ruin probability.

## Expected checked output

```text
test_aggregate: PASS
test_continuous: PASS
test_discrete: PASS
test_grouped_risk: PASS
test_phase_credibility: PASS
test_supplements: PASS
```

The demo and both examples are also built and run in checked and optimized
configurations.

## Release audits

The release process verifies:

- valid TOML manifest;
- translated text is ASCII;
- free-form Fortran lines do not exceed 132 columns;
- every translated Fortran file has an SPDX identifier;
- every module/program uses `implicit none`;
- original and translated SHA-256 manifests;
- the exact final ZIP rebuilds and passes all tests from a clean extraction.

## Known numerical limitations

- Native special functions target practical actuarial ranges, not every extreme
  tail or ill-conditioned parameter combination.
- PIG probabilities use a direct Bessel-K integral and may be slower than a
  specialist library.
- Phase-type matrix routines are designed for modest dimensions.
- The selected phase-type ruin routine has narrower scope than upstream `ruin`.
- RNG streams are deterministic within this project but differ from R.
