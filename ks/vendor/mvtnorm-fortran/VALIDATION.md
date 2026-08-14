# Validation report

## Compiler

GNU Fortran 14.2.0 was used for direct validation.

Checked build flags:

```text
-std=f2018 -pedantic -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -fcheck=all -fbacktrace -O0
```

Optimized build flags:

```text
-std=f2018 -pedantic -Wall -Wextra -Wconversion-extra
-Wimplicit-interface -Werror -O2
```

## Test suites

```text
test_conditioning_likelihood: PASS
test_distributions: PASS
test_probabilities: PASS
test_quantiles_scores: PASS
test_triangular: PASS
```

The demo and both examples compile and run in the optimized build.

## Independent references

The suites include independent SciPy and NumPy references for:

- Multivariate normal and Student-t log densities.
- Bivariate normal probabilities.
- A three-dimensional normal rectangle probability.
- A three-dimensional Student-t rectangle probability.
- Conditional Gaussian means and covariance matrices.
- Exact, interval-censored, and mixed-data log likelihoods.
- Distribution CDF/quantile inversion.
- Cholesky, covariance, precision, correlation, and partial-correlation
  identities.
- Independent-component simultaneous quantiles.
- Random-effects mixture probabilities.

Simulation tests check sample means and covariance with deterministic seeds and
Monte Carlo tolerances.

## Probability reference cases

For the three-dimensional reference model,

```text
mean = [0.2, -0.1, 0.3]

covariance =
[ 1.0   0.4  -0.2  ]
[ 0.4   1.5   0.25 ]
[-0.2   0.25  0.8  ]

lower = [-1.0, -0.5, -0.8]
upper = [ 0.7,  1.2,  0.6]
```

validated probabilities are approximately:

```text
multivariate normal: 0.138861665
multivariate t, df=7: 0.13325511
```

The checked Fortran run used 400,000 antithetic evaluations and passed the
specified integration tolerances.

## Scripts

```text
scripts/validate.sh
scripts/validate_optimized.sh
scripts/validate.bat
```

FPM was not installed in the validation environment. The manifest was parsed as
TOML and follows FPM automatic source, application, example, and test discovery.
