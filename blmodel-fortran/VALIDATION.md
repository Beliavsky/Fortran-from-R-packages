# Validation

## Environment

- GNU Fortran 14.2.0
- Fortran 2018 mode
- Linux x86-64

FPM was not installed in the validation runtime. The manifest was parsed as TOML
and the same source, test, application, and example graph was compiled directly.

## Checked compiler flags

```text
-std=f2018
-Wall -Wextra -Wpedantic
-Wconversion-extra -Wimplicit-interface
-Werror
-fcheck=all -fbacktrace
-O0 -g
```

## Test programs

- `test_distributions`: fixed independent normal, Student-t, and
  power-exponential density references; invalid covariance rejection
- `test_equilibrium`: annualized moments, diagonal helpers, MAD/CVaR inverse
  optimization, market-return identity, and elliptical equilibrium
- `test_posterior`: return recentering, diagonal view covariance, fixed posterior
  probabilities, and normalization
- `test_blmodel`: complete high-level elliptic workflow and invalid-prior handling

## Results

```text
test_distributions: PASS
test_equilibrium: PASS
test_posterior: PASS
test_blmodel: PASS
```

The demo and both examples compile and execute in checked and optimized builds.
The exact release archive is extracted into a clean directory and validated
again before publication.
