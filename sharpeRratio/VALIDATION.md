# Validation

The validation suite contains five programs:

1. `test_records`: exact record counts, constant-sign identities, seeded
   permutation reproducibility, and confidence ordering.
2. `test_calibration`: exact spline-knot values, branch behavior, and sign
   symmetry of the SNR map.
3. `test_statistics`: R type-7 quantiles, sample variance, Gaussian
   classification, and heavy-tail classification.
4. `test_estimator_fixed_nu`: fixed-tail positive/negative estimates,
   deterministic intervals, and filtering of non-finite observations.
5. `test_estimator_automatic`: Gaussian sentinel selection and end-to-end
   Student-t simulation, automatic tail fitting, permutation estimation, and
   SNR calculation.

The demonstration simulates skewed-Student-compatible heavy-tailed returns and
prints the SNR, confidence range, fitted tail exponent, mean record imbalance,
and observation count.

Run `scripts/validate.sh` or `scripts/validate.bat` to execute checked and
optimized builds.

Validated in this translation environment with:

- GNU Fortran (Debian 14.2.0-19) 14.2.0
- checked flags: `-std=f2018 -O0 -Wall -Wextra -Wconversion-extra -Werror -fcheck=all -fbacktrace`
- optimized flags: `-std=f2018 -O2 -Wall -Wextra -Wconversion-extra -Werror`

FPM was not installed in the validation environment. The manifest was parsed
with Python's TOML parser, and the same source/module graph was compiled through
the included GNU Fortran scripts.
