# Validation

The translation is validated with four independent test programs:

- `test_moments`: sample skewness/kurtosis and both model recursions
- `test_likelihood_forecast`: fixed likelihood references and forecast recursion
- `test_constraints`: constraint vectors, feasibility, and invalid-model penalties
- `test_estimation`: deterministic estimation, feasibility, and objective improvement

The fixed recursion and likelihood values were calculated independently with
NumPy formulas matching the documented model equations.

## Compiler settings

Checked build:

```text
-std=f2018 -O0 -g -Wall -Wextra -Werror -Wconversion-extra
-Wimplicit-interface -fcheck=all -fbacktrace
```

Optimized build:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wconversion-extra
-Wimplicit-interface
```

Both builds compile all library modules, tests, the demo, and both examples.

## Expected test output

```text
test_constraints: PASS
test_estimation: PASS
test_likelihood_forecast: PASS
test_moments: PASS
```

The final release archive is extracted into a clean directory and the complete
checked validation script is rerun before packaging.

## Final release audit

The release contains 1,090 lines of Fortran across library, test, demo, and
example units. The FPM manifest parses as TOML. All translated text is ASCII,
all free-form Fortran lines are at most 132 columns, and every Fortran unit has
`implicit none` and a GPL SPDX identifier.

The exact ZIP was extracted into an empty directory, its original and translated
checksum manifests were verified, and `scripts/validate.sh` completed with:

```text
clean_archive_validation: PASS
```

FPM itself was unavailable in the validation environment, so the numerical
build used the included direct GNU Fortran scripts.
