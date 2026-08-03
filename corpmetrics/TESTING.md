# Testing

The release is tested in two clean configurations with GNU Fortran 14.2:

```text
./scripts/test_gfortran.sh strict
./scripts/test_gfortran.sh release
```

The strict configuration uses:

```text
-std=f2018 -O0 -g -Wall -Wextra -Werror -pedantic
-fcheck=all -ffpe-trap=invalid,zero,overflow -fbacktrace
```

The release configuration uses `-O3` with warnings treated as errors.

## Test programs

- `test_metrics`: all balance-sheet ratios, CAPM, all three DDM branches,
  annual/semiannual fixed-income calculations, and income-statement metrics
- `test_investment`: term-specific NPV, two IRR references, public low-level
  functions, dimension errors, and an unbracketed root
- `test_loan`: amortization references, zero-rate limit, final balance, and
  ties-to-even two-decimal rounding

Every application and example is also compiled and run by the build script.
Reference values were calculated independently from the equations in the R
source using Python/NumPy and, for the nontrivial IRR case, SciPy Brent root
solving. See `REFERENCE_GENERATION.md`.
