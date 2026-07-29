# Testing

## Toolchain used

- GNU Fortran 14.2.0
- system BLAS and LAPACK
- Fortran 2018 mode

The environment used to prepare the port did not contain the `fpm` executable.
The FPM manifest was parsed independently, and every source, test, application,
and example target was compiled directly with `gfortran`.

## Strict debug flags

```text
-std=f2018
-Wall -Wextra -Werror -pedantic
-fimplicit-none
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

## Test coverage

`test/test_strategies.f90` checks:

- fixed CPPI positions, portfolio prices, and churn;
- DPPI moving targets on a falling market path;
- Black-76 OBPI target premium, deltas, trades, and portfolio prices for buyers and sellers;
- SLPI trigger behavior for buyers and sellers;
- SHPI scheduled hedge behavior; and
- invalid volume and expiry inputs.

The fixed values were generated independently from the equations rather than
copied from the Fortran output.

`test/test_msfc.f90` checks:

- the ten-contract example documented by the R package;
- polynomial and daily-curve dimensions;
- continuous repricing of every input contract;
- five independent fixed daily curve values;
- nonconstant-prior handling and exact repricing after prior adjustment; and
- contract inclusion filtering.

An independent NumPy assembly and solve of the KKT system matched the first five
MSFC values to the displayed double-precision digits. The Fortran curve repriced
the ten example contracts with maximum absolute error below `7e-12` in the
strict build.

## Commands

With FPM:

```text
fpm test
fpm run
fpm run --example msfc_example
fpm run --example strategies_example
```

Manual linking requires `-llapack -lblas`.
