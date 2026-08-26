# Testing

## Test programs

`test/test_all.f90`
: Core special-function checks, common distribution inversions, exact normal,
  exponential, and uniform results, Student-t and F log-tail references, flag
  behavior, and elemental array calls.

`test/test_reference.f90`
: 292 independently generated finite scalar references for unaffected original
  density, CDF, and quantile formulas. Expected values were calculated with
  independent SciPy distribution and special-function implementations.

`test/test_inversion.f90`
: CDF/quantile inversion for all 110 families that provide both procedures,
  using nondefault valid parameters.

`test/test_corrections.f90`
: Targeted checks for every corrected source family, including inversion,
  lower/upper complements, log probabilities, and density/CDF derivative
  identities.

`test/test_es_reference.f90`
: 25 independently integrated lower-tail expected-shortfall references,
  including common distributions and corrected package-specific families.

## Strict build

The strict validation build uses GNU Fortran 14.2 with:

```text
-std=f2018
-Wall -Wextra -Werror -pedantic
-fcheck=all -fbacktrace
-ffpe-trap=invalid,zero,overflow
-O0 -g
```

Modules are compiled in this order:

```text
vares_kinds.f90
vares_special.f90
vares_quadrature.f90
vares_distributions_01.f90 through vares_distributions_10.f90
vares.f90
```

FPM first builds the `rfortran-core` source dependency used by the central
Student-t and F compatibility wrappers.

The optimized validation build uses `-O3` with the same language and warning
checks.

The execution environment used for this port did not contain an `fpm`
executable. The manifest was parsed independently, and all FPM application,
example, and test targets were compiled directly with `gfortran`.

## Results

Both the strict and optimized configurations passed all five test programs,
the demonstration application, and the vector example with GNU Fortran 14.2.
The final source archive was also extracted into a clean directory and rebuilt
there to verify that no unarchived build products were required.
