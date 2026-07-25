# Validation

## Environment

```text
GNU Fortran 14.2.0
LAPACK and BLAS from the system linker
```

## Debug configuration

```text
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface
-Werror -fcheck=all -fbacktrace
```

Command:

```text
make debug
```

Result:

```text
Term-structure tests passed.
GPL-2.0-or-later source license checks passed.
debug build, tests, and applications passed.
```

## Optimized configuration

```text
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface
-Werror -fbacktrace
```

Command:

```text
make release
```

Result:

```text
Term-structure tests passed.
GPL-2.0-or-later source license checks passed.
release build, tests, and applications passed.
```

## Numerical coverage

The regression suite checks:

- Nelson-Siegel and Svensson curve values at maturity zero
- Exact Nelson-Siegel recovery on noiseless synthetic data
- Svensson fitting with both SSE and source-compatible L1 objectives
- Positive fitted decay constants
- The 48-observation yield example published in the original package manual
- Nelson-Siegel improvement over a constant-rate fit
- Svensson improvement over Nelson-Siegel on the published example
- CSV parsing and all command-line model/objective modes
- Demo and example execution

The published example produced approximately:

```text
Nelson-Siegel RMSE: 1.827167e-03
Svensson SSE RMSE:  7.910187e-04
```

These figures document the tested Fortran build; they are not asserted to match
R's `nlminb` endpoint exactly.

## Scope statement

All non-plotting computational routines present in the attached `fBonds` source
are implemented and tested. There are no additional pricing, duration,
convexity, or cash-flow algorithms in the source package.
