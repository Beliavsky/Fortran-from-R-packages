# Validation

## Environment

- Compiler: GNU Fortran (Debian 14.2.0-19) 14.2.0
- Platform: Linux x86-64
- Fortran standard: Fortran 2018
- R available: no
- fpm available: no

## Strict debug configuration

`make debug` invokes:

```
-std=f2018 -O0 -g -Wall -Wextra -Wimplicit-interface -Werror
-fcheck=all -fbacktrace
```

It compiles all library modules, two numerical test programs, the demonstration,
and the constrained example. It then runs all executables and the license audit.

## Optimized configuration

`make release` invokes:

```
-std=f2018 -O2 -Wall -Wextra -Wimplicit-interface -Werror -fbacktrace
```

The same tests and applications are rebuilt and executed.

## Executed numerical tests

### jDE and synchronous population jDE

- deterministic midpoint/bounce-back bound correction
- scalar and vector equality-tolerance conversion
- median calculation
- asynchronous jDE on the two-dimensional sphere
- dither-only jDE with maximum-based stopping reference
- equality-constrained quadratic optimization
- inequality-constrained quadratic optimization
- synchronous SPJDE on the Aluffi-Pentini problem
- synchronous SPJDE equality and inequality constraints
- injected initial optimum with a three-vector random population
- convergence codes and saved final-population output

### NCDE

- all four Becker-Lago global minima and quadrant coverage
- both minima of an inequality-constrained one-dimensional problem
- both minima of an equality-constrained one-dimensional problem
- strict archive feasibility
- fixed and automatically identified niche radii
- enabled and disabled archive-neighbor reinitialization
- sorted final population and archive outputs

## Observed output

```
JDE and SPJDE tests passed.
NCDE tests passed.
JDE sphere:  -4.374198E-06   4.108426E-06 f=   3.601277E-11
SPJDE Aluffi:  -1.046682E+00  -3.557022E-06 f=  -3.523861E-01
NCDE archived solutions: 4
parameters:   4.99994546E-01   4.99995454E-01
objective:   4.99990000E-01
equality residual:  -9.99999867E-06
GPL-2.0-or-later source license checks passed.
debug build, tests, and applications passed.
release build, tests, and applications passed.
```

Because the algorithms are stochastic and use Fortran's random-number stream,
exact R iteration histories are not claimed. Tests use fixed Fortran seeds and
verify known optima, constraints, archives, and result structure.
