# Validation status

## Completed

The project was compiled with gfortran 14.2.0 using Fortran 2018, warnings, explicit-interface checks, bounds/runtime checks, and backtraces.

The smoke suite covers:

- Quadratic map dimensions and initial value
- Jarque-Bera test sanity
- Runs-test probability range
- AR fitting and positive CSS
- GARCH positivity and persistence constraint
- BDS execution
- Stationary bootstrap sample range

The extended suite covers:

- Fourier and amplitude-adjusted surrogates
- ADF, PP, and KPSS tests
- Terasvirta and White tests
- Phillips-Ouliaris test
- Equality-constrained portfolio optimization

Both suites pass in the supplied environment.

## Not completed

- Broad numerical comparison against R `tseries` 0.10-62
- Exact replication of R optimizer trajectories and warnings
- Large-sample performance testing
- Testing with ifx, nvfortran, Flang, or multiple operating systems
- Exhaustive invalid-input and singular-system tests
- Exact random-number equivalence with R
- Validation of every combination of optional arguments

The repository should therefore be treated as experimental.
