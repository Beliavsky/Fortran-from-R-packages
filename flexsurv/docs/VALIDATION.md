# Validation

The direct validation scripts compile every source module and test with GNU
Fortran using:

```text
-std=f2008
-Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

The v0.3.0 suite contains 23 test programs. Coverage includes distribution
identities and quantiles, censoring/truncation likelihoods, source-style starts
and optimizer alternatives, custom models, fractional polynomials, basic and
full mixture models, analytic-versus-numerical Louis information, mixture
multi-state bootstrap summaries, both Royston-Parmar and splines2ns spline
bases, spline fitting/interactions, QP initialization, deSolve-driven Markov
prediction, semi-Markov simulation with predictable time-dependent covariates,
shared-regression cross-transition covariance, bootstrap uncertainty,
final-state summaries, AJ comparisons, advanced rate-table standardization, and
fixed-parameter right-truncation inference.

The example is run after every test program. Project Fortran source/test/example
files are additionally checked to contain no tabs and no lines longer than 132
columns.
