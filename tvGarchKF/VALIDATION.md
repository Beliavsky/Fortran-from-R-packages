# Validation

The package contains five independent checked/optimized test programs.

1. `test_functions`
   - Polynomial, nonlinear, and trigonometric paths
   - The documented `3*(1-log(u))` argument
   - Invalid specification handling
2. `test_filter`
   - Exact hand-computed state variance, MSE, gain, state, variance, and
     objective fixtures
   - Missing observations and appended forecasts
   - Corrected constraint rejection
3. `test_simulation`
   - Exact deterministic recursion from supplied innovations
   - Seed reproducibility
4. `test_fit`
   - Simulation/fitting workflow
   - Criterion improvement and stationary positive solution
   - Five-decimal compatibility output
5. `test_tv_parameter`
   - Global and overlapping-window GARCH fitting
   - Midpoint construction and valid local estimates

Validated with GNU Fortran 14.2.0 using:

- checked: `-std=f2018 -Wall -Wextra -Werror -pedantic -O0 -g -fcheck=all -fbacktrace`
- optimized: `-std=f2018 -Wall -Wextra -Werror -pedantic -O3`

Both configurations passed all tests and the demonstration program.
