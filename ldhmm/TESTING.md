# Testing

The regression suite covers:

- Lambda-distribution normalization, CDF symmetry, and analytic moments
- Transition initialization and stationary probabilities
- Natural/working parameter round trips
- Forward/backward likelihood identities
- Posterior, forecast-state, conditional-density, and forecast-density normalization
- Viterbi state ranges and pseudo-residual finiteness
- Decoded statistic histories
- Moving averages, outlier removal, ACF utilities, and simulation
- Likelihood improvement under BFGS and Nelder-Mead fitting

Run with FPM:

```sh
fpm test
```

Or with GNU Fortran:

```sh
./scripts/test_gfortran.sh
```

The direct script compiles with Fortran 2018 mode, warnings, and runtime checks.
An optimized build script is also supplied.
