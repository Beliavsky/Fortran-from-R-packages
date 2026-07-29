# Testing

The release script builds the complete source tree twice:

1. strict debug configuration with bounds checking, floating-point traps, backtraces, warnings, and implicit-interface checks;
2. optimized `-O3` configuration.

It runs:

- `test_core_algorithms`
  - normal and Student-t probability checks
  - CDF/quantile inversion
  - hidden `arfilt` coefficient reference
  - Kupiec and Christoffersen likelihood-ratio references
  - all four loss-function references
- `test_traffic`
  - violation counts
  - binomial traffic-light probabilities
  - ES breach severity
  - weighted absolute deviation
- `test_varcast_models`
  - end-to-end parametric execution of all six models
  - positive finite volatility
  - finite VaR and ES
  - ES/VaR ordering
- `test_semiparametric`
  - short-memory local-polynomial scale estimation
  - long-memory scale estimation
  - Student-t fitting
  - retained smoothing results

The script also builds and runs `ufrisk_demo` and `semiparametric_figarch`.

The known GNU linker note about an executable stack is caused by the supplied rugarch module's use of an internal procedure as an optimizer callback. It is not a test failure and does not affect numerical output.
