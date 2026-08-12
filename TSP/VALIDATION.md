# Validation

The source was compiled and tested using GNU Fortran 14.2.0.

Strict compilation used:

```text
gfortran -std=f2018 -Wall -Wextra -Werror -Wimplicit-interface \
  -fcheck=all ...
```

The checked test suite covers:

- tour lengths and Euclidean coordinate distances
- positive/negative infinity replacement
- insertion-cost calculation
- all native construction heuristics
- nearest-neighbor and repetitive nearest-neighbor
- simulated annealing returning valid tours
- asymmetric 2-opt behavior
- dummy-city insertion
- ATSP-to-TSP reformulation and filtering
- single and multiple tour cuts
- TSPLIB ATT and GEO reference examples from the upstream tests
- TSPLIB explicit symmetric round trip
- 100 deterministic randomized asymmetric 2-opt cases, checking both
  non-worsening behavior and local 2-optimality by direct enumeration of every
  allowed reversal after convergence
- 100 deterministic randomized insertion cases, checking each insertion delta
  against direct cycle-length recomputation and checking all insertion
  heuristics return valid permutations

Additional development stress runs used 500 randomized asymmetric 2-opt cases
and 500 randomized insertion cases.

The example program also compiles with the same strict flags and runs
successfully.

The sandbox did not provide an `fpm` executable, so the FPM manifest could not
be invoked directly here. `fpm.toml` is included and the exact source/test/example
layout was compiled directly with gfortran.
