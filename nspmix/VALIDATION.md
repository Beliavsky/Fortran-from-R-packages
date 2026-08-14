# Validation

Validated with GNU Fortran 14.2.0 using:

```text
-std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

Permanent tests:

1. `test_hcnm` - constrained two-component likelihood versus an independent dense scalar search.
2. `test_hcnm_random` - 80 randomized two-component problems versus an independent golden-section likelihood maximizer.
3. `test_families` - analytical normal, geometric and negative-binomial component-density references.
4. `test_cnm` - weighted Poisson-mixture NPMLE support/mass constraints.
5. `test_semiparametric` - CVPS semiparametric and grouped-logistic fixed-support paths.

All tests pass with runtime bounds/allocation checking and implicit external
interfaces disabled. The `poisson_mixture` example also compiles and runs under
the same flags.
