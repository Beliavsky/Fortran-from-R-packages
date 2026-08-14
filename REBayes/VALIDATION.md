# Validation

Validation environment: GNU Fortran 14.2.0.

All permanent tests are compiled with free-source Fortran 2018, warnings as
errors, explicit-interface warnings, and runtime bounds/allocation checks.
The optimized release pass uses:

    -std=f2018 -O2 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all

Permanent tests cover:

- an analytically symmetric KW problem;
- 100 deterministic randomized two-component KW problems versus the exact
  one-dimensional derivative root;
- Gaussian and binomial mixture fits;
- KW quantiles/smoothing and Huber utilities;
- general-D regularized logistic regression;
- MEDDE nonnegativity and quantiles;
- weighted/repeated-measures variance and joint mean/variance mixtures;
- a fixed independent noncentral-t density reference.

Additional development validation:

1. 60 randomized finite-grid mixture problems with 3--9 components were compared
   against SciPy SLSQP.  All 60 succeeded.  Maximum absolute objective difference
   was approximately 8.45e-11; maximum fitted-density difference was about
   2.90e-7; maximum Fortran KKT gap was about 1.0e-10.

2. General-D RLR was compared on 20 randomized problems against an independent
   split-variable constrained SLSQP formulation.  Maximum objective discrepancy
   was approximately 1.16e-9.

3. The noncentral-t density was compared on 120 randomized cases against
   `scipy.stats.nct.pdf`.  With the cached 96-point quadrature rule, maximum
   absolute error was approximately 2.20e-7 and maximum relative error among
   densities above 1e-8 was approximately 1.52e-6.
