# Validation

The automated test program covers:

- all 12 Gauss-Kronrod rules on a degree-six polynomial;
- the Beta(3,10) vignette probability on `[0.3, 0.8]`;
- reversed finite intervals;
- semi-infinite and doubly infinite one-dimensional integrals;
- Cuhre default rules in two, three, and four dimensions;
- a correlated bivariate-normal rectangle probability;
- a mixed semi-infinite/doubly-infinite integral;
- invalid multidimensional bounds;
- unconstrained Rosenbrock optimization;
- bounded Rosenbrock optimization;
- an active-bound solution;
- polymorphic optimizer user data; and
- logistic-regression coefficients and log likelihood.

Reference values were independently calculated with analytic formulas and
SciPy. In particular:

- `P(-1<X1<1, -1<X2<1)` for correlation 0.5 is
  `0.49797177783920804`;
- the logistic-regression coefficients are approximately
  `(-0.27184836, 1.08739343)`; and
- its maximized log likelihood is `-4.941579983434302`.

Both full runtime-check and optimized GNU Fortran builds are exercised by the
provided scripts.
