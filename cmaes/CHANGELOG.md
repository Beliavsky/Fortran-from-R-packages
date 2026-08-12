# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of cmaes 1.0-12.
- Ported the complete CMA-ES computational loop and control defaults.
- Added box constraints, multiplicative penalty, minimization/maximization,
  stop rules, flat-fitness escape, and best-solution retention semantics.
- Added scalar and vectorized objective interfaces.
- Added sigma/eigenvalue/population/value diagnostics and population extraction.
- Added sphere, random, shifted Rosenbrock, Rastrigin, shift, rotation, and bias
  helpers.
- Replaced R RNG with a deterministic local generator and R symmetric eigen
  decomposition with LAPACK DSYEV.
- Added five regression test programs and two examples.
- Preserved GPL-2 licensing and the complete original attached source tree.
