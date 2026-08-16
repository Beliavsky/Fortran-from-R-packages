# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `orthopolynom` 1.0-6.1.
- Ported all exported orthogonal-polynomial families, recurrence tables,
  norms, weights, and normalized polynomial constructors.
- Ported monic recurrences, Jacobi matrices, polynomial-list calculus,
  values, roots, change-of-basis powers, Pochhammer helpers, and scaling.
- Used the supplied `polynom-fortran-v0.1.0` as a local FPM dependency.
- Added strict runtime and independent numerical regression tests.
