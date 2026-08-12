# Changelog

## 0.1.0

- Initial modern Fortran/FPM translation of `linprog` 0.9-6.
- Port original two-phase tableau simplex and sensitivity calculations.
- Add optional explicit dual solve.
- Integrate user-supplied `lpSolve-fortran` 0.1.0 as an FPM path dependency.
- Port restricted MPS reader/writer.
- Preserve documented equality-constraint bug in the internal solver and test
  the correct lpSolve-backed behavior separately.
- Add five regression tests and two examples.
