# Reference generation

Independent references were generated in Python/NumPy/SciPy.

- ER bootstrap references independently implement the project's documented xorshift stream, resampling, studentization, centering, and one-/two-sided comparisons.
- CC references use direct matrix formulas and SciPy normal/chi-square survival functions.
- ESR coefficient references minimize the same Fissler-Ziegel objective with SciPy Nelder-Mead, starting from the Fortran solution and tightening tolerances.
- Core loss references are evaluated directly from the published algebra.

The reference constants are embedded in the test programs; no Python runtime is required to build or test the package.
