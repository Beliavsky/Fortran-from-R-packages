# Independent reference generation

Fixed reference values in the tests were generated with Python, NumPy, and
SciPy by implementing the likelihood equations directly and minimizing them
with `scipy.optimize.minimize(method='Nelder-Mead')`.

The reference implementation did not call the Fortran library. In particular,
it independently evaluated:

- the GEV likelihood on a deterministic GEV quantile grid;
- the GPD likelihood on a deterministic GPD quantile grid;
- fixed GEV, GPD, and point-process likelihood vectors;
- the logistic bivariate exponent measure and marginal tail formulas.

The expected values are embedded in `test/test_fitting.f90` and
`test/test_bivariate.f90` to make the tests reproducible without Python.
