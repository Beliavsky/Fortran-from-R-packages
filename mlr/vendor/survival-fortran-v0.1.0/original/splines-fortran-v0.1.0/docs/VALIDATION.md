# Validation

The test program covers:

- Cubic B-spline values at interior points and both endpoints.
- First and second derivatives against independent numerical references.
- Partition of unity.
- R type-7 quantiles used for automatic knots.
- Cubic regression-basis extrapolation values.
- Natural interpolation of a six-point data set, including an independently
  computed reference value and zero second derivatives at both endpoints.
- Agreement between B-spline and local polynomial representations.
- Periodic interpolation, data reproduction, wrapping, and sine accuracy.
- Natural regression-basis dimensions and linear extrapolation.
- Monotone inverse-spline construction.
- Linear interpolation with unsorted input pairs.

The release scripts compile with Fortran 2018, warnings enabled, and full GNU
Fortran runtime checking. A separate optimized build is also validated.

A separate irregular-knot comparison was run against SciPy's cubic
`BSpline` evaluator at six points for derivative orders 0 through 3. The
largest absolute difference was approximately `3.6e-15`.
