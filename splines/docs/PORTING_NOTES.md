# Porting notes

## Representation

The R package uses lists, attributes, and S3 classes. The Fortran port uses two
explicit derived types:

- `b_spline_t`: knot sequence, B-spline coefficients, order, and natural or
  periodic behavior.
- `poly_spline_t`: breakpoints and local polynomial coefficients in ascending
  powers of `x - knot(i)`.

Both types provide scalar and vector `evaluate` bindings.

## B-spline kernel

The local basis and derivative routines are direct modern Fortran translations
of the algorithms in `src/splines.c`. Difference tables and work arrays are
local allocatables, making evaluation reentrant and thread-safe; the original C
file used static work pointers.

## Linear algebra

The original R code delegated square solves and QR operations to R. This port
uses a dependency-free partial-pivoting linear solver. Natural regression bases
are formed by an orthonormal completion of the two endpoint second-derivative
constraint directions. The resulting columns span the same natural-spline
space, but column signs and rotations need not match R's QR basis exactly.

## Behavior differences

- `spline_design` sorts a copy of the knot vector, matching the R routine.
- `linear_interp` accepts query points in arbitrary order. The old R wrapper
  computed a sorted query vector but passed the unsorted vector to its C loop;
  the Fortran routine implements the intended order-independent behavior.
- Invalid inputs are reported through optional integer `status` arguments.
  Functions that cannot return a value use IEEE quiet NaN where appropriate.
- R missing-value and name/attribute propagation are not reproduced.

## Bounds

Natural interpolation objects extrapolate linearly. Periodic objects wrap
arguments by the stored period. Ordinary B-spline objects return NaN and a
nonzero status outside their valid knot interval. The `bs_basis` routine follows
R's Taylor-polynomial extrapolation behavior outside boundary knots.
