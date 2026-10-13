# splines2: modern Fortran / FPM

A translation of the computational API of Wenjie Wang and Jun Yan's R package
**splines2 0.5.4**, licensed GPL-3.0-or-later. It provides B-, M-, I-, C-,
periodic B/M, natural cubic, cardinal natural, and generalized Bernstein bases.
No R, Rcpp, Armadillo, separately installed BLAS, or LAPACK is needed at runtime.

## Build and use

Place this directory at the repository root alongside `splines` and
`rfortran-core`, then run from this directory:

```sh
fpm build
fpm test
fpm run --example
```

The package uses sibling FPM path dependencies. They are intentionally not
included in this directory. GNU Fortran and FPM are sufficient, including on
Windows: no platform-specific shell invocation, POSIX API, or system-library
link is required. The verified build platform was Linux with GNU Fortran 14.2.0
and FPM 0.13.0; Windows was not available for execution.

```fortran
use splines2, only : dp, spline_basis, b_spline, derivative_spline
implicit none
type(spline_basis) :: basis, slope
basis = b_spline([0.0_dp, 0.5_dp, 1.0_dp], &
   knots=[0.25_dp, 0.75_dp], intercept=.true.)
if (basis%status /= 0) error stop basis%message
slope = derivative_spline(basis)
print *, basis%values
```

## Public API

All constructors return `type(spline_basis)`. Its `values` matrix has one row
per input coordinate. Check `status == SPLINES2_OK` before using it; invalid
constructor arguments return an empty matrix and a descriptive `message`.
The result retains `x`, `knots`, `boundary_knots`, `degree`, `intercept`,
`periodic`, `derivs`, `integral`, and `scale`. Natural bases expose `transform`;
C-splines expose normalization `scales`. `equidistant_fallback` reports when
quantile-based knots needed an equally spaced replacement.

- `b_spline(x, df, knots, degree, intercept, boundary_knots, periodic, derivs, integral)`
- `m_spline(x, df, knots, degree, intercept, boundary_knots, periodic, derivs, integral)`
- `i_spline(x, df, knots, degree, intercept, boundary_knots, derivs)`
- `c_spline(x, df, knots, degree, intercept, boundary_knots, derivs, scale)`
- `natural_spline(x, df, knots, intercept, boundary_knots, trim, derivs, integral)`
- `nsk(x, df, knots, intercept, boundary_knots, trim, derivs, integral)`
- `bernstein_poly(x, degree, intercept, boundary_knots, derivs, integral)`
- `dbs(x, derivs, df, knots, degree, intercept, boundary_knots, periodic)`
- `ibs(x, df, knots, degree, intercept, boundary_knots, periodic)`
- `predict_spline(object, newx)` reuses stored knots and boundaries
- `derivative_spline(object, derivs)` adds a derivative order, default one

Except `x` (or `object`/`newx`), constructor arguments are optional. Use keyword
arguments to avoid positional ambiguity. Degree defaults to 3, derivatives to
0, periodic/integral to false, C scaling to true, and natural trim to zero.
`dbs` defaults to one derivative. Intercept defaults to true for I/C and false
for the other families. Aliases `bsp`, `msp`, `isp`, `csp`, `nsp`, `bpoly`
are public Fortran generics with the same signatures as their target constructors.
`dp` is imported from shared `r_kinds` and re-exported.

## Numerical behavior and boundaries

- Nonempty explicit knots override `df`; otherwise type-7 sample quantiles
  determine the requested column count. Quantiles exclude NaNs and coordinates
  outside the boundaries; periodic quantiles use wrapped coordinates.
- Interior knots are sorted and must be strictly inside the boundaries. B/M/I/C
  support repeated interior knots through multiplicity `degree+1`; natural and
  periodic bases require distinct knots. Periodic B/M require at least
  `degree-1` interior knots.
- Default boundaries span nonmissing coordinates; natural `trim` replaces them
  with quantiles. Explicit boundaries override trim. Bernstein boundaries must
  be increasing; other families sort a supplied pair.
- Isolated NaNs produce entire NaN output rows. Empty/all-NaN coordinates,
  infinities, nonfinite knots, degenerate boundaries, invalid degree/derivative,
  insufficient df, and zero remaining basis columns return status errors.
- B/M/I/C use endpoint-polynomial continuation outside the boundaries. Natural
  splines continue linearly. Bernstein evaluation requires in-bound coordinates.
- Integrals follow upstream conventions. For ordinary B/M and unscaled C
  (the integral of I), endpoint polynomials continue algebraically. Natural
  integrals are zero left of the left boundary, and have quadratic continuation
  to the right. Periodic integrals are zero below the left boundary and
  accumulate whole-period areas to the right. Therefore differentiating an
  integral outside the lower boundary need not reproduce its separately
  requested basis value for natural or periodic families.
- Setting `integral=.true.` and `derivs=k` requests the kth derivative of the
  integrated basis: `k=1` returns the original basis. I/C degree names refer to
  their associated M-spline degree, giving actual degrees `degree+1`/`degree+2`.
- Natural-spline columns use the upstream closed-form, column-normalized
  transform, not the unrelated QR basis from `splines::ns`. `nsk` is cardinal at
  the boundary/interior knots when the first column is retained.

`bernstein_poly` and its arithmetic kernel are pure. Knot construction and
validation helpers are pure where possible. Other basis constructors call the
existing impure `splines` evaluator (and `nsk` its solver); they cannot honestly
be declared pure without modifying that independent dependency.

## Translation coverage

Package status: **substantial**. **15 of 15 (100%)** exported computational R
functions have meaningful mappings. This fraction measures mapped computational
functions, not complete R compatibility. It includes distinct exported aliases.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `bSpline` | `splines2` | `b_spline` | substantial |
| `bsp` | `splines2` | `bsp` | substantial |
| `dbs` | `splines2` | `dbs` | substantial |
| `ibs` | `splines2` | `ibs` | substantial |
| `mSpline` | `splines2` | `m_spline` | substantial |
| `msp` | `splines2` | `msp` | substantial |
| `iSpline` | `splines2` | `i_spline` | substantial |
| `isp` | `splines2` | `isp` | substantial |
| `cSpline` | `splines2` | `c_spline` | substantial |
| `csp` | `splines2` | `csp` | substantial |
| `naturalSpline` | `splines2` | `natural_spline` | substantial |
| `nsp` | `splines2` | `nsp` | substantial |
| `nsk` | `splines2` | `nsk` | substantial |
| `bernsteinPoly` | `splines2` | `bernstein_poly` | substantial |
| `bpoly` | `splines2` | `bpoly` | substantial |

The exact export inventory comes from upstream `NAMESPACE` and was independently
checked with `getNamespaceExports("splines2")`. Registered S3 methods for
`deriv`, `predict`, `update`, `knots`, `$`, and `makepredictcall` are not namespace
exports and are excluded from this denominator. Plotting/printing methods are
presentation-only. Typed prediction and differentiation are provided as useful
Fortran APIs; no additional coverage entries are claimed for unexported methods.

R classes, attributes, dimnames, coercion/recycling, formula/model-frame handling,
warning emission, and generic dispatch are omitted. NaN represents R missing
numeric values, without preserving a distinction between NA and NaN. The
Fortran derivative helper also permits a zero additional derivative order.
See [API_COVERAGE.md](API_COVERAGE.md) for exact scope and differences.

## Validation

`fpm test` runs 1,190 deterministic cases generated by the actual compiled
upstream R package plus 231 independent identity/API checks. No R installation
is needed to run the committed tests. Reference generation is documented in
[VALIDATION.md](VALIDATION.md); [NOTICE.md](NOTICE.md) records source pins,
licenses, and dependency reuse decisions.
