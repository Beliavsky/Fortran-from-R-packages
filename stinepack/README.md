# stinepack

Modern free-form Fortran translation of the computational code in the R package
`stinepack` 1.5, *Stineman, a Consistently Well Behaved Method of Interpolation*.
The package implements Stineman interpolation, the original and scaled Stineman
slope estimators, parabola-based slope estimation, and NaN gap filling for vectors
and matrices.

This is an unofficial translation for the
[Beliavsky/Fortran-from-R-packages](https://github.com/Beliavsky/Fortran-from-R-packages)
collection. It is not endorsed by the upstream authors, CRAN, or the R Foundation.

## Build

From the extracted top-level `stinepack` directory:

```text
fpm build
fpm test
fpm run --example interpolate_sine
```

The code has no external package or BLAS/LAPACK dependency. It uses one package
working precision, `dp = real64`, and only free-form Fortran.

## Public API

```fortran
use stinepack_api, only : dp, parabola_slopes, stineman_slopes, &
                          stinterp, stinterp_result, na_stinterp
```

- `parabola_slopes(x, y)` estimates knot slopes from parabolas through triples
  of points.
- `stineman_slopes(x, y, scale)` implements the original Stineman arc-based
  slope estimator, optionally after range scaling.
- `stinterp(x, y, xout, yp, method)` evaluates the piecewise rational Stineman
  interpolant. The optional `method` accepts a unique prefix of
  `scaledstineman`, `stineman`, or `parabola`, following the relevant behavior
  of R's `match.arg`. Explicit `yp` slopes and `method` are mutually exclusive.
- `na_stinterp` is a rank-generic interface. For a vector it fills NaNs and can
  remove unresolved edge NaNs. For a matrix it interpolates each column and,
  when `na_rm=.true.`, removes rows that still contain a NaN.

`stinterp` returns `type(stinterp_result)`, containing `x`, `y`, `status`, and
`message`. This represents the R result list while avoiding termination from a
library routine. Out-of-domain interpolation values are quiet NaNs.
`vector_interp_result` and `matrix_interp_result` similarly carry gap-filled
values and diagnostics.

Plain vectors and matrices default to coordinates `1, 2, ..., n`, matching the
upstream documented `time.default` behavior. R time-series/class attributes are
not represented; callers can pass an explicit `along` coordinate vector.

## Numerical behavior

The interpolation formulas follow `upstream/R/stinterp.R` directly, including:

- strictly increasing knots;
- supplied or estimated slopes;
- scaled Stineman, unscaled Stineman, and parabola slope methods;
- R-style unique-prefix matching for method names;
- the endpoint extrapolation tolerance of five machine epsilons times the knot
  range;
- quiet NaNs for requested points outside the accepted domain.

NaN checks use `ieee_is_nan`/IEEE arithmetic. Do not build this package with
`-ffast-math`, `-Ofast`, `-ffinite-math-only`, or equivalent options when NaN
behavior is required.

## Validation

Deterministic Fortran tests cover the three slope/interpolation methods,
explicit slopes, two-point behavior, method-prefix validation, endpoint
handling, vector gap filling, matrix gap filling, and `na_rm` behavior.

The implementation was additionally compared against an independent direct
translation of the retained upstream R formulas on deterministic randomized
fixtures:

- 120 interpolation cases across the three estimated-slope methods, 1,618
  finite evaluated values: maximum absolute difference `0.0` in the test
  environment;
- 60 interpolation cases with explicitly supplied slopes, 674 finite evaluated
  values: maximum absolute difference `0.0` in the test environment.

These comparisons check the numerical formulas, not R object/class semantics.
See `VALIDATION.md` for the build and validation notes.

## Translation coverage

Package status: **substantial**. Coverage is **4 of 4 (100%)** exported
computational R functions. This fraction measures computational functions with
a meaningful Fortran mapping; it does **not** claim complete compatibility with
R's object model, S3 dispatch, argument coercion, attributes, or `...` behavior.
The dispatch-only R generic `na.stinterp` is not counted; its explicitly
registered computational default method is counted.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `na.stinterp.default` | `stinepack_api` | `na_stinterp` | substantial |
| `parabolaSlopes` | `stinepack_api` | `parabola_slopes` | complete |
| `stinemanSlopes` | `stinepack_api` | `stineman_slopes` | complete |
| `stinterp` | `stinepack_api` | `stinterp` | substantial |

The machine-readable mapping is in `[extra.translation]` and its
`[[extra.translation.function]]` entries in `fpm.toml`. Additional details are
in `API_COVERAGE.md`.

## Differences from the R package

The computational formulas are retained, but several R-specific facilities are
intentionally outside this port:

- S3 dispatch and arbitrary class preservation;
- `ts`, `zoo`, `its`, data-frame, and `time()` attribute semantics;
- general R coercion, recycling, named/partial argument matching, and arbitrary
  `...` forwarding;
- R's `stop()` exception mechanism. Fortran routines report invalid inputs with
  `status` and `message` fields where validation is part of the public API.

The Fortran `na_stinterp` interface directly supports numeric rank-1 and rank-2
arrays, which covers the upstream vector/matrix numerical operation.

## License and provenance

The upstream R package declares `GPL (>= 2)`. This translation is distributed
under GPL-2.0-or-later. `LICENSE` contains the GNU GPL version 2 text, while
`NOTICE.md` and `PROVENANCE.md` record upstream authorship and source lineage.
Retained upstream metadata and R computational source files are under
`upstream/` for traceability; they are not vendored Fortran dependencies.
