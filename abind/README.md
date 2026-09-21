# abind

Modern free-form Fortran translation of the computational code in the R package
`abind` 1.4-8. The upstream package combines and subsets multidimensional arrays;
this translation provides a general-rank numeric `array_value` type backed by
Fortran column-major storage so the same API works beyond fixed rank-2/rank-3
wrappers.

## Build

```text
fpm build
fpm test
fpm run --example basic
```

The package has no external dependencies and does not link BLAS or LAPACK.

## Public API

The facade module is `abind`.

- `array_value`: arbitrary-rank `real(dp)` data plus dimensions and optional labels.
- `abind_arrays`: combine arrays along an existing or newly inserted dimension.
- `asub_array`: subset arbitrary dimensions using resolved one-based indices.
- `afill_array`: fill an array using RHS dimension labels and explicit selectors.
- `adrop_array`: remove selected singleton dimensions.
- `acorn_array`: take leading or trailing slices from each dimension.
- `indices_from_names`, `indices_from_mask`: convert common R-style selectors to
  one-based integer positions.

All maintained real data use the single public kind `dp` defined from
`iso_fortran_env:real64` in `abind_types`.

## Translation coverage

Package status: **substantial**.

**5 of 5 (100.0%)** exported computational R functions are mapped. This fraction
measures functions with meaningful computational Fortran mappings; it does **not**
mean complete compatibility with R's S3/object/interface semantics.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
| --- | --- | --- | --- |
| `abind` | `abind_core` | `abind_arrays` | substantial |
| `asub` | `abind_core` | `asub_array` | substantial |
| `afill<-` | `abind_core` | `afill_array` | substantial |
| `adrop` | `abind_core` | `adrop_array` | substantial |
| `acorn` | `abind_core` | `acorn_array` | substantial |

The main differences from R are intentional interface differences: Fortran does
not reproduce S3 dispatch, data-frame coercion, expression-derived argument names,
R replacement-function syntax, or caller-environment mutation. Numeric indexing,
dimension insertion, array layout, name-driven filling, and rank-preserving corner
selection are implemented directly.

See `docs/API_COVERAGE.md` for details and `docs/VALIDATION.md` for validation.

## Upstream provenance

The original package is `abind` by Tony Plate and Richard Heiberger, distributed
under the MIT license. Original R sources and package metadata used for this
translation are retained under `upstream/` for provenance. See `NOTICE.md`.
