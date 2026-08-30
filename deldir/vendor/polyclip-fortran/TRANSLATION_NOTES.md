# Translation notes

## Upstream target

- R package: `polyclip` 1.10-7 (2024-07-23)
- Package note: built from Clipper C++ 6.4.0
- License: Boost Software License 1.0

The complete source snapshot used for the translation is retained in
`upstream/`.

## Translation strategy

The offsetting formulas and Minkowski construction follow the retained Clipper
implementation closely. Closed polygon Boolean operations are implemented in
native modern Fortran as a planar segment arrangement:

1. collect and split segments at intersections,
2. classify each atomic edge on both sides with the requested winding rules,
3. retain and orient boundary edges for the requested Boolean operation,
4. cancel duplicate/opposite edges, and
5. stitch the surviving edges into output paths.

Open subject paths are split at clip-boundary intersections and classified by
segment midpoint, preserving Clipper's open-path contribution semantics.

This gives a fully native Fortran implementation rather than a C++ wrapper.

## Coordinate scaling

The R wrapper scales coordinates to signed 64-bit integer coordinates. The
source code casts `(x-x0)/eps` and `(y-y0)/eps` to the integer type, which
truncates toward zero. The Fortran translation intentionally preserves this
source behavior even though user-facing package documentation sometimes uses
the word "round".

Default `eps`, `x0`, and `y0` choices follow the R-side wrappers, including the
special Minkowski origin convention.

## Deterministic parity checks

A standalone harness built from the retained Clipper 6.4.0 source was used as
a reference. The Fortran implementation was checked against it for:

- intersection, union, difference, and XOR of overlapping/touching/identical
  polygons,
- self-intersecting polygon simplification with even-odd and nonzero rules,
- positive and negative winding fill rules,
- open-line intersection and difference,
- positive and negative polygon offsets,
- polygons containing holes,
- square, miter, and round joins,
- open round line offsets,
- Minkowski sums, and
- point-in-polygon inside/boundary/outside classifications.

These deterministic fixtures match the Clipper geometry, including open-path
vertex preservation/traversal in the tested cases.

## Important parity qualification

The Boolean engine is not a source-for-source implementation of Clipper's
Vatti scanline topology. With pathological near-coincident intersections,
independently computed intersections can quantize to neighboring integer grid
points differently from Clipper's scanline intersection ordering. Randomized
stress tests found occasional one-grid-unit/tiny-edge differences in such
cases.

At the R package's usual default scale, one integer grid unit is `eps`
(approximately the coordinate span divided by 1e9), so these differences are
at the selected geometric resolution. Nevertheless this release does **not**
claim bit-for-bit/topological identity with Clipper for every degenerate or
near-degenerate input.

## Maintained Fortran source rules

- `dp = real64` is defined once in `src/polyclip_kinds.f90` and re-exported by
  the public module.
- Maintained real variables use `real(dp)`.
- Maintained real literals use `_dp` suffixes where appropriate.
- No `double precision`, `real*8`, `kind(0.0d0)`, or D-exponent constants.
- `implicit none` and explicit module/procedure interfaces are used throughout.

## Release validation

The release tests exercise Boolean operations/fill rules, self-intersections,
open paths, offsets and holes, join/cap styles, Minkowski sums, and
point-in-polygon.

A strict GNU Fortran build uses:

```text
-std=f2018
-Wall
-Wextra
-Werror
-Wimplicit-interface
-fimplicit-none
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

FPM metadata is provided in `fpm.toml`.
