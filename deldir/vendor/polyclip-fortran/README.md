# polyclip-fortran

Modern free-form Fortran translation of the computational geometry in the R
package `polyclip` 1.10-7, which is based on Angus Johnson's Clipper library.
The project is packaged for the Fortran Package Manager (FPM).

## Implemented computational API

The public module is `polyclip` and re-exports the package real kind `dp`.

- `polyclip_apply` -- intersection, union, difference, and XOR
- `polysimplify` -- resolve self-intersections using a selected fill rule
- `polyoffset` -- offsets/buffers for closed polygons
- `polylineoffset` -- offsets/buffers for open or closed polygonal lines
- `polyminkowski` -- Minkowski sum
- `pointinpolygon` -- outside/inside/boundary classification

Boolean filling supports even-odd, nonzero, positive, and negative winding
rules. Offsetting supports square, round, and miter joins and the Clipper end
styles for closed/open lines.

The lightweight public geometry containers are `poly_path` and `poly_set`.
Constants for clip operations, fill rules, joins, and end types are exported
from `use polyclip`.

The module is named `polyclip`, so the clipping routine is named
`polyclip_apply` rather than `polyclip`.

## Numeric representation

Like the R package interface, input coordinates are mapped to signed 64-bit
integer coordinates before geometry operations and mapped back afterwards.
If scaling parameters are omitted, defaults follow the R wrappers. In
particular, the upstream C++ interface converts scaled floating-point
coordinates to integer coordinates by truncation toward zero; this translation
preserves that source behavior.

All maintained floating-point source uses the single public kind
`dp = real64` defined in `polyclip_kinds`.

## Example

```fortran
program example
   use polyclip, only: dp, poly_set, make_path, make_set, &
      polyclip_apply, clip_intersection
   implicit none

   type(poly_set) :: a, b, answer

   a = make_set([make_path([0.0_dp, 4.0_dp, 4.0_dp, 0.0_dp], &
                           [0.0_dp, 0.0_dp, 4.0_dp, 4.0_dp])])
   b = make_set([make_path([2.0_dp, 6.0_dp, 6.0_dp, 2.0_dp], &
                           [2.0_dp, 2.0_dp, 6.0_dp, 6.0_dp])])

   call polyclip_apply(a, b, answer, op=clip_intersection)
end program example
```

A complete executable example is in `example/basic_example.f90`.

## Build and test with FPM

```text
fpm build
fpm test
fpm run --example basic_example
```

The release was also tested directly with GNU Fortran using strict Fortran
2018, warnings-as-errors, interface checking, bounds checking, and
floating-point traps; see `TRANSLATION_NOTES.md`.

## Scope

The maintained Fortran project translates the numerical geometry, not the R
list/S3/native-registration interface. The full upstream package source is
retained under `upstream/`.

See `TRANSLATION_NOTES.md` for parity details and the known near-coincident
intersection qualification.
