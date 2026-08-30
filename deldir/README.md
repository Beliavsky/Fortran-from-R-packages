# deldir-fortran

Modern Fortran translation of the computational code in R package **deldir 2.0-4** by Rolf Turner. The package computes planar Delaunay triangulations and Dirichlet/Voronoi tessellations and related geometric summaries.

Version 0.2.0 adds native arbitrary-polygon clipping through the translated `polyclip-fortran` dependency.

## Build

With Fortran Package Manager:

```sh
fpm build
fpm test
fpm run --example basic
fpm run --example polygon_clip
```

The maintained deldir source is free-form Fortran 2018 and uses one package-wide real kind:

```fortran
use deldir, only: dp
```

`dp` is defined once from `iso_fortran_env::real64` in `deldir_kinds` and re-exported by the public `deldir` module. The vendored `polyclip-fortran` dependency is a separate package with its own `dp = real64` kind.

## Basic use

```fortran
use deldir
implicit none

real(dp) :: x(5) = [0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp, 0.5_dp]
real(dp) :: y(5) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp, 0.5_dp]
type(deldir_result) :: fit
type(voronoi_tile), allocatable :: tiles(:)

call deldir_compute(x, y, fit, rw=[0.0_dp, 1.0_dp, 0.0_dp, 1.0_dp])
call deldir_tiles(fit, tiles)
```

## Polygon-clipped tiles

The upstream R package implements `tile.list(..., clipp=...)` by calling the separate `polyclip` package after constructing the ordinary rectangular tessellation. Version 0.2.0 follows the same architecture.

```fortran
use deldir
implicit none
type(poly_set) :: clipp
type(voronoi_tile), allocatable :: tiles(:)

allocate(clipp%path(1))
clipp%path(1) = make_path([0.1_dp,0.9_dp,0.7_dp,0.3_dp], &
                          [0.1_dp,0.1_dp,0.9_dp,0.9_dp])
call deldir_tiles(fit, tiles, clipp=clipp)
```

`poly_path`, `poly_set`, `make_path`, `make_set`, `append_path`, and `clear_set` are re-exported from `deldir`, so users do not need a second `use polyclip` statement for ordinary clipping setup.

Clipped tiles with no intersection are removed, matching upstream `tile.list()`. A clipped cell may have more than one polygon component. `voronoi_tile%n_components()` reports the number of components. For one-component cells the historical `tile%x`, `tile%y`, and `tile%boundary_point` arrays remain populated. For multi-component cells use `tile%components(i)%x`, `%y`, and `%boundary_point`. Signed component areas allow holes to contribute with the correct sign.

The `bp`/`boundary_point` semantics match upstream `doClip()`: vertices inherited from the original rectangular-window tile retain their old boundary flag; vertices created only by the arbitrary clipping polygon are not marked as rectangular-boundary points.

## Native Fortran API

The public API includes:

- `deldir_compute`: Lee-Schacter Delaunay triangulation and rectangular-window Voronoi tessellation.
- `deldir_tiles`: explicit Voronoi polygons, optionally `clipp=` by an arbitrary `poly_set`.
- `deldir_triangles`, `deldir_triangle_matrix`.
- `deldir_get_neighbors`.
- `tile_area`, `polygon_area`, `polygon_signed_area`, `tile_perimeter`, `tile_centroid`, `tile_centroids`.
- `centroidal_voronoi`: Lloyd iteration; v0.2.0 can reapply a polygon `clipp=` on every iteration.
- `deldir_tile_info`: `tileInfo()` summaries, including clipped/multi-component cells.
- `mean_nearest_neighbor_distance` corresponding to `mnndR()`.
- `which_tile`, `inside_rect`, `inside_polygon`.
- `deldir_dividing_chain` corresponding to `divchain()` for integer tags.
- `deldir_law_summary` corresponding to `lawSummary()`.
- `deldir_bin_sort`, `duplicated_xy`, `corner_indices`, `midpoint_inside`, and `find_new_in_old`.

The primary `deldir_result` retains the package's Delaunay segments, Dirichlet segments, point summaries, convex-hull area, rectangular-window area, and original-index information in typed Fortran derived types.

## Translation approach

The original Lee-Schacter triangulation/tessellation kernel supplied by `deldir` was already Fortran and was modernized rather than replaced. The standard rectangular tile path continues to use direct half-plane clipping and therefore does not route ordinary calculations through `polyclip`. Only the optional arbitrary-polygon clipping layer uses `polyclip-fortran`.

The translated `polyclip-fortran` v0.1.0 project is vendored under `vendor/polyclip-fortran` so this release is self-contained for FPM builds.

## Numerical parity qualification

The underlying `polyclip-fortran` Boolean engine is a native planar-segment arrangement implementation rather than a source-for-source Vatti scanline translation. Normal and hard deterministic fixtures agree geometrically with Clipper 6.4.0, but pathological near-coincident intersections can occasionally differ by one integer clipping-grid unit. `deldir` polygon-clipping results inherit that narrow qualification. The Lee-Schacter triangulation and ordinary rectangular tessellation are unaffected.

## Deliberate omissions

Presentation and R-runtime functionality is not part of the native library: plotting methods, print/S3 methods, formula/data-frame handling, package-version helpers, and colour utilities are omitted.

## Licensing and provenance

The upstream `deldir` package declares `GPL (>= 2)`. This translation is distributed under **GPL-2.0-or-later**. See `COPYING`, `NOTICE.md`, and the retained `upstream/deldir-2.0-4/` sources. The vendored `polyclip-fortran` dependency is distributed under its upstream-compatible **BSL-1.0** terms; its `LICENSE` and provenance files are retained under `vendor/polyclip-fortran/`.
