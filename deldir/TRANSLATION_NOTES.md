# Translation notes

## Upstream target

`deldir` 2.0-4, dated 2024-02-27.

## Directly modernized upstream kernel

The maintained `src/deldir_kernel.f90` is derived from all Fortran routines under upstream `src/`, including `master`, incremental point insertion, adjacency maintenance, Lee-Schacter local optimality tests, Delaunay/Dirichlet segment output, circumcentres, rectangular-boundary intersection, bin sorting, triangle membership, Stokes/polygon area support, and mean nearest-neighbour distance.

The transformation preserves the algorithm while replacing implicit typing and R-runtime diagnostics. On a deterministic five-point fixture, the modern kernel and retained original Fortran produced identical formatted numerical output for all Delaunay segments, all Dirichlet segments, segment/status counts, convex-hull area, and rectangular-window area.

## R-side computational code translated

The native API translates or subsumes the computational content of:

- `deldir()`
- `binsrtR()`
- `duplicatedxy()`
- `tile.list()`, including v0.2.0 `clipp=` arbitrary-polygon clipping
- `doClip()` through `polyclip-fortran`
- `tileArea()` / `tilePerim()` / `tilePerim0()`
- `tile.centroids()`, including native multi-component/hole support
- `cvt()`; v0.2.0 additionally supports persistent polygon clipping through Lloyd iterations
- `triang.list()` / `triMat()`
- `getNbrs()`
- `which.tile()`
- `insideRect()` / `insidePoly()`
- `mnndR()`
- `divchain()` for native integer tags
- `tileInfo()` computational summaries, including `clipp=`
- `lawSummary()`
- `findNewInOld()`, `get.cnrind()`, and `mid.in()` equivalents

R-specific identifiers/tags are represented through native integer indices or explicit integer tag arrays rather than arbitrary R objects.

## Polygon clipping in v0.2.0

Upstream does not replace the Delaunay rectangular window with the arbitrary polygon. It first creates the normal tile list and then calls `polyclip::polyclip(tile, clipp)`. The Fortran implementation follows this ordering:

1. construct each Voronoi cell inside `result%rw` with the existing half-plane path;
2. intersect that cell with the supplied `poly_set`;
3. retain inherited rectangular-boundary flags for vertices found in the old tile;
4. drop cells with empty intersections;
5. represent disconnected intersections as multiple `voronoi_component` objects.

Even-odd clipping supports nested paths/holes. Component signed areas are retained so `tile_area()` and `tile_centroid()` handle holes correctly. `deldir_tile_info()` combines edges across all components, corresponding to upstream `tileInfo()`'s `getEdges()` logic.

The translated `polyclip-fortran` v0.1.0 package is vendored as an FPM path dependency. Its maintained source is a separate package and retains its own kind module and BSL-1.0 licensing/provenance.

## Not translated

- plotting and printing methods;
- S3/list/data-frame/formula dispatch;
- R package startup/version helpers;
- colour utilities.

## Numerical kinds

`dp = real64` is defined exactly once in the maintained deldir source, in `deldir_kinds`, and re-exported from `deldir`. Maintained deldir sources use `real(dp)` and `_dp` real literals. Upstream source and the vendored dependency are separate source trees and are not part of the deldir source-hygiene count.

## Parity qualification inherited from polyclip-fortran

The translated Boolean polygon engine is not source-for-source Clipper/Vatti. Pathological near-coincident intersections may occasionally move an intersection by one integer clipping-grid unit. This qualification applies only when `clipp=` is used; the original triangulation and rectangular Voronoi path do not depend on it.
