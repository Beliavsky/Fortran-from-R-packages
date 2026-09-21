# geometry

Modern free-form Fortran translation of the computational layer of the R package `geometry` 0.5.2.

The package provides coordinate transforms, matrix/array helpers, signed-distance primitives, barycentric conversions and simplex search, convex hulls, Delaunay triangulations, halfspace intersections, convex-hull intersection, boundary extraction from tetrahedral meshes, and restricted pure-Fortran DistMesh relaxation. The implementation is Fortran-only and does not vendor or link the Qhull C sources embedded by the upstream R package.

## Build

From this directory as a top-level repository package:

```text
fpm build
fpm test
fpm run --example geometry_example
```

No BLAS, LAPACK, ARPACK, R, Rcpp, Qhull, `magic`, `lpSolve`, or `linprog` installation is required by this translation.

## Computational-geometry scope

`convhulln_points`, `delaunayn_points`, and `halfspacen_vertices` are independent Fortran implementations, not bindings to Qhull. The hull routines use supporting-hyperplane enumeration; Delaunay triangulation uses lower facets of paraboloid-lifted points with deterministic nonlinear tie breaking; halfspace intersections enumerate active hyperplane sets. These choices are deterministic and useful for small and moderate problems, but their combinatorial cost is not suitable for the thousands-of-points workloads where Qhull is intended. 2D and 3D hull geometry receives the strongest compatibility testing.

`distmesh2d_circle` and `distmeshnd_sphere` retain the DistMesh edge-force relaxation but replace R callback distance/size functions with explicit circle/sphere plus uniform-size APIs so that all Fortran dummy arguments can carry standard `INTENT` attributes.

## Translation coverage

Package status: **substantial**.

**35 of 35 (100.0%)** exported computational R functions have a meaningful Fortran mapping. The fraction measures mapped computational functions, not complete R compatibility. Plotting/presentation exports `trimesh`, `tetramesh`, and `to.mesh3d` are excluded by the stated coverage convention. Qhull option strings, R object/S3 behavior, callback dispatch, file-output options, and exact RNG streams are not generally reproduced.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
| `entry.value<-` | `geometry_basic` | `entry_set_real` | substantial |
| `Unique` | `geometry_basic` | `unique_rows` | complete |
| `bary2cart` | `geometry_delaunay` | `bary2cart_points` | complete |
| `cart2bary` | `geometry_delaunay` | `cart2bary_points` | substantial |
| `cart2pol` | `geometry_basic` | `cart2pol_points` | substantial |
| `cart2sph` | `geometry_basic` | `cart2sph_points` | substantial |
| `convhulln` | `geometry_hull` | `convhulln_points` | partial |
| `delaunayn` | `geometry_delaunay` | `delaunayn_points` | partial |
| `distmesh2d` | `geometry_mesh` | `distmesh2d_circle` | partial |
| `distmeshnd` | `geometry_mesh` | `distmeshnd_sphere` | partial |
| `dot` | `geometry_basic` | `dot_columns`, `dot_rows` | substantial |
| `entry.value` | `geometry_basic` | `entry_value_real` | substantial |
| `extprod3d` | `geometry_basic` | `extprod3d_rows` | complete |
| `feasible.point` | `geometry_hull` | `feasible_point_hulls` | partial |
| `halfspacen` | `geometry_hull` | `halfspacen_vertices` | substantial |
| `inhulln` | `geometry_hull` | `inhulln_points` | substantial |
| `intersectn` | `geometry_hull` | `intersectn_points` | substantial |
| `matmax` | `geometry_basic` | `matmax_rows` | complete |
| `matmin` | `geometry_basic` | `matmin_rows` | complete |
| `matorder` | `geometry_basic` | `matorder_rows` | complete |
| `matsort` | `geometry_basic` | `matsort_rows` | complete |
| `mesh.dcircle` | `geometry_basic` | `mesh_dcircle` | complete |
| `mesh.diff` | `geometry_basic` | `mesh_diff_values` | substantial |
| `mesh.drectangle` | `geometry_basic` | `mesh_drectangle` | complete |
| `mesh.dsphere` | `geometry_basic` | `mesh_dsphere` | complete |
| `mesh.hunif` | `geometry_basic` | `mesh_hunif` | complete |
| `mesh.intersect` | `geometry_basic` | `mesh_intersect_values` | substantial |
| `mesh.union` | `geometry_basic` | `mesh_union_values` | substantial |
| `pol2cart` | `geometry_basic` | `pol2cart_points` | substantial |
| `polyarea` | `geometry_basic` | `polyarea_vector` | partial |
| `rbox` | `geometry_basic` | `rbox_points` | substantial |
| `sph2cart` | `geometry_basic` | `sph2cart_points` | substantial |
| `surf.tri` | `geometry_mesh` | `surf_triangles` | complete |
| `tsearch` | `geometry_delaunay` | `tsearch_points` | substantial |
| `tsearchn` | `geometry_delaunay` | `tsearchn_points` | substantial |

Detailed compatibility notes for every mapping are in `docs/API_COVERAGE.md`, and `fpm.toml` is the machine-readable source of truth for coverage.

## Provenance and license

The upstream R sources used for translation are retained under `upstream/R/` so their copyright and attribution comments remain visible. See `NOTICE.md`, `upstream/LICENSE-NOTES`, and `LICENSE`. This translation is distributed under GPL-3.0-or-later.
