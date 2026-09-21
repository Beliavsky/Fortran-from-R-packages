# API coverage

Coverage basis: distinct exported computational R functions. Plotting, rendering, and presentation-only functions are excluded. Package status is **substantial** with **35 of 35 (100.0%)** mapped computational exports. A mapped function may still be partial with respect to R/Qhull interface behavior.

| R function | Fortran procedure(s) | Status | Compatibility notes |
|---|---|---|---|
| `entry.value<-` | `entry_set_real` | substantial | Flat real(dp) storage plus an explicit dimension vector replaces R replacement-function and arbitrary-type semantics. |
| `Unique` | `unique_rows` | complete | Numeric matrix row sorting and duplicate removal, including rows.are.sets behavior, are translated. |
| `bary2cart` | `bary2cart_points` | complete | Matrix barycentric-to-Cartesian conversion is translated. |
| `cart2bary` | `cart2bary_points` | substantial | Numerical conversion is translated; degeneracy is returned through a success flag instead of an R warning/NULL. |
| `cart2pol` | `cart2pol_points` | substantial | Matrix coordinate conversion is translated; R scalar recycling and argument-shape dispatch are omitted. |
| `cart2sph` | `cart2sph_points` | substantial | Matrix coordinate conversion is translated; R scalar recycling and argument-shape dispatch are omitted. |
| `convhulln` | `convhulln_points` | partial | Independent supporting-hyperplane implementation. Accurate 2D/3D area-volume paths are provided; Qhull option strings, file output, and Qhull-scale performance are omitted. |
| `delaunayn` | `delaunayn_points` | partial | Independent lifted-paraboloid lower-hull triangulation with deterministic tie breaking. Qhull options and large/high-dimensional performance are omitted. |
| `distmesh2d` | `distmesh2d_circle` | partial | DistMesh force relaxation is translated for a circular signed-distance region with uniform target edge length; R callback regions and plotting are omitted. |
| `distmeshnd` | `distmeshnd_sphere` | partial | N-dimensional force relaxation is translated for supplied start points in a spherical region with uniform target edge length; callback region/size functions are omitted. |
| `dot` | `dot_columns`, `dot_rows` | substantial | Row-wise and column-wise matrix dot products are translated; arbitrary-rank R apply semantics are omitted. |
| `entry.value` | `entry_value_real` | substantial | General-rank one-based subscript lookup is translated for flat real(dp) array storage. |
| `extprod3d` | `extprod3d_rows` | complete | Row-wise three-dimensional cross products are translated. |
| `feasible.point` | `feasible_point_hulls` | partial | Finds a common interior point by projection and halfspace vertices; lpSolve scaling-option retries are not reproduced. |
| `halfspacen` | `halfspacen_vertices` | substantial | Bounded halfspace vertices are computed by active-hyperplane enumeration; Qhull option strings and unbounded-region reporting are omitted. |
| `inhulln` | `inhulln_points` | substantial | Strict interior testing against translated hull facet inequalities is provided; it consumes the Fortran hull_result rather than an R object. |
| `intersectn` | `intersectn_points` | substantial | Convex-hull intersection is translated through facet inequalities and halfspace vertices; autoscale/Qhull-option behavior is omitted. |
| `matmax` | `matmax_rows` | complete | Numeric row maxima are translated. |
| `matmin` | `matmin_rows` | complete | Numeric row minima are translated. |
| `matorder` | `matorder_rows` | complete | Lexicographic numeric row ordering is translated. |
| `matsort` | `matsort_rows` | complete | Independent ascending sort of each numeric row is translated. |
| `mesh.dcircle` | `mesh_dcircle` | complete | Signed distance to a centered circle is translated. |
| `mesh.diff` | `mesh_diff_values` | substantial | The signed-distance difference operator max(dA,-dB) is translated for precomputed distance vectors; R function callbacks are omitted. |
| `mesh.drectangle` | `mesh_drectangle` | complete | Signed Euclidean distance to an axis-aligned rectangle is translated. |
| `mesh.dsphere` | `mesh_dsphere` | complete | Signed distance to a centered sphere/hypersphere is translated. |
| `mesh.hunif` | `mesh_hunif` | complete | Uniform desired edge length is translated. |
| `mesh.intersect` | `mesh_intersect_values` | substantial | The signed-distance intersection operator max(dA,dB) is translated for precomputed distance vectors; callbacks are omitted. |
| `mesh.union` | `mesh_union_values` | substantial | The signed-distance union operator min(dA,dB) is translated for precomputed distance vectors; callbacks are omitted. |
| `pol2cart` | `pol2cart_points` | substantial | Matrix coordinate conversion is translated; R scalar recycling and argument-shape dispatch are omitted. |
| `polyarea` | `polyarea_vector` | partial | Polygon-vector area is translated; R two-dimensional array/dimension dispatch is not exposed. |
| `rbox` | `rbox_points` | substantial | Hypercube corners and uniform random points are generated; random streams are not expected to match R. |
| `sph2cart` | `sph2cart_points` | substantial | Matrix coordinate conversion is translated; R scalar recycling and argument-shape dispatch are omitted. |
| `surf.tri` | `surf_triangles` | complete | Boundary triangles are extracted from tetrahedral meshes and oriented outward. |
| `tsearch` | `tsearch_points` | substantial | Triangle search and barycentric coordinates are translated by deterministic simplex scanning; quadtree acceleration is omitted. |
| `tsearchn` | `tsearchn_points` | substantial | N-dimensional simplex search and barycentric coordinates are translated; the experimental Qhull-backed R-object branch is omitted. |

## Excluded exported functions

`trimesh`, `tetramesh`, and `to.mesh3d` are plotting/rendering or presentation-oriented interfaces and are excluded from the computational denominator under the requested coverage convention. Registered plot methods are excluded for the same reason.
