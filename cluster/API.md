# API

All observations are stored by rows and variables by columns. Real values use
`dp = kind(1.0d0)`.

## Result types

- `partition_result`: hard clustering, medoids, optional fuzzy memberships,
  objective, iteration count, and status.
- `hierarchy_result`: merge matrix, heights, leaf order, method, and hierarchy
  coefficient.
- `mona_result`: binary monothetic partition, order, and split variables.
- `silhouette_result`: individual widths, neighboring clusters, and average.
- `gap_result`: observed dispersion, gap values, standard errors, selected k.
- `ellipsoid_result`: center, covariance, inverse shape, volume, and status.

Each result type has an `ok()` method.

## Dissimilarities

### `daisy(x, distances [, metric, weights, status, message])`

Numeric distances. Metrics are `euclidean`, `manhattan`, `maximum`, `canberra`,
`binary`, `minkowski`, and `gower`. Missing real values represented by IEEE NaN
are omitted pairwise.

### `daisy_mixed(x, variable_types, distances [, weights, status, message])`

Gower-style mixed dissimilarity. Variable type constants are:

- `variable_numeric`
- `variable_binary_symmetric`
- `variable_binary_asymmetric`
- `variable_nominal`
- `variable_ordinal`

Nominal and binary variables are represented by numeric category codes.

## Partitioning

### `pam(x, k, result [, metric, max_iter])`
### `pam_distance(distances, k, result [, max_iter, initial_medoids])`

Partitioning around medoids.

### `clara(x, k, result [, samples, sample_size, metric, seed])`

Repeated sampled PAM. Sampling uses a portable deterministic generator when a
seed is supplied.

### `fanny(x, k, result [, membership_exponent, metric, tolerance, max_iter])`

Fuzzy c-medoids clustering. `membership_exponent` must exceed one.

## Hierarchical clustering

### `agnes(x, result [, method, metric])`
### `agnes_distance(distances, result [, method])`

Supported methods: `single`, `complete`, `average`, `weighted`, and `ward`.

### `diana(x, result [, metric])`
### `diana_distance(distances, result)`

Top-down divisive analysis using the splinter-group rule.

### `mona(x, result [, max_clusters])`

Monothetic analysis for integer 0/1 data.

### `coef_hier(result)`

Returns the stored agglomerative or divisive coefficient.

## Diagnostics and helpers

- `silhouette(labels, distances, result)`
- `sort_silhouette(labels, widths, order)`
- `medoids(labels, distances, indices [, objective, status])`
- `meanabsdev(x [, center])`
- `size_diss(length_value)`
- `lower_to_upper_tri_inds(n, indices)`
- `upper_to_lower_tri_inds(n, indices)`
- `within_cluster_dispersion(labels, distances)`
- `max_se(f, se [, method, se_factor])`
- `clus_gap(x, k_max, b_references, cluster_fun, result [, seed])`

A `clus_gap` callback has the interface:

```fortran
subroutine clustering_callback(x, k, labels, status)
  real(dp), intent(in) :: x(:, :)
  integer, intent(in) :: k
  integer, allocatable, intent(out) :: labels(:)
  integer, intent(out) :: status
end subroutine
```

## Ellipsoids

- `ellipsoidhull(x, result [, tolerance, max_iter])`
- `predict_ellipsoid(result, points, squared_distance [, inside, status])`
- `volume_ellipsoid(result)`
- `ellipsoid_points(result, n_points, points [, status])`

Boundary-point generation supports dimensions two and three.
