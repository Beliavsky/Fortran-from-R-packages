# API

All public symbols are available from the `fastcluster` module.

## Kinds and status codes

- `dp`: double-precision real kind (`kind(1.0d0)`)
- `fc_success`
- `fc_invalid_argument`
- `fc_nan_distance`
- `fc_numerical_failure`
- `fc_allocation_failure`

## `hclust_result`

```fortran
type :: hclust_result
  integer :: n
  integer :: status
  character(len=:), allocatable :: message
  character(len=:), allocatable :: method
  character(len=:), allocatable :: metric
  integer, allocatable :: merge(:, :)
  real(dp), allocatable :: height(:)
  integer, allocatable :: order(:)
contains
  procedure :: ok
end type hclust_result
```

For `n` observations, `merge` has shape `(n-1,2)`, `height` has length `n-1`,
and `order` has length `n`.

Merge entries follow R's `hclust` convention:

- `-i` means original observation `i`.
- `j > 0` means the cluster formed in merge row `j`.

## `hclust`

Generic interface covering full and condensed distance input.

### Full matrix

```fortran
call hclust(distances, method, result [, members])
```

`distances` is a symmetric `n x n` matrix. The diagonal is ignored.

### Condensed distances

```fortran
call hclust(d, n, method, result [, members])
```

`d` has `n*(n-1)/2` values in R `dist` order:

```text
d(2,1), d(3,1), ..., d(n,1), d(3,2), ..., d(n,n-1)
```

Accepted method names are:

```text
single, complete, average, mcquitty, weighted,
ward, ward.D, ward.D2, centroid, median
```

`ward` is an alias for `ward.D`, matching the R wrapper.

`members`, when present, contains positive initial cluster sizes. It affects
average, Ward, and centroid updates.

The specific procedures `hclust_matrix` and `hclust_condensed` are also public.

## `hclust_vector`

```fortran
call hclust_vector(x, method, result [, members] [, metric] [, p])
```

Rows of `x` are observations and columns are variables.

Accepted methods are:

```text
single, ward, centroid, median
```

For single linkage, accepted metrics are:

```text
euclidean, maximum, manhattan, canberra, binary, minkowski
```

`p` is required conceptually for Minkowski distance and defaults to `2` when
omitted. Ward, centroid, and median require Euclidean distance.

## Distance helpers

```fortran
call pairwise_distances(x, metric, distances [, p] [, status] [, message])
call condensed_to_matrix(d, n, distances [, status] [, message])
call matrix_to_condensed(distances, d [, status] [, message])
```

`pairwise_distances` also has an optional `squared_euclidean` argument used by
the clustering implementation. It must not be enabled for non-Euclidean
metrics.

## Missing values and infinities

IEEE NaN coordinates are omitted pairwise. Euclidean, Manhattan, Canberra, and
Minkowski distances are rescaled when only a subset of coordinates is usable,
matching R `dist`. A pair with no usable coordinates returns
`fc_nan_distance`.

Positive infinity is accepted as a dissimilarity. NaN dissimilarities are
rejected. Methods whose extended-real Lance-Williams expression would be
indeterminate retain an infinite updated distance rather than executing an
invalid floating-point operation.
