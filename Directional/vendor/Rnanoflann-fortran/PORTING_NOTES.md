# Porting notes

## Source basis

This port was made from the attached CRAN-style Rnanoflann 0.0.3 source tree. The numerical definitions are taken from `src/knn.cpp`, `R/nn.R`, and `inst/include/internal/dists.hpp` rather than from the prose formulas when those differ.

## Public API mapping

The R package has one exported R function:

- `nn()` -> Fortran `nn()` returning `type(nn_result)`

The Fortran port additionally exposes `metric_distance()` and `metric_code()` because direct numerical distance access is useful outside R's wrapper model.

## Matrix orientation

Rnanoflann accepts observations in rows, transposes them for Armadillo/nanoflann, and by default transposes the returned `k x n_query` arrays back to `n_query x k`. The Fortran API accepts observations in rows directly and reproduces the final R-facing orientation when `trans=.true.`.

## Exact source formulas retained

Several names are not conventional mathematical distances, but the port preserves the implementation:

- `cosine` returns cosine similarity, not `1-cosine`; nearest-neighbor search therefore minimizes the similarity just as upstream does.
- `harmonic_mean` likewise returns `2 dot(x,y) / sum(x+y)` and is minimized.
- `kullback_leibler` computes `(y-x)*(log(y)-log(x))` summed over finite terms, i.e. a symmetric Jeffreys-type divergence rather than one-sided KL.
- `itakura_saito` is directional and preserves upstream's data/query ordering `x/y - log(x/y) - 1`.
- `sorensen` is the sum of coordinate-wise ratios `abs(x-y)/(x+y)`, matching the C++ implementation.
- `motyka` is `1 - sum(min(x,y))/sum(x+y)`, matching the C++ implementation.

## Corrections to upstream search behavior

### Radius query pointer

In both serial and parallel branches of upstream `nn_helper`, radius search calls:

```cpp
mat_index.index_->radiusSearch(points.memptr(), ...)
```

inside the loop over query columns. This always passes the first query point. Standard search correctly uses `points.colptr(i)`. The Fortran port uses each query row for radius search.

### Radius output overflow

The upstream result arrays have exactly `k` slots per query, but `radiusSearch()` has no `k` cap. The code copies every returned result into those arrays. If more than `k` points are inside the radius, this writes past the allocated column. The Fortran port stores at most `k` results and reports the total number found in `result%counts`.

### `sorted` and `eps`

The C++ code creates `SearchParameters(eps, sorted)` but never passes that object to `knnSearch()` or `radiusSearch()`. Consequently `eps` is ignored and radius search uses nanoflann's default `sorted=.true.` regardless of the R argument. The Fortran port keeps exact search (so `eps` remains compatibility-only) but honors `sorted` for radius output.

## Why the active search kernel is not a copied KD-tree

Every metric adaptor in `dists.hpp`, including the Euclidean adaptor, implements `accum_dist()` as zero. Nanoflann uses `accum_dist()` to accumulate coordinate-wise lower bounds for branch pruning. With a zero contribution, the custom tree cannot establish the usual metric lower bound and exact searches traverse branches that a standard L1/L2 adaptor could prune.

The Fortran v0.1 implementation therefore computes the exact objective directly and keeps a bounded top-k result set. This preserves the numerical contract while removing the C++/Rcpp/Armadillo/nanoflann dependency. A future performance-oriented version could add native pruning-enabled KD-tree specializations for metrics with valid coordinate lower bounds without changing the public API.

## Hellinger preprocessing

`R/nn.R` replaces both matrices by their element-wise square roots before entering C++. The Fortran `nn()` and `metric_distance()` accept raw nonnegative values and perform the same transformation internally. Negative Hellinger input is rejected explicitly instead of propagating NaNs from `sqrt()`.

## Parallelism

The query loop is decorated with OpenMP directives. Without `-fopenmp` they are standard Fortran comments. With OpenMP enabled, different queries are independent and are processed in parallel. `cores=0` leaves thread selection to the OpenMP runtime.
