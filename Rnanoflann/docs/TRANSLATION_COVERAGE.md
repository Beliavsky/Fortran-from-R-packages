# Translation coverage

## Exported R API

| Upstream export | Fortran | Status |
| --- | --- | --- |
| `nn` | `nn` | translated |

## Search modes

| Feature | Status |
| --- | --- |
| standard k-NN | translated |
| radius search | translated, with upstream pointer/overflow defects corrected |
| squared Euclidean | translated |
| squared Hellinger | translated |
| R-style output transpose | translated |
| OpenMP query parallelism | available when compiled with OpenMP |
| `eps` approximate search | compatibility argument only; upstream does not pass it to nanoflann |
| `leafs` tree tuning | compatibility argument only in the exact Fortran search kernel |

## Distance methods

All methods present in `src/knn.cpp` are translated: Euclidean, Hellinger, Manhattan, Canberra, Kullback-Leibler-labelled symmetric divergence, Jensen-Shannon, Itakura-Saito, Bhattacharyya, Jeffries-Matusita, minimum, maximum, total variation, Sorensen, cosine, Gower, Minkowski, Soergel, Kulczynski, Wave Hedges, Motyka, and harmonic mean.

## R/C++ infrastructure intentionally not reproduced

- Rcpp registration and `.Call` interface
- RcppArmadillo matrix wrappers
- `LinkingTo` export of the C++ nanoflann header
- R documentation/S3 conventions

The original package is preserved under `upstream/` for users who need those interfaces.
