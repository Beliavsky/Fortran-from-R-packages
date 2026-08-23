# Rnanoflann-fortran

Modern Fortran/FPM computational port of CRAN **Rnanoflann 0.0.3**.

The R package exports one public operation, `nn()`, which searches a data matrix for neighbors of each query point using a collection of distance/dissimilarity functions. This port exposes the same numerical operation through a Fortran-native API and also exposes the individual distance formulas for direct use.

## Quick start

```fortran
program demo
    use rnanoflann, only: dp, nn_result, nn
    implicit none
    real(dp) :: data(5,2), query(2,2)
    type(nn_result) :: ans

    data(1,:) = [0.0_dp, 0.0_dp]
    data(2,:) = [1.0_dp, 0.0_dp]
    data(3,:) = [0.0_dp, 1.0_dp]
    data(4,:) = [1.0_dp, 1.0_dp]
    data(5,:) = [2.0_dp, 2.0_dp]

    query(1,:) = [0.9_dp, 0.1_dp]
    query(2,:) = [1.8_dp, 1.9_dp]

    ans = nn(data, query, k=2, method="euclidean")
    print *, ans%indices
    print *, ans%distances
end program demo
```

Fortran matrices use the natural row-observation convention: `data(n_data,n_dim)` and `points(n_query,n_dim)`. The default `trans=.true.` therefore returns `indices(n_query,k)` and `distances(n_query,k)`, matching the default R-facing result after its transpose.

## Supported methods

All distance names implemented by the attached package are available:

- `euclidean`
- `hellinger`
- `manhattan`
- `canberra`
- `kullback_leibler`
- `jensen_shannon`
- `itakura_saito`
- `bhattacharyya`
- `jeffries_matusita`
- `minimum`
- `maximum`
- `total_variation`
- `sorensen`
- `cosine`
- `gower`
- `minkowski`
- `soergel`
- `kulczynski`
- `wave_hedges`
- `motyka`
- `harmonic_mean`

`metric_distance(x,y,method,...)` provides the same formulas independently of the neighbor search.

For `hellinger`, raw nonnegative vectors are accepted and square-root transformed internally, matching the R wrapper. `square=.true.` returns the squared Euclidean/Hellinger form. `minkowski` requires `p > 0`.

## Standard search

```fortran
ans = nn(data, query, k=5, method="manhattan")
```

Neighbor indices are 1-based and sorted by increasing value of the selected upstream metric, with index used as the deterministic tie breaker.

## Radius search

```fortran
ans = nn(data, query, k=10, method="euclidean", &
    search="radius", radius=0.5_dp, sorted=.true.)
```

The radius comparison follows nanoflann's strict `distance < radius` rule. At most `k` neighbors are stored. Empty/unused entries contain index 0 and distance `sqrt(huge(0.0_dp))`, corresponding to approximately `1.34078e154`, as documented by the R package. `ans%counts(i)` reports the total number found within the radius, which may exceed `k` when output is truncated.

## Search-engine note

The upstream package constructs a nanoflann KD-tree, but each of its custom metric adaptors returns zero from `accum_dist`. That provides no coordinate-wise lower bound for pruning. In addition, the package constructs `SearchParameters(eps,sorted)` but does not pass it to either `knnSearch` or `radiusSearch`, so `eps` and the caller's `sorted` option do not affect the C++ search as written.

This port therefore uses an exact native-Fortran scan/search kernel. It produces the same exact nearest-neighbor objective without retaining a C++ tree whose custom metric bounds are inactive. `eps` and `leafs` remain accepted for API/source compatibility. Radius `sorted` is intentionally honored in the Fortran version.

## Parallel queries

`parallel=.true.` enters an OpenMP-decorated query loop. The source remains valid without OpenMP, in which case those directives are comments and execution is serial. To enable OpenMP with GNU Fortran, compile with `-fopenmp`.

## Build

With FPM:

```text
fpm build
fpm test
fpm run --example basic
```

GNU Fortran validation for this release used:

```text
gfortran -std=f2018 -O2 -Wall -Wextra -Werror -fcheck=all
```

and a second OpenMP build used `-fopenmp`.

## License

The translated package code is **GPL-3.0-or-later**, following Rnanoflann's `License: GPL (>= 3)`. The upstream source snapshot also contains Jose Luis Blanco's/Muja/Lowe's BSD-licensed `nanoflann.hpp`; see `LICENSES.md` and the original package under `upstream/`.
