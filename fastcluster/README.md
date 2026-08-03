# fastcluster-fortran

Modern Fortran translation of the computational interface of `fastcluster`
1.3.0, packaged for the Fortran Package Manager (FPM).

The library performs hierarchical agglomerative clustering from either a
precomputed dissimilarity matrix or a matrix of observations. It preserves the
R `hclust` merge convention: negative merge entries identify original
observations and positive entries identify earlier merge rows.

## Implemented functionality

- `hclust` from a full symmetric distance matrix
- `hclust` from an R-compatible condensed distance vector
- `hclust_vector` from observations stored by rows
- Linkage methods:
  - single
  - complete
  - average/UPGMA
  - McQuitty/weighted/WPGMA
  - Ward.D
  - Ward.D2
  - centroid/UPGMC
  - median/WPGMC
- Vector methods:
  - single linkage with Euclidean, maximum, Manhattan, Canberra, binary, or
    Minkowski distance
  - Ward, centroid, and median linkage with Euclidean distance
- Optional initial cluster membership weights
- Missing-coordinate handling compatible with R's `dist` scaling rules
- Infinity-safe validation and explicit status reporting

## Build

```text
fpm build
fpm test
fpm run demo_fastcluster
```

The project has no external library dependencies.

## Minimal example

```fortran
use fastcluster
implicit none

real(dp) :: x(4, 2)
type(hclust_result) :: fit

x = reshape([0.0_dp, 1.0_dp, 0.0_dp, 5.0_dp, &
             0.0_dp, 0.0_dp, 2.0_dp, 5.0_dp], shape(x))

call hclust_vector(x, 'ward', fit)
if (.not. fit%ok()) error stop fit%message

print *, fit%merge
print *, fit%height
print *, fit%order
```

## Numerical implementation

The upstream C++ package selects specialized minimum-spanning-tree,
nearest-neighbor-chain, and generic algorithms. This portable Fortran version
uses a common Lance-Williams engine. It has `O(n^2)` memory use and `O(n^3)`
worst-case running time. The linkage formulas and output representation are
preserved, but very large clustering jobs will be slower than the upstream
`fastcluster` implementation.

See `PORTING.md`, `API.md`, and `TRANSLATION_COVERAGE.md` for details.

## License

The upstream package is offered under a BSD-2-Clause/FreeBSD or GPL-2 option.
This translation is distributed under BSD-2-Clause. The upstream license,
GPL-2 text, metadata, R interface, C++ implementation, documentation, and tests
are retained under `original/` and at the project root.
