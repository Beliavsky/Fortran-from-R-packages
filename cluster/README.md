# cluster-fortran

Modern Fortran 2018 translation of the computational core of R package
`cluster` 2.1.8.3, packaged for the Fortran Package Manager (FPM).

## Included algorithms

- `daisy` numeric distances and Gower-style mixed-variable dissimilarities
- `pam` partitioning around medoids with BUILD and SWAP phases
- `clara` deterministic sampled PAM for larger data sets
- `fanny` fuzzy medoid clustering
- `agnes` agglomerative hierarchical clustering
- `diana` divisive hierarchical clustering
- `mona` monothetic clustering for binary data
- silhouette widths and medoid extraction
- gap statistic with a typed clustering callback
- hierarchy coefficients and triangular-distance index helpers
- minimum-volume enclosing ellipsoids, prediction, volume, and boundary points

R plotting, S3 printing, menus, package datasets as live objects, and graphical
helpers are omitted. The original package sources and data remain under
`original/cluster/` for reference and license preservation.

## Build

```text
fpm build
fpm test
fpm run demo_cluster
```

GNU Fortran users can also run:

```text
scripts/test_gfortran.sh
scripts/test_gfortran_optimized.sh
```

On Windows, use the corresponding `.bat` scripts.

## Basic use

```fortran
use cluster, only: dp, partition_result, pam
real(dp) :: x(100, 4)
type(partition_result) :: fit

call pam(x, 3, fit)
if (.not. fit%ok()) error stop fit%message
print *, fit%medoids
print *, fit%clustering
```

See `API.md`, `PORTING.md`, and the programs in `example/`.
