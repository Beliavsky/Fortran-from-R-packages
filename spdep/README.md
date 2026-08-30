# spdep — modern Fortran computational translation

This directory contains a modern free-form Fortran translation of major
computational functionality from R package `spdep` 1.4-2. It is intended to be
placed as the top-level `spdep` directory in
`Beliavsky/Fortran-from-R-packages`.

The translation keeps the upstream GPL-2-or-later licensing and attribution,
uses FPM, contains no C or R runtime dependency, and does not require system
BLAS/LAPACK libraries. It uses one real kind (`dp = real64`) throughout.

## Build and test

```text
fpm build
fpm test
```

Run the example with:

```text
fpm run --example moran_example
```

## Public data types

- `neighbor_list`: compact one-based neighbor vectors
- `spatial_weights`: neighbor topology plus one real weight vector per region
- `knn_result`: KNN indices and distances
- `weights_constants`: `n`, `S0`, `S1`, `S2`
- `spatial_test_result`: statistic, expectation, variance, z score, p value
- `local_stat_result`: vector-valued local statistic results
- `eb_result`: raw and empirical-Bayes rates
- `mst_result`: minimum-spanning-tree edge list and costs
- `spatial_delta_result`: Bavaud spatial-delta moments and inference

## Small example

```fortran
program demo
   use spdep, only : dp, neighbor_list, spatial_weights, spatial_test_result, &
      cell2nb, nb2listw, moran_test
   implicit none

   type(neighbor_list) :: nb
   type(spatial_weights) :: listw
   type(spatial_test_result) :: res
   real(dp) :: x(6)

   nb = cell2nb(1, 6)
   listw = nb2listw(nb, "W")
   x = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 8.0_dp, 13.0_dp]
   res = moran_test(x, listw)

   print *, res%statistic, res%p_value
end program demo
```

See `API_COVERAGE.md` for the translation boundary and `PROVENANCE.md` for
algorithm/source mappings. Exact upstream metadata and citations are retained
in `UPSTREAM_DESCRIPTION` and `UPSTREAM_CITATION.R`.

## Numerical and portability notes

The package uses no `-ffast-math`, `-Ofast`, `-ffinite-math-only`, or related
finite-only assumptions. Invalid numerical results are represented with IEEE
quiet NaNs; explicit NaN tests in maintained source use `ieee_is_nan`.

KNN, graph shortest paths, distance-decay weights/autocovariates, graph
complements and cumulative lags, IPFP balancing, and the small spectral
calculation needed by linearized diffusive weights are implemented directly in
free-form Fortran. This makes the current package self-contained and avoids copying or
vendoring packages already represented elsewhere in the target repository.
