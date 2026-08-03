# Testing

## FPM

```text
fpm test
fpm run demo_fastcluster
fpm run --example example_distance_matrix
fpm run --example example_vector_metrics
fpm run --example example_linkage_methods
```

## Direct GNU Fortran validation

Run:

```text
sh scripts/test_gfortran.sh
sh scripts/test_gfortran_optimized.sh
```

Windows command scripts are also supplied.

## Test coverage

- all six vector metrics
- R-compatible missing-coordinate scaling
- condensed/full distance conversion
- single, complete, average, McQuitty, and Ward.D2 distance clustering
- Ward, centroid, and median vector clustering
- optional membership weights
- merge numbering and dendrogram leaf order
- invalid methods, asymmetric matrices, NaNs, invalid members, and invalid
  method/metric combinations

Random tie-free validation during translation compared vector single, Ward,
centroid, and median merge sequences and heights with SciPy 1.17.0 for data
sets of dimensions 7-by-3 and 9-by-5.
