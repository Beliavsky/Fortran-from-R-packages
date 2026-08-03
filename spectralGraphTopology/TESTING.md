# Testing

Run all FPM tests:

```text
fpm test
```

Or use the supplied scripts:

```text
./run_gfortran_tests.sh
run_gfortran_tests.bat
```

The tests cover:

1. Graph operators, inverses, adjoints, matrix utilities, distances, block
   diagonal construction, and recovery metrics.
2. k-component and cospectral recovery on exactly generated covariance
   matrices, plus direct spectral-symmetry checks for bipartite and joint
   bipartite learning. Exact edge-support assertions are intentionally avoided
   because equivalent bipartite optima can vary with numerical platform.
3. Constrained Laplacian Rank, smooth graph learning, signal-representation
   learning, GLE-MM, GLE-ADMM, and CGL.
4. Eigenvector/eigenvalue updates, multi-component graphs, `nu=0` consistency,
   initial affinity construction, and the data-matrix input path.

The shell script compiles with GNU Fortran using strict diagnostics, runtime
bounds checking, and floating-point traps before running the tests, demo, and all examples. It also performs an
optimized build and smoke test.
