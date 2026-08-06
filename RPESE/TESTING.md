# Testing

Five standalone test programs are included:

1. `test_measures`: deterministic checks for point estimates and compatibility modes.
2. `test_timeseries`: periodogram frequency selection, sinusoidal peak recovery, AR(1), and polynomial design.
3. `test_if_iid_matrix`: iid influence-function standard errors and matrix-column dispatch.
4. `test_correlated`: exponential/Gamma spectral fits, prewhitening, and adaptive blending.
5. `test_bootstrap_api`: iid/block bootstrap methods, named wrappers, and string parsers.

Run checked tests with:

```text
make MODE=checked test
```

Run optimized tests with:

```text
make clean
make MODE=optimized test
```
