# Build report

Validated compiler: GNU Fortran 14.2.

Commands:

```sh
make clean
make test-check
make test-opt
make example
```

Both test configurations pass all six test programs. The example completes
with finite simulation, likelihood, dependence, CDF, and quantile output.

`-Wno-maybe-uninitialized` is included because GNU Fortran can emit false
positives for allocatable derived-type descriptors. All other enabled warnings
are errors. The checked configuration retains full runtime bounds and shape
checking.
