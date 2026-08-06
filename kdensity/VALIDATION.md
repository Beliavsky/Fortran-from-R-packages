# Validation

The deterministic test suite covers:

1. Symmetric and asymmetric kernel values.
2. Built-in parametric-start estimation.
3. nrd0, nrd, JH, RHE, and HS bandwidth selectors.
4. Normal-start Gaussian-kernel fitting, normalization, and infinite-bandwidth
   behavior.
5. Beta and gamma kernels on bounded and positive supports.

The beta-kernel example specifically checks the endpoint-singularity path that
requires transformed quadrature.

Validated compiler:

```text
GNU Fortran 14.2.0
```

Both checked and optimized builds pass with warnings treated as errors. Checked
builds use bounds, allocation, pointer, and floating-point runtime checks. The
Makefile uses heap trampolines for callback procedures, and the generated
objects do not request an executable stack.
