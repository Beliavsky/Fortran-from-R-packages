# Validation

The translation was compiled and tested with GNU Fortran using:

```text
-std=f2008
-Wall -Wextra
-Wimplicit-interface -Werror=implicit-interface
-fcheck=all
-ffpe-trap=invalid,zero,overflow
```

The regression suite checks the upstream identities for both real and complex matrices, including:

- Hermitian transpose and conjugate cross-products
- quadratic and three-argument forms
- both upstream multiplication-order implementations
- transposed-orientation forms
- diagonal-only and trace-only calculations
- inverse forms against an independent 2x2 analytic solve
- Cholesky-factor evaluation
- singular-matrix reporting

Observed result:

```text
test_quadform: PASS
```

The example program also runs successfully.
