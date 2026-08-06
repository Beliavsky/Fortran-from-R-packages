# Testing

The test suite contains five programs:

1. `test_gradients`: finite-difference checks for exponential and Gamma
   gradients plus elastic-net penalty/proximal checks.
2. `test_exponential_fit`: recovery of a known exponential log-link model and
   verification that the fitted likelihood improves on zero coefficients.
3. `test_gamma_fit`: recovery of Gamma coefficients and shape, followed by a
   fixed-shape penalized fit.
4. `test_cross_validation`: repeated balanced K-fold selection and a complete
   regularization path.
5. `test_validation`: rejection of nonpositive responses and invalid alpha.

Run:

```text
make MODE=checked clean test
make MODE=optimized clean test
```

The example is exercised with:

```text
make MODE=checked example
```
