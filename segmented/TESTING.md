# Testing

Four test programs are included:

- `test_segmented_lm`: two-breakpoint continuous regression, prediction, slopes,
  confidence intervals, AAPC, and broken-line evaluation.
- `test_stepmented_glm`: exact mean-shift estimation plus binomial and Poisson
  segmented models.
- `test_inference_selection`: Davies/score tests, power, and BIC breakpoint
  selection including the zero-break model.
- `test_segmented_lme`: random-intercept segmented mixed model with BLUP output.

The GNU Fortran validation scripts use warnings as errors. Debug mode adds bounds,
allocation, array-temporary, and floating-point runtime checks. Release mode uses
`-O3` while retaining warning-as-error compilation.

```text
./run_gfortran_tests.sh debug
./run_gfortran_tests.sh release
```
