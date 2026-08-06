# Validation

Five deterministic test programs exercise the translated library.

1. `test_duration_convexity`
   - direct parity with the upstream formulas;
   - elemental array evaluation;
   - source-compatible zero-yield NaNs;
   - corrected analytic zero-yield limits.
2. `test_total_return`
   - exact one-period formula fixture;
   - default calculated sensitivities;
   - user-supplied duration and convexity;
   - corrected numerical mode.
3. `test_matrix`
   - independent column processing;
   - vector/matrix interface equivalence;
   - missing first row.
4. `test_preprocessing`
   - last observation carried forward;
   - preservation of leading missing values;
   - percentage conversion;
   - matrix preprocessing.
5. `test_validation`
   - invalid maturity and scale;
   - optional-array shape validation;
   - invalid yield-domain handling.

Both checked and optimized builds are required to pass before packaging. The
archive-validation process extracts the source-only ZIP into a separate
directory and reruns both suites and the example.
