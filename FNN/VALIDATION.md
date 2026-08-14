# Validation

Validation was performed with GNU Fortran 14.2.0 on the translated source.

## Strict compilation

All library modules, all permanent tests, and the example compile with:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all
```

No implicit interfaces, bounds errors, allocation errors, or runtime-checking
warnings remain in the final pass.

## Permanent tests

- `test_core`: deterministic checks for kd/cover/brute search, external-query
  search, correlation distance, duplicated observations, information measures,
  classification, leave-one-out classification, regression, explicit OWNN,
  automatic-k OWNN, and optional accuracy output.
- `test_randomized`: 200 deterministic randomized data sets; exact kd-tree and
  cover-tree indices/distances are compared against the independent brute-force
  kernel for both self and external-query search.
- `test_information_reference`: entropy, cross-entropy, KL divergence and KSG
  mutual information are compared against independently calculated SciPy
  reference values.  The tolerance is `3e-11`.

All tests pass.

## Development stress/benchmark checks

On a 2,000 x 5 random data set with k=10 in the validation container, an
optimized development build produced identical neighbor indices/distances for
all methods.  Approximate CPU times were 0.228 s for brute force, 0.066 s for
the kd-tree, and 0.120 s for the cover hierarchy.  These timings are included
only as a sanity check, not as portable performance claims.

## FPM

The validation environment did not contain the `fpm` executable.  `fpm.toml`
was parsed separately and the exact FPM source/test/example layout was compiled
directly with gfortran.
