# Verification

Run strict GNU Fortran tests with:

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
```

Run independent behavior checks against installed forcats 1.0.1 with:

```text
Rscript reference/reference_invariants.R
```

The Fortran tests cover construction, explicit missingness, level code
remapping, collapse behavior, stable ordering and ties, level-set operations,
concatenation, crossing, matching, counts, proportions, lumping, and numeric
summary reordering.
