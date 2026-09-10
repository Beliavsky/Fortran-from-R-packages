# Verification

The deterministic FPM tests cover:

- all four atomic vector types and explicit missingness;
- logical/integer/real common types and promotion;
- checked casts and size-one recycling, including recycling to size zero;
- integer, real, and deferred-length character slicing and concatenation;
- stable sorting with missing values last;
- whole-vector and per-element repetition;
- unique values, duplicate detection, matching, and membership;
- dense group and consecutive-run identifiers;
- missing-aware equality and three-way comparison;
- type-stable conditional selection;
- compatible vector insertion at front, middle, and end positions;
- row binding with type promotion; and
- column binding with one-row recycling.

Run the ordinary test suite with:

```text
fpm test
```

The implementation was also verified under GNU Fortran with warnings promoted to
errors and runtime checking enabled:

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
```

`reference/reference_invariants.R` checks the corresponding upstream vctrs behavior
for numeric common types, recycling to size zero, missing-value equality, and matching:

```text
Rscript reference/reference_invariants.R
```
