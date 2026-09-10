# Verification

Run strict GNU Fortran tests with:

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
```

Run independent behavior checks against installed dplyr 1.2.1 with:

```text
Rscript reference/reference_invariants.R
```

The Fortran suite covers typed and string filtering, expression precedence,
arithmetic mutation, selection verbs, stable ordering, all implemented join
types, unmatched masks, set operations, stable grouping, built-in summaries,
and vector helpers.

