# Verification

The deterministic Fortran suite covers fixed matching, bounds, counts,
locations, positive and negative substrings, ASCII case conversion, whitespace,
padding, truncation, literal replacement, wildcards, identifier case, joining,
ordering, uniqueness, subsets, and splitting.

Run strict checks with:

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all"
```

Independent base-R and optional installed-stringr checks are provided by:

```text
Rscript reference/reference_invariants.R
```
