# Verification

Run the Fortran tests with strict GNU Fortran diagnostics and runtime checks:

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
```

Run the independent R behavior checks with an installed tidyr 1.3.2 or later:

```text
Rscript reference/reference_invariants.R
```

The Fortran tests cover both pivot orderings, a round trip from wide to long and
back, absent wide cells, missing-value transformations, row expansion, literal
string separation and union, Cartesian products, type preservation, and masks.

`test_tidyverse_pipeline` reads a CSV through readr and then exercises tibble,
tidyselect, dplyr, tidyr, and stringr in one pipeline.
