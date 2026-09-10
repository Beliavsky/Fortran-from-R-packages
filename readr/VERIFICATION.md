# Verification

Run strict GNU Fortran tests with:

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror -fcheck=all -fbacktrace"
```

Run independent behavior checks against installed readr 2.2.0 with:

```text
Rscript reference/reference_invariants.R
```

The Fortran suite covers quotes, escaped quotes, embedded newlines, CRLF,
comments, ragged records, missing tokens, every implemented scalar parser,
numeric locales, factors, type inference, explicit specifications, structured
problem coordinates, output quoting, and CSV round trips.
