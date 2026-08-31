# Verification

Verification was performed on 2026-08-30 with GNU Fortran 14.2.0.

The library, all tests, and the example were compiled directly with strict
settings equivalent to the FPM source layout:

```text
-std=f2018
-Wall
-Wextra
-Werror=implicit-interface
-Werror=implicit-procedure
-fcheck=all
-fbacktrace
-O0
```

Six deterministic tests passed:

- `test_errors`: invalid confidence-level and imputation-index status handling.
- `test_mi_combine`: multivariate `MIcombine` with finite complete-data degrees
  of freedom, checked against independently calculated Rubin-rule references.
- `test_mi_summary`: standard errors and finite-df Student-t confidence limits.
- `test_mi_scalar`: scalar overload and the zero-between-imputation-variance
  limiting case.
- `test_imputation`: numeric imputation-list construction, extraction,
  dimensions, row binding, and column binding.
- `test_pv`: multi-variable plausible-value selection and materialization.

The example `mi_example` also compiled and ran successfully.

For local verification of the package source in this isolated runtime, an
interface-compatible temporary test shim was used outside the `mitools`
directory for the three `rfortran-core` procedures imported by this package.
The public signatures were checked against the current repository sources for
`r_kinds`, `r_descriptive`, and `r_distributions`. The shim is not included in
the package or ZIP. In the target repository, FPM resolves the declared sibling
`../rfortran-core` dependency instead.

The current execution environment does not contain the Fortran Package Manager
executable. The required `fpm build`, `fpm test`, and `fpm clean --all` commands
were explicitly attempted and each returned shell exit status 127 (`fpm: command
not found`). The exact transcript is retained in `FPM_ATTEMPTS.txt`. No substitute
or fake `fpm` executable is used.
Final static audit: 12 maintained Fortran files (library, tests, and example),
with zero policy issues and no duplicate source, build products, vendored
dependency source, or prohibited compiler/link settings.

