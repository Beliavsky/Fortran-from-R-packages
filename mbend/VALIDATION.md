# Validation

The final package manifest is configured for the repository layout required by
Fortran-from-R-packages:

```toml
[dependencies]
rfortran-linalg = { path = "../rfortran-linalg" }
```

The execution environment used for this translation did not contain an `fpm`
executable and had no outbound shell network access from which to install one.
Therefore the literal commands `fpm build` and `fpm test` could not be invoked in
this environment.

The maintained Fortran sources were instead compiled and exercised directly
with `gfortran` using:

```text
-std=f2018 -Wall -Wextra -Werror -fcheck=all
```

For local validation only, an API-compatible test shim supplied
`r_linalg:symmetric_eigen` and `symmetric_eigenvalues`; that shim is outside the
package and is not included in the ZIP. The production call signature was
checked against the current `rfortran-linalg` public API, whose implementation
uses LAPACK `syev` and supports the optional `descending` argument.

The regression test passed for:

- upstream canonical unweighted HJ bending and published diagnostics;
- upstream canonical reciprocal-weighted HJ bending, including the published
  428-iteration result and weighted diagnostics;
- zero-weight fixed elements;
- equivalence of reciprocal weights supplied directly vs. `reciprocal=.true.`;
- weighted LRS bending of a correlation matrix;
- already-positive-definite and invalid-method error paths.

The example program also compiled and ran successfully. Its unweighted canonical
matrix matched the upstream README values to the tested tolerances.

A final source audit checked free-form source, maximum line length 132, absence
of semicolon-separated code statements, absence of legacy real-kind syntax, and
explicit `INTENT` plus trailing FORD `!!` documentation on every dummy argument.
