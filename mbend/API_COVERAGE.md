# API coverage

## Exported computational API

Upstream `NAMESPACE` exports one computational function: `bend`.

### `bend` -> `mbend_api:bend`

Translated numerical behavior:

- validates square symmetric input and compatible nonnegative symmetric weights;
- unweighted and weighted Jorjani et al. (`method="hj"`) bending;
- unweighted and weighted Schaeffer (`method="lrs"`) bending;
- reciprocal nonzero weights;
- zero weights that keep selected matrix elements fixed;
- correlation-matrix detection and fixed unit diagonal;
- upstream default iteration limit and HJ eigenvalue floor;
- initial/final eigenspectra and convergence iteration count;
- minimum/maximum deviation and their one-based locations;
- average deviation, AAD, Pearson correlation, and RMSD;
- weighted positive-count, weighted AAD, weighted correlation, and weighted RMSD.

Material interface differences:

- R row names and column names are not represented.
- R messages and warnings are replaced by status/result fields.
- R's named list return value is represented by the `bend_result` derived type.
- An effective all-zero weight matrix receives an explicit error instead of
  propagating NaNs from division by zero.

Coverage status is **substantial**, not complete, because the computational
operation is mapped while R-specific interface semantics are intentionally not.
