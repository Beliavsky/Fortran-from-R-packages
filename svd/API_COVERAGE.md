# API coverage

Coverage basis: distinct exported computational R functions plus explicitly registered
computational S3 methods in `upstream/NAMESPACE`. Printing, plotting, presentation, R type
coercion-only behavior, and S4 container-dispatch conveniences are excluded from the
denominator under the task's instruction to skip R-specific interface code.

Mapped: 12 / 12 exported computational functions.

The matrix-free API is represented by `extmat_operator`, which extends the shared
`RSpectra` `linear_operator` type. `make_extmat` stores caller-supplied forward/transpose
callbacks together with persistent polymorphic context state, replacing R closures and
environments without requiring C interoperability.

Real leading SVD/eigen computations reuse sibling `RSpectra`. Complex dense SVD reuses
sibling `rfortran-linalg`. The bundled upstream PROPACK and nuTRLan implementations are
not translated or vendored because compatible shared top-level numerical packages already
exist in the target repository.

See the exact per-function mapping records in `fpm.toml` and the matching table in
`README.md`.
