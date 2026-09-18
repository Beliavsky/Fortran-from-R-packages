# NOTICE and provenance

## Upstream project

This work is a modern Fortran translation of **poLCA 1.6.0.2**, "Polytomous Variable Latent Class Analysis," authored by Drew A. Linzer and Jeffrey Lewis. The attached upstream `DESCRIPTION` identifies the package license as `GPL (>= 2)`.

Upstream project and documentation: https://github.com/dlinzer/poLCA

The numerical design follows the upstream R and C sources preserved under `upstream/`, including the EM response-probability update, posterior class calculation, multinomial-logit prior update, observed-likelihood Newton derivatives, and score-based standard-error construction.

## License

The upstream package is GPL version 2 or later. The translated Fortran source is distributed under the same terms. `LICENSE` contains the GNU General Public License version 2 text; the "or later" option follows the upstream package declaration.

## Dependency/provenance decisions

The upstream R package declares `MASS` and `scatterplot3d`. `scatterplot3d` is presentation-only for the translated scope and is omitted. Upstream uses `MASS::ginv` for an information-matrix generalized inverse. Rather than vendor MASS, BLAS, or LAPACK, this translation implements the required symmetric Moore-Penrose pseudoinverse with a local Jacobi eigendecomposition and validates the resulting standard errors independently.

No source from MASS, scatterplot3d, BLAS, LAPACK, `rfortran-compat`, or another translated R package is copied into this package. No system BLAS/LAPACK link is used.

A repository search was made for an existing compatible `poLCA` translation/shared dependency before implementation; no package-specific reusable dependency was required for this self-contained numerical core.

## Material differences from R

- Callers provide numeric matrices directly rather than R formulas, data frames, factors, or S3 objects.
- Plotting, printing, formatting, and interactive console behavior are omitted.
- Missing manifest values are represented by integer zero.
- Random-number streams are not intended to match R bit-for-bit.
- `polca_simdata` uses independent cell-wise missingness when a missing fraction is requested.
- Stable log-sum-exp calculations are used for posterior normalization instead of the upstream C routine's `DBL_MAX` likelihood scaling trick; these are mathematically equivalent apart from floating-point edge behavior.
