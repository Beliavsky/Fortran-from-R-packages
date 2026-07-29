# Porting notes

## Source package

- Package: `mvtnorm`
- Version: 1.4-2
- License: GPL version 2
- Supplied source date: July 2026

The upstream package combines R interfaces, C matrix and likelihood kernels,
Alan Genz's Fortran probability code, TVPACK, and Tetsuhisa Miwa's C algorithm.
The port consolidates the numerical API into modern Fortran modules.

## Probability engine

The default implementation follows the conditional-normal transformation used
by Genz-style multivariate probability algorithms:

1. Standardize covariance matrices to correlations.
2. Reorder variables by marginal interval probability.
3. Form a lower Cholesky factor.
4. Integrate sequential conditional probability widths over a unit cube.
5. Use shifted Halton points, antithetic pairing, and independent randomized
   batches for an error estimate.

For Student-t laws, one additional unit-cube coordinate generates the
chi-square mixing variable. Shifted and Kshirsagar noncentrality conventions
are kept distinct.

This is a modern reimplementation, not a mechanical conversion of the
fixed-form `mvt.f` common-block and entry-point structure. It eliminates global
common blocks, fixed dimension-1000 work arrays, and implicit interfaces.

## TVPACK and Miwa selectors

The original package offers three user-visible method selectors. The Fortran
port preserves those selectors, but only the bivariate normal TVPACK-compatible
route is a dedicated deterministic calculation. It uses the Plackett
correlation integral with adaptive Simpson subdivision.

Higher-dimensional TVPACK and Miwa requests use more conservative settings of
the common conditional-integration engine. This avoids embedding historical
fixed-form/common-block code or a large C translation, but means results and
evaluation counts are not bit-for-bit identical to R `mvtnorm`.

## Matrix classes

R's `ltMatrices` and `syMatrices` hold many matrices plus metadata in one
object. Fortran callers use ordinary rank-two matrices and, when needed,
rank-three batches. Lower-triangular packing is column-major by default, with a
row-major option matching the upstream constructor.

The port explicitly symmetrizes covariance, precision, correlation, and partial
correlation outputs to prevent one-ulp cross-platform asymmetry.

## Likelihoods and scores

Exact normal log densities use triangular solves. Interval probabilities use
the same rectangle-probability engine as `pmvnorm`. Mixed exact/interval data
use a marginal density plus a conditional rectangle probability.

The upstream package contains specialized analytical and simulation-based score
kernels. The Fortran port returns deterministic central-difference scores for
means and packed lower-Cholesky parameters. This broadens consistency and
reduces implementation risk, at the cost of speed.

`deperma_score` and `destandardize_score` similarly use numerical Jacobians.

## Random numbers

The port uses Fortran's intrinsic pseudorandom generator with deterministic
integer seeding and Box-Muller normal generation. Seeded values are reproducible
within a compiler/runtime family but are not expected to match R's RNG stream.

## Portability safeguards

- No expression relies on short-circuit evaluation.
- Matrix symmetry is enforced explicitly.
- Probability points are clipped away from exact zero and one before inverse
  CDF evaluation.
- Singular probability correlation matrices produce a typed failure rather
  than an out-of-bounds or divide-by-zero path.
- Quantile inversion uses bracketing and bisection rather than a noisy R GLM
  helper.
- Source is ASCII-only and stays within the standard 132-column free-form
  limit.
