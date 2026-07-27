# Porting notes

## Array orientation

Return observations are rows and funds are columns. Factor observations are rows
and factor series are columns. Pairwise screening arrays are indexed as
`(coefficient, focal fund, peer fund)`.

## Typed control and results

R control lists are replaced by `peer_control`. Important fields include:

- `test_type = 1` for asymptotic inference or `2` for bootstrap inference
- `ttype = 1` for the direct ratio difference or `2` for the studentized
  cross-product form used by the upstream default
- `hac` for long-run covariance estimation
- `n_boot`, `block_length`, and `p_boot` for bootstrap inference
- `has_lambda` and `lambda` for fixed-lambda peer ratios
- `gamma_pos` and `gamma_neg` for the positive/negative split
- `screen_beta` for factor-coefficient screening

All procedures return a status code and diagnostic message rather than raising an
R exception.

## HAC implementation

The R package delegates HAC covariance to `sandwich::vcovHAC`. The Fortran port
uses a self-contained automatic Parzen-kernel estimator. The same statistical
purpose is preserved, but exact numerical equality with R's dependency stack is
not promised.

## Bootstrap implementation

The IID and circular block algorithms follow the upstream methodology. Fortran's
intrinsic generator and explicit integer/geometric samplers replace R's RNG, so
identical seeds do not imply identical R and Fortran resamples. Results are fully
reproducible within the Fortran implementation.

The block-size selectors fit a bivariate VAR(1), resample its residuals with a
stationary bootstrap, simulate under the fitted model, and select the candidate
whose rejection frequency is closest to the requested size.

## Peer-ratio adjustment

`adjust_pi` numerically inverts the upstream finite-sample forward adjustment.
The default uses 50 bisection iterations. The `fast` control path uses 15.
Data-driven lambda selection evaluates the upstream grid
`[0.3, 0.4, 0.5, 0.6, 0.7]` with bootstrap MSE selection.

## Missing values and degenerate pairs

IEEE NaNs are removed pairwise. A pair is omitted when it has too few complete
observations, a singular regression, zero variance, or another non-finite
statistic. Peer counts are therefore based on valid p-values, not simply the
number of columns.

## Deliberate safety changes

- Logical expressions do not depend on short-circuit evaluation.
- Reused allocatable result objects are explicitly initialized, avoiding
  optimizer-dependent uninitialized-component diagnostics.
- Invalid block lengths, probability levels, control values, and dimensions are
  rejected before numerical work.
- Bootstrap screening requires an explicit positive block length; expensive
  data-driven block selection is exposed as a separate routine.

## Source naming

The top-level `peerperformance` module re-exports the public API. Internal module
names are prefixed `peerperformance_` to avoid namespace collisions.
