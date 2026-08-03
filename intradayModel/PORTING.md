# Porting notes

## Direct translations

The following numerical structure follows the R implementation:

- two-dimensional periodic state-space model on log volume;
- `[1, 1]` observation vector;
- daily-state transition and innovation only at day boundaries;
- dynamic-state AR(1) transition at every bin;
- Kalman filtering and Rauch-Tung-Striebel smoothing;
- closed-form EM updates for all eight parameter groups;
- centered seasonal profile update;
- ordinary and accelerated EM choices;
- one-bin-ahead forecast components and burn-in removal;
- MAE, MAPE, and RMSE output.

Default starting values match `specify_uniss` in the R source.

## Fortran interfaces

R lists and S3 objects are represented by derived types. Parameter presence is explicit
through `parameter_mask`; callers do not need to encode unspecified values as `NA`.
Input is an `n_bin x n_day` matrix rather than an `xts` object.

## Numerical safeguards

- Observation and state variances have a `1e-12` lower floor after EM updates.
- Near-singular 2-by-2 predicted covariance matrices receive a small diagonal ridge in
  the smoother inversion.
- Invalid accelerated candidates fall back to the second ordinary EM update.
- The ordinary EM path stores the newly computed converged iterate. The upstream R code
  breaks before assigning that final iterate to its internal object.
- Logical checks involving nonfinite data are written without relying on short-circuit
  evaluation, which Fortran does not guarantee.

These safeguards do not alter the model equations.

## Omitted infrastructure

- `generate_plots` and all `ggplot2`/`patchwork` code;
- `xts`/`zoo` calendar conversion and warning formatting;
- bundled serialized R datasets in the compiled library;
- R S3 printing and attribute conventions.

The original package tree is retained under `original/` for provenance.
