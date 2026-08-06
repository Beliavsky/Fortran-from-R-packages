# Porting notes

## Representation

The R package stores one sparse cash-flow matrix and one sparse price matrix
for each quotation date. The Fortran port uses a `bond_panel_t` containing
parallel arrays sorted by `(day, id, tupq)`. This retains only the nonzero
entries and avoids a dependency on a sparse-matrix library.

The price must be repeated on every cash-flow row for a given bond and day,
matching the matrices produced by the R routine `get_cfp_slist`.

## Estimator fidelity

`calc_dbar` and `calc_hhat_numerator` retain the original kernel factors
exactly. In particular, the cross-product term contains both maturity-kernel
weights, whereas the `dbar` denominator contains one. The kernels have a peak
of 0.75 and are not normalized in these routines. Consequently, even a small
exact-grid coupon example need not equal an ordinary unweighted least-squares
solution. The tests include an independently calculated reference for this
behaviour.

The interpolation weights used by `estimate_yield` preserve the orientation in
the R implementation. They are not silently changed to a more conventional
linear-interpolation convention.

## Prediction

R's `stats::loess` is replaced by a self-contained local quadratic regression
using nearest-neighbour tricube weights. It follows the important numerical
choices of the R path (quadratic local fit and smoothing log discounts), but it
is not intended to be bit-for-bit identical to `stats::loess`, whose fitting
and edge-handling implementation is R-specific.

The Fortran interpolation routines clamp requests outside the fitted grid to
the nearest boundary, corresponding to `rule = 2` in the R code.

## Simulation

The original `ycevo_data` routine relies heavily on R date classes and
`lubridate`. `simulate_bond_panel` preserves its computational purpose: it
creates coupon bonds, semiannual payments, a Nelson-Siegel/cubic evolving true
curve, discounted prices, and repeated price/cash-flow rows. It uses integer
quotation days rather than reproducing Gregorian business-day and leap-day
semantics.

For exact real-data work, construct `bond_panel_t` directly or use the CSV
reader.

## Omitted infrastructure

The following components are not numerical algorithms and are intentionally
not translated:

- `ggplot2` and base-R plotting.
- Tibbles, tidy selection, unnesting, and S3 class methods.
- `future.apply` parallel dispatch and `progressr` progress reporting.
- R assertion-message and column-renaming helpers.

The computational equivalents are exposed through explicit Fortran arguments,
derived types, status codes, and optional error messages.

## Licensing

The source package declares GPL-3. The Fortran translation is distributed
under GPL-3. The original computational source files are retained under
`original/` and remain under their original copyright.
