# yieldcurves 0.1.0

* Initial release.
* Nelson-Siegel (1987) and Svensson (1994) yield curve fitting with
  multi-start optimization and optional observation weights.
* Cubic spline interpolation via `stats::splinefun`.
* Forward rate extraction (analytical for NS/Svensson, numerical for splines).
* Discount factor computation with continuous, annual, and semi-annual
  compounding.
* Duration and convexity for zero-coupon bonds via `yc_duration()`.
* Coupon bond duration, modified duration, and convexity via
  `yc_bond_duration()` (annual, semi-annual, or continuous compounding).
* Z-spread computation via `yc_zspread()`.
* Key rate durations via `yc_key_rate_duration()`.
* Par-to-zero and zero-to-par rate conversions via bootstrap stripping,
  supporting annual and semi-annual coupon frequencies.
* Principal component decomposition of yield curve time series.
* Carry and roll-down analysis.
* Slope measures (2s10s, 2s30s, butterfly) and level-slope-curvature
  decomposition.
