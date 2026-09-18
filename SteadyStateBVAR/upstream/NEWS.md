# SteadyStateBVAR 0.2.0

* `IRF()` argument `shock` renamed to `impulse`, for better clarity.
  If you have `shock = ` in existing calls to `IRF()`, rename it to `impulse = `.

* `forecast()`, `conditional_forecast()`, and `IRF()` now transform each 
  posterior draw before summarizing, rather than summarizing and then transforming.
  This fixes incorrect prediction/credible intervals and median point estimates/predictions
  when using `growth_rate_idx` argument. Mean point estimates/predictions are unaffected.
  
* Now `d_pred` is automatically created in `fit()` to make the
  user-experience smoother. However, it is still left as an argument.

* Fixed `IRF()` plotting: specifying only one of `response`/`impulse` 
  (leaving the other `NULL`) previously ignored the specified index and 
  plotted all response-impulse combinations. Now the specified index is 
  correctly held fixed while the other dimension varies.

* New `steady_state_priors_plot()` function for visualizing steady-state priors.

* `forecast()` now overlays posterior steady-state estimates on forecast 
  plots.

* Stan code has been optimized. This includes assuming that all prior covariance
  matrices are diagonal.

* Added a vignette `vignette("SteadyStateBVAR-intro")`

* Function documentations, vignettes, and the README have been revamped quite
  extensively.

# SteadyStateBVAR 0.1.1

* Stan code now uses the new array syntax (@andrjohns, #4).

* Small rework of plotting in `forecast()` and corrected documentation. 

# SteadyStateBVAR 0.1.0

* Initial CRAN release.

# SteadyStateBVAR 0.0.0.9000

* Initial development version.
* Preparing first CRAN submission.
