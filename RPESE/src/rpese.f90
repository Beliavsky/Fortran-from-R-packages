! SPDX-License-Identifier: GPL-3.0-or-later
module rpese
  use rpese_kinds, only : dp, pi
  use rpese_types
  use rpese_measures
  use rpese_timeseries, only : fit_periodogram, lag1_correlation, fit_ar1, polynomial_design
  use rpese_core, only : estimate_se, estimate_se_matrix, influence_values, spectral_variance
  use rpese_api
  use rpese_compat
  implicit none
  public
end module rpese
