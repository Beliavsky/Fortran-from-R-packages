! SPDX-License-Identifier: GPL-2.0-or-later
module flexsurv
  use flexsurv_kinds, only : dp
  use flexsurv_distributions
  use flexsurv_spline
  use flexsurv_splines2ns
  use flexsurv_fit
  use flexsurv_spline_fit
  use flexsurv_spline_interactions
  use flexsurv_mixture
  use flexsurv_mixture_full
  use flexsurv_fmixmsm
  use flexsurv_multistate
  use flexsurv_multistate_uncertainty
  use flexsurv_shared_multistate
  use flexsurv_final_states
  use flexsurv_standardize
  use flexsurv_standardize_advanced
  use flexsurv_ajfit
  use flexsurv_rtrunc
  use flexsurv_fracpoly
  use flexsurv_diagnostics
  use flexsurv_custom
  use relsurv_ratetable, only : ratetable_type, make_ratetable
  implicit none
  public
end module flexsurv
