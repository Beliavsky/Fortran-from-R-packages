! SPDX-License-Identifier: GPL-2.0-or-later
module yieldcurve
   use yieldcurve_kinds, only : dp
   use yieldcurve_status, only : yc_success, yc_invalid_argument, yc_dimension_error, &
      yc_rank_deficient, yc_no_solution
   use yieldcurve_factors, only : beta1_spot, beta2_spot, beta1_forward, beta2_forward, &
      factor_beta1, factor_beta2
   use yieldcurve_models, only : nelson_siegel_fit, svensson_fit, ns_rates, svensson_rates
   implicit none
   public
end module yieldcurve
