! Public-domain numerical translations; see LICENSE and NOTICE.md.
module burstfin
   use r_kinds, only: dp
   use burstfin_portfolio, only: alpha_proxy, alpha_proxy_grid, alpha_surface, slide_weight
   use burstfin_portfolio, only: var_add_benchmark, var_relative_benchmark
   use burstfin_factor, only: statistical_factor_model, factor_model_stat, factor_covariance
   use burstfin_shrink, only: var_shrink_eqcor
   implicit none
   public
end module burstfin
