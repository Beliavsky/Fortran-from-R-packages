! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from OptionPricing 0.1.2 by Wolfgang Hormann and Kemal Dingec.
module optionpricing
   use optionpricing_kinds, only : dp, pi
   use optionpricing_types, only : european_result, greeks_result, moments_result, conditional_result
   use optionpricing_european, only : bs_european_call, bs_european_put, bs_ec, bs_ep
   use optionpricing_asian_analytic, only : bs_a, covariance_conditional_log_prices, &
      conditional_average_moments, eval_ecv_a, asian_call_app_lord, eval_ecv, &
      find_bcv, eval_lb, eval_equad, eval_eqcv
   use optionpricing_asian_mc, only : asian_call_naive_mc, asian_call_ncv_lr_mc, &
      asian_call_cmc_cv, asian_call_best_mc, simulate_asian_call_z, &
      asian_call_naive_greeks_z, conditional_estimates_z, build_conditional_samples
   use optionpricing_asian_qmc, only : korobov_lattice, randomized_korobov_normals, &
      naive_pca_matrix, conditional_generation_matrix, asian_call_naive_qmc, &
      asian_call_best_qmc, asian_call
   implicit none
   public
end module optionpricing
