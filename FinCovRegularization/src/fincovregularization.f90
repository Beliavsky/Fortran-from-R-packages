! SPDX-License-Identifier: GPL-2.0-only
module fincovregularization
   use fincov_kinds, only : dp
   use fincov_status
   use fincov_types, only : cv_result
   use fincov_norms, only : f_norm2, o_norm2
   use fincov_regularization, only : banding, tapering, hard_thresholding, soft_thresholding, &
      ind_cov, threshold_min
   use fincov_factor_models, only : macro_factor_cov, fundamental_factor_cov, stat_factor_cov
   use fincov_portfolio, only : gmvp, risk_parity, risk_parity_objective
   use fincov_cv, only : banding_cv, tapering_cv, threshold_cv
   implicit none
   public
end module fincovregularization
