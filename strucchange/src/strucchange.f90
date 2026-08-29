! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange
   use r_kinds, only : dp
   use strucchange_breakpoints, only : breakpoint_confint_result
   use strucchange_breakpoints, only : breakpoint_path_result
   use strucchange_breakpoints, only : breakpoint_result
   use strucchange_breakpoints, only : best_break_count
   use strucchange_breakpoints, only : breakpoint_confidence_intervals
   use strucchange_breakpoints, only : compute_breakpoint_path
   use strucchange_breakpoints, only : compute_breakpoints
   use strucchange_breakpoints, only : segmented_fit
   use strucchange_efp, only : moving_estimates_process
   use strucchange_efp, only : ols_cusum, ols_mosum
   use strucchange_efp, only : process_max, process_max_l2
   use strucchange_efp, only : process_mean_l2, process_range
   use strucchange_efp, only : recursive_cusum
   use strucchange_efp, only : recursive_estimates_process
   use strucchange_efp, only : recursive_mosum
   use strucchange_efp, only : score_cusum, score_mosum
   use strucchange_fstats, only : ave_f_statistic, compute_fstats
   use strucchange_fstats, only : exp_f_statistic, fstats_result
   use strucchange_fstats, only : sup_f_statistic
   use strucchange_functionals, only : cat_l2_bb_critical_value
   use strucchange_functionals, only : cat_l2_bb_pvalue, cat_l2_bb_statistic
   use strucchange_functionals, only : max_mosum_critical_value
   use strucchange_functionals, only : max_mosum_pvalue, max_mosum_statistic
   use strucchange_functionals, only : sup_lm_critical_value
   use strucchange_functionals, only : sup_lm_pvalue, sup_lm_statistic
   use strucchange_gefp, only : generalized_fluctuation_process
   use strucchange_monitoring, only : log_plus
   use strucchange_monitoring, only : monitor_me_critical_value
   use strucchange_monitoring, only : monitor_ols_cusum_boundary
   use strucchange_monitoring, only : monitor_power_boundary
   use strucchange_monitoring, only : monitor_re_boundary
   use strucchange_monitoring, only : monitor_re_critical_value
   use strucchange_monitoring, only : mre_critical_value, pargmax_v
   use strucchange_pvalues, only : efp_pvalue, fstats_pvalue
   use strucchange_recresid, only : recursive_residuals
   use strucchange_regression, only : inverse_crossprod, ols_fit, root_matrix
   implicit none
   private
   public :: dp
   public :: best_break_count
   public :: breakpoint_confidence_intervals
   public :: breakpoint_confint_result
   public :: cat_l2_bb_critical_value
   public :: cat_l2_bb_pvalue
   public :: cat_l2_bb_statistic
   public :: ave_f_statistic
   public :: breakpoint_path_result
   public :: breakpoint_result
   public :: compute_breakpoint_path
   public :: compute_breakpoints
   public :: compute_fstats
   public :: efp_pvalue
   public :: exp_f_statistic
   public :: fstats_pvalue
   public :: fstats_result
   public :: generalized_fluctuation_process
   public :: inverse_crossprod
   public :: log_plus
   public :: max_mosum_critical_value
   public :: max_mosum_pvalue
   public :: max_mosum_statistic
   public :: monitor_me_critical_value
   public :: monitor_ols_cusum_boundary
   public :: monitor_power_boundary
   public :: monitor_re_boundary
   public :: monitor_re_critical_value
   public :: moving_estimates_process
   public :: mre_critical_value
   public :: ols_cusum
   public :: ols_fit
   public :: ols_mosum
   public :: pargmax_v
   public :: process_max
   public :: process_max_l2
   public :: process_mean_l2
   public :: process_range
   public :: recursive_cusum
   public :: recursive_estimates_process
   public :: recursive_mosum
   public :: recursive_residuals
   public :: root_matrix
   public :: score_cusum
   public :: score_mosum
   public :: segmented_fit
   public :: sup_lm_critical_value
   public :: sup_lm_pvalue
   public :: sup_lm_statistic
   public :: sup_f_statistic
end module strucchange
