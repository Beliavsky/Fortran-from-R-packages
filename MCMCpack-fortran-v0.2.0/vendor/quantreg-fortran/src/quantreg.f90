! SPDX-License-Identifier: GPL-2.0-or-later
module quantreg
  use quantreg_kinds, only : dp
  use quantreg_types, only : rq_result, rq_multi_result, nlrq_result, lprq_result
  use quantreg_dense, only : rq_fit_fnb, rq_fit_fnc, rq_fit_qfnb, rq_wfit_fnb
  use quantreg_dense, only : rq_fit_lasso, rq_fit_scad, check_loss, check_loss_sum
  use quantreg_select, only : qselect, kuantiles
  use quantreg_local, only : lprq
  use quantreg_nonlinear, only : nlrq_fit, nlrq_model_fn
  use quantreg_bootstrap, only : rq_bootstrap_xy, seed_rng
  use quantreg_utils, only : rq_fit_pfn, recursive_least_squares, combinations, random_exponential
  implicit none
  public
end module quantreg
