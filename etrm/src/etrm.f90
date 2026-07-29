! SPDX-License-Identifier: MIT
! Derived from etrm 1.0.2, Copyright (c) 2021 etrm authors.
module etrm
   use etrm_kinds, only : dp
   use etrm_status, only : etrm_ok, etrm_err_size, etrm_err_argument, &
      etrm_err_allocation, etrm_err_linear_solve
   use etrm_dates, only : civil_to_day, day_offset
   use etrm_types, only : strategy_result, strategy_summary, msfc_result, summarize_strategy
   use etrm_strategies, only : cppi, dppi, obpi, shpi, slpi
   use etrm_msfc, only : msfc, maximum_smoothness_forward_curve
   implicit none
   public
end module etrm
