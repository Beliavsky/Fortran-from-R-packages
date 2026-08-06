! SPDX-License-Identifier: MIT
! Copyright (c) 2021 treasuryTR authors

module treasurytr
  use treasurytr_kinds, only : dp
  use treasurytr_math, only : convexity, mod_duration, period_total_return
  use treasurytr_series, only : carry_forward, percent_to_decimal, prepare_yields, &
    total_return, tt_success, tt_err_size, tt_err_maturity, tt_err_scale
  implicit none
  public
end module treasurytr
