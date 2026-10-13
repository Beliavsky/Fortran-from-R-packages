! SPDX-License-Identifier: GPL-2.0-only
! Derived from matrixTests 0.2.3.1 by Karolis Koncevicius.
module matrixtests_types
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   implicit none
   private
   integer, parameter, public :: mt_ok = 0, mt_invalid = 1, mt_insufficient = 2
   integer, parameter, public :: mt_degenerate = 3, mt_limit = 4, mt_nonfinite = 5
   type, public :: test_result
      !! Union of the numerical columns returned by the upstream test families.
      integer :: obs = 0, obs_x = 0, obs_y = 0, obs_tot = 0, obs_paired = 0, obs_groups = 0
      integer :: status = mt_invalid
      real(dp) :: statistic, pvalue, conf_low, conf_high, df, df_num, df_denom
      real(dp) :: mean, mean_x, mean_y, mean_diff, var, var_x, var_y, var_diff, var_pooled, var_ratio, stderr
      real(dp) :: cor, skewness, kurtosis, sumsq_between, sumsq_within, meansq_between, meansq_within
      real(dp) :: mesor, amplitude, acrophase, rsquared, period, null, conf_level
      character(9) :: alternative = 'two.sided'
      logical :: exact = .false., corrected = .false.
   end type test_result
   public :: dp, empty_result
contains
   pure function empty_result() result(r)
      type(test_result) :: r
      real(dp) :: na
      na = ieee_value(0.0_dp, ieee_quiet_nan)
      r%statistic = na
      r%pvalue = na
      r%conf_low = na
      r%conf_high = na
      r%df = na
      r%df_num = na
      r%df_denom = na
      r%mean = na
      r%mean_x = na
      r%mean_y = na
      r%mean_diff = na
      r%var = na
      r%var_x = na
      r%var_y = na
      r%var_diff = na
      r%var_pooled = na
      r%var_ratio = na
      r%stderr = na
      r%cor = na
      r%skewness = na
      r%kurtosis = na
      r%sumsq_between = na
      r%sumsq_within = na
      r%meansq_between = na
      r%meansq_within = na
      r%mesor = na
      r%amplitude = na
      r%acrophase = na
      r%rsquared = na
      r%period = na
      r%null = na
      r%conf_level = na
   end function empty_result
end module matrixtests_types
