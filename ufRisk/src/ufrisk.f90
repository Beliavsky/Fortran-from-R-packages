! SPDX-License-Identifier: GPL-3.0-only
module ufrisk
   use kind_mod, only : dp
   use ufrisk_types
   use ufrisk_varcast, only : varcast
   use ufrisk_backtests, only : trafftest, covtest, lossfunc
   use ufrisk_loggarch, only : arfilt_coefficients, estimate_student_df
   use ufrisk_smoothing, only : long_memory_smooth
   implicit none
   public
end module ufrisk
