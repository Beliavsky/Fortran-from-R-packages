! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on SharpeR, copyright 2012-2025 Steven E. Pav.
module sharper_types
   use sharper_kinds, only: dp
   implicit none
   private

   type, public :: sr_result
      real(dp), allocatable :: value(:)
      integer, allocatable :: df(:)
      real(dp) :: c0 = 0.0_dp
      real(dp) :: ope = 1.0_dp
      real(dp), allocatable :: rescal(:)
      real(dp), allocatable :: cumulants(:, :)
   end type sr_result

   type, public :: sropt_result
      real(dp) :: value = 0.0_dp
      integer :: df1 = 0
      integer :: df2 = 0
      real(dp) :: drag = 0.0_dp
      real(dp) :: ope = 1.0_dp
      real(dp) :: t2 = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: covariance(:, :)
   end type sropt_result

   type, public :: del_sropt_result
      real(dp) :: value = 0.0_dp
      real(dp) :: sub_value = 0.0_dp
      real(dp) :: delta_value = 0.0_dp
      integer :: df1 = 0
      integer :: df2 = 0
      integer :: df1_sub = 0
      real(dp) :: drag = 0.0_dp
      real(dp) :: ope = 1.0_dp
      real(dp) :: t2 = 0.0_dp
      real(dp) :: t2_sub = 0.0_dp
      real(dp) :: t2_delta = 0.0_dp
   end type del_sropt_result

   type, public :: test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: estimate = 0.0_dp
      real(dp) :: null_value = 0.0_dp
      real(dp) :: conf_low = 0.0_dp
      real(dp) :: conf_high = 0.0_dp
      integer :: df1 = 0
      integer :: df2 = 0
      integer :: parameter = 0
      character(len=24) :: alternative = 'two.sided'
      character(len=72) :: method = ''
      integer :: status = 0
   end type test_result

   type, public :: power_result
      real(dp) :: n = 0.0_dp
      real(dp) :: n_epoch = 0.0_dp
      real(dp) :: effect = 0.0_dp
      real(dp) :: sig_level = 0.05_dp
      real(dp) :: power = 0.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      integer :: status = 0
   end type power_result

   type, public :: moment_vcov_result
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: covariance(:, :)
      integer :: n = 0
      integer :: p = 0
      integer :: status = 0
   end type moment_vcov_result

end module sharper_types
