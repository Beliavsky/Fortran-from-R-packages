! SPDX-License-Identifier: GPL-2.0-only
! Derived from vrtest 1.2 by Jae H. Kim.
module vrtest_types
   use vrtest_kinds, only : dp
   implicit none
   private

   type, public :: ar1_result
      real(dp) :: coefficient = 0.0_dp
      real(dp) :: standard_error = 0.0_dp
      real(dp) :: innovation_variance = 0.0_dp
   end type ar1_result

   type, public :: lm_result
      real(dp) :: homoskedastic = 0.0_dp
      real(dp) :: heteroskedastic = 0.0_dp
      real(dp) :: variance_ratio = 0.0_dp
   end type lm_result

   type, public :: lmcd_result
      real(dp), allocatable :: homoskedastic(:)
      real(dp), allocatable :: heteroskedastic(:)
      real(dp) :: cd_homoskedastic = 0.0_dp
      real(dp) :: cd_heteroskedastic = 0.0_dp
   end type lmcd_result

   type, public :: auto_vr_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: variance_ratio_sum = 0.0_dp
      real(dp) :: bandwidth = 0.0_dp
      real(dp) :: ar1_coefficient = 0.0_dp
   end type auto_vr_result

   type, public :: auto_q_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: selected_lag = 1
   end type auto_q_result

   type, public :: bootstrap_result
      integer, allocatable :: holding_periods(:)
      real(dp), allocatable :: lm_p_values(:)
      real(dp) :: cd_p_value = 1.0_dp
      real(dp), allocatable :: confidence_intervals(:,:)
   end type bootstrap_result

   type, public :: auto_bootstrap_result
      real(dp) :: test_statistic = 0.0_dp
      real(dp) :: variance_ratio_sum = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp), allocatable :: statistic_interval(:)
      real(dp), allocatable :: vr_sum_interval(:)
   end type auto_bootstrap_result

   type, public :: chow_denning_result
      integer, allocatable :: holding_periods(:)
      real(dp) :: cd_homoskedastic = 0.0_dp
      real(dp) :: cd_heteroskedastic = 0.0_dp
      real(dp) :: critical_values(3) = 0.0_dp
   end type chow_denning_result

   type, public :: wald_result
      integer, allocatable :: holding_periods(:)
      real(dp) :: statistic = 0.0_dp
      real(dp) :: critical_values(3) = 0.0_dp
      integer :: solve_info = 0
   end type wald_result

   type, public :: wright_result
      real(dp), allocatable :: statistics(:,:)
   end type wright_result

   type, public :: joint_wright_result
      integer, allocatable :: holding_periods(:)
      real(dp) :: rank_uniform = 0.0_dp
      real(dp) :: rank_normal = 0.0_dp
      real(dp) :: sign = 0.0_dp
   end type joint_wright_result

   type, public :: wright_critical_result
      integer, allocatable :: holding_periods(:)
      real(dp), allocatable :: critical_values(:,:)
   end type wright_critical_result

   type, public :: subsample_result
      integer, allocatable :: holding_periods(:)
      integer, allocatable :: block_lengths(:)
      real(dp), allocatable :: p_values(:)
   end type subsample_result

   type, public :: panel_vr_result
      real(dp) :: max_absolute_statistic = 0.0_dp
      real(dp) :: sum_square_statistic = 0.0_dp
      real(dp) :: mean_statistic = 0.0_dp
      real(dp) :: max_absolute_p_value = 1.0_dp
      real(dp) :: sum_square_p_value = 1.0_dp
      real(dp) :: mean_p_value = 1.0_dp
   end type panel_vr_result

   type, public :: vr_minus_one_result
      real(dp) :: automatic = 0.0_dp
      integer, allocatable :: holding_periods(:)
      real(dp), allocatable :: values(:)
   end type vr_minus_one_result

   type, public :: vr_curve_result
      integer, allocatable :: holding_periods(:)
      real(dp), allocatable :: variance_ratios(:)
      real(dp), allocatable :: lower_95(:)
      real(dp), allocatable :: upper_95(:)
   end type vr_curve_result

   type, public :: average_exponential_result
      real(dp) :: exponential_lm = 0.0_dp
      real(dp) :: exponential_lr = 0.0_dp
   end type average_exponential_result

   type, public :: spectral_shape_result
      real(dp) :: anderson_darling = 0.0_dp
      real(dp) :: cramer_von_mises = 0.0_dp
      real(dp) :: maximum = 0.0_dp
   end type spectral_shape_result

   type, public :: generalized_spectral_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: bootstrap_critical_values(3) = 0.0_dp
   end type generalized_spectral_result

   type, public :: dl_result
      real(dp) :: cp_statistic = 0.0_dp
      real(dp) :: kp_statistic = 0.0_dp
      real(dp) :: cp_p_value = 1.0_dp
      real(dp) :: kp_p_value = 1.0_dp
   end type dl_result

   type, public :: chen_deo_result
      integer, allocatable :: holding_periods(:)
      real(dp) :: variance_ratio_sum = 0.0_dp
      real(dp) :: qp_statistic = 0.0_dp
      real(dp) :: chi_square_upper_quantiles(5) = 0.0_dp
      integer :: solve_info = 0
   end type chen_deo_result

end module vrtest_types
