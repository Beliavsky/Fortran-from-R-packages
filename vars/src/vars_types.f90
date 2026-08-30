! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_types
   use r_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: vars_success = 0
   integer, parameter, public :: vars_invalid_argument = -1
   integer, parameter, public :: vars_singular = -2
   integer, parameter, public :: vars_not_identified = -3
   integer, parameter, public :: vars_no_convergence = 1

   integer, parameter, public :: var_none = 0
   integer, parameter, public :: var_const = 1
   integer, parameter, public :: var_trend = 2
   integer, parameter, public :: var_both = 3

   type, public :: var_model
      integer :: p = 0
      integer :: k = 0
      integer :: nobs = 0
      integer :: totobs = 0
      integer :: nreg = 0
      integer :: deterministic = var_const
      integer :: season = 0
      integer :: exogen_cols = 0
      integer, allocatable :: df_resid(:)
      real(dp), allocatable :: y(:, :)
      real(dp), allocatable :: response(:, :)
      real(dp), allocatable :: x(:, :)
      real(dp), allocatable :: coef(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: resid(:, :)
      real(dp), allocatable :: sigma_u(:, :)
      real(dp), allocatable :: a(:, :, :)
      logical, allocatable :: active(:, :)
   end type var_model

   type, public :: var_selection_result
      real(dp), allocatable :: criteria(:, :)
      integer :: selection(4) = 0
   end type var_selection_result

   type, public :: vars_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
      real(dp) :: p_value = 1.0_dp
   end type vars_test_result

   type, public :: forecast_result
      real(dp), allocatable :: point(:, :)
      real(dp), allocatable :: lower(:, :)
      real(dp), allocatable :: upper(:, :)
      real(dp), allocatable :: se(:, :)
   end type forecast_result

   type, public :: svar_result
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: b(:, :)
      real(dp), allocatable :: a_se(:, :)
      real(dp), allocatable :: b_se(:, :)
      real(dp), allocatable :: sigma_u(:, :)
      real(dp), allocatable :: gamma(:)
      integer :: iterations = 0
      logical :: converged = .false.
      real(dp) :: convergence = huge(1.0_dp)
      real(dp) :: lr_statistic = 0.0_dp
      real(dp) :: lr_df = 0.0_dp
      real(dp) :: lr_p_value = 1.0_dp
   end type svar_result

   type, public :: svec_result
      real(dp), allocatable :: sr(:, :)
      real(dp), allocatable :: lr(:, :)
      real(dp), allocatable :: sigma_u(:, :)
      real(dp), allocatable :: gamma(:)
      integer :: restrictions_lr = 0
      integer :: restrictions_sr = 0
      integer :: iterations = 0
      logical :: converged = .false.
      real(dp) :: convergence = huge(1.0_dp)
      real(dp) :: lr_statistic = 0.0_dp
      real(dp) :: lr_df = 0.0_dp
      real(dp) :: lr_p_value = 1.0_dp
   end type svec_result

end module vars_types
