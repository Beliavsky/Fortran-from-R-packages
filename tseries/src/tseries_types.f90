! Part of the experimental modern Fortran translation of tseries 0.10-62.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original tseries authors retain copyright; see NOTICE and ORIGIN.md.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only

module tseries_types
   use tseries_kinds, only : dp
   implicit none
   private

   public :: test_result
   public :: arma_result
   public :: garch_result
   public :: drawdown_result
   public :: portfolio_result
   public :: bds_result

   type :: test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp), allocatable :: parameters(:)
      integer :: status = 0
      character(len=:), allocatable :: method
      character(len=:), allocatable :: message
   end type test_result

   type :: bds_result
      real(dp), allocatable :: statistic(:, :)
      real(dp), allocatable :: p_value(:, :)
      real(dp), allocatable :: eps(:)
      integer :: max_embedding = 0
      integer :: status = 0
      character(len=:), allocatable :: message
   end type bds_result

   type :: arma_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: css = huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: p = 0
      integer :: q = 0
      integer :: iterations = 0
      integer :: status = 0
      logical :: include_intercept = .true.
      character(len=:), allocatable :: message
   end type arma_result

   type :: garch_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: conditional_variance(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp) :: negative_log_likelihood = huge(1.0_dp)
      integer :: p = 0
      integer :: q = 0
      integer :: iterations = 0
      integer :: status = 0
      character(len=:), allocatable :: message
   end type garch_result

   type :: drawdown_result
      real(dp) :: maximum = 0.0_dp
      integer :: from_index = 0
      integer :: to_index = 0
   end type drawdown_result

   type :: portfolio_result
      real(dp), allocatable :: weights(:)
      real(dp) :: expected_return = 0.0_dp
      real(dp) :: standard_deviation = 0.0_dp
      integer :: status = 0
      character(len=:), allocatable :: message
   end type portfolio_result

end module tseries_types
