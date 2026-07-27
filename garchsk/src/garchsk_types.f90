! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from GARCHSK 0.1.0, Copyright (C) 2021 Kei Nakagawa.
module garchsk_types
   use garchsk_kinds, only : dp
   implicit none
   private

   type, public :: moment_path
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: skewness(:)
      real(dp), allocatable :: kurtosis(:)
      logical :: success = .false.
      character(len=160) :: message = ''
   end type moment_path

   type, public :: forecast_result
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: skewness(:)
      real(dp), allocatable :: kurtosis(:)
      logical :: success = .false.
      character(len=160) :: message = ''
   end type forecast_result

   type, public :: estimate_result
      real(dp), allocatable :: params(:)
      real(dp), allocatable :: standard_errors(:)
      real(dp), allocatable :: upstream_standard_errors(:)
      real(dp), allocatable :: t_statistics(:)
      real(dp) :: negative_log_likelihood = huge(1.0_dp)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp) :: upstream_aic = huge(1.0_dp)
      real(dp) :: upstream_bic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      logical :: converged = .false.
      logical :: covariance_available = .false.
      character(len=160) :: message = ''
   end type estimate_result

end module garchsk_types
