! SPDX-License-Identifier: LGPL-3.0-or-later
! Based on MarkowitzR, copyright 2014-2020 Steven E. Pav.
module markowitzr_types
   use markowitzr_kinds, only: dp
   implicit none
   private

   type, public :: theta_result
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: covariance(:, :)
      integer :: n = 0
      integer :: pp = 0
      integer :: status = 0
      character(len=160) :: message = ''
   end type theta_result

   type, public :: markowitz_result
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: w(:, :)
      real(dp), allocatable :: w_covariance(:, :)
      integer, allocatable :: w_indices(:)
      integer :: n = 0
      integer :: p = 0
      integer :: ff = 0
      integer :: status = 0
      character(len=160) :: message = ''
   end type markowitz_result

end module markowitzr_types
