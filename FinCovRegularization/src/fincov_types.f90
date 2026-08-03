! SPDX-License-Identifier: GPL-2.0-only
module fincov_types
   use fincov_kinds, only : dp
   use fincov_status, only : fincov_ok
   implicit none
   private

   type, public :: cv_result
      character(len=:), allocatable :: regularization
      character(len=:), allocatable :: method
      real(dp) :: parameter_opt = 0.0_dp
      real(dp), allocatable :: cv_error(:)
      real(dp), allocatable :: parameter_grid(:)
      integer :: parameter_index = 0
      integer :: n_cv = 0
      integer :: seed = 0
      character(len=1) :: norm = 'F'
      real(dp) :: h = 0.5_dp
      integer :: status = fincov_ok
   end type cv_result
end module fincov_types
