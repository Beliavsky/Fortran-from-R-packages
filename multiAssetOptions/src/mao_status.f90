! This file is part of multiAssetOptions-fortran.
! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module mao_status
   implicit none
   private

   integer, parameter, public :: mao_success = 0
   integer, parameter, public :: mao_invalid_argument = 1
   integer, parameter, public :: mao_allocation_error = 2
   integer, parameter, public :: mao_solver_failure = 3
   integer, parameter, public :: mao_step_failure = 4

   type, public :: status_type
      integer :: code = mao_success
      character(len=:), allocatable :: message
   contains
      procedure :: ok => status_ok
   end type status_type

   public :: clear_status, set_status

contains

   logical function status_ok(self)
      class(status_type), intent(in) :: self
      status_ok = self%code == mao_success
   end function status_ok

   subroutine clear_status(status)
      type(status_type), intent(out) :: status
      status%code = mao_success
      status%message = ''
   end subroutine clear_status

   subroutine set_status(status, code, message)
      type(status_type), intent(out) :: status
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      status%code = code
      status%message = trim(message)
   end subroutine set_status

end module mao_status
