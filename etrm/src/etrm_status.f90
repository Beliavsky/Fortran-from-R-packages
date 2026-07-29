! SPDX-License-Identifier: MIT
! Derived from etrm 1.0.2, Copyright (c) 2021 etrm authors.
module etrm_status
   implicit none
   private

   integer, parameter, public :: etrm_ok = 0
   integer, parameter, public :: etrm_err_size = 1
   integer, parameter, public :: etrm_err_argument = 2
   integer, parameter, public :: etrm_err_allocation = 3
   integer, parameter, public :: etrm_err_linear_solve = 4

   public :: set_status

contains

   subroutine set_status(status, message, code, text)
      integer, intent(out) :: status
      character(len=:), allocatable, intent(out) :: message
      integer, intent(in) :: code
      character(len=*), intent(in) :: text

      status = code
      message = trim(text)
   end subroutine set_status

end module etrm_status
