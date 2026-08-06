! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_status
  implicit none
  private
  integer, parameter, public :: waveslim_ok = 0
  integer, parameter, public :: waveslim_invalid_input = 1
  integer, parameter, public :: waveslim_invalid_filter = 2
  integer, parameter, public :: waveslim_invalid_level = 3
  integer, parameter, public :: waveslim_singular = 4
  integer, parameter, public :: waveslim_not_converged = 5
  integer, parameter, public :: waveslim_not_supported = 6

  type, public :: status_type
    integer :: code = waveslim_ok
    character(len=:), allocatable :: message
  contains
    procedure :: ok => status_ok
  end type status_type

  public :: set_status, clear_status
contains
  logical function status_ok(self)
    class(status_type), intent(in) :: self
    status_ok = self%code == waveslim_ok
  end function status_ok

  subroutine clear_status(status)
    type(status_type), intent(out) :: status
    status%code = waveslim_ok
    status%message = 'ok'
  end subroutine clear_status

  subroutine set_status(status, code, message)
    type(status_type), intent(out) :: status
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    status%code = code
    status%message = trim(message)
  end subroutine set_status
end module waveslim_status
