! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_status
  implicit none
  private

  integer, parameter, public :: sn_ok = 0
  integer, parameter, public :: sn_invalid_argument = 1
  integer, parameter, public :: sn_dimension_mismatch = 2
  integer, parameter, public :: sn_not_positive_definite = 3
  integer, parameter, public :: sn_singular_matrix = 4
  integer, parameter, public :: sn_no_convergence = 5
  integer, parameter, public :: sn_not_supported = 6
  integer, parameter, public :: sn_allocation_error = 7

  type, public :: sn_error
    integer :: code = sn_ok
    character(len=:), allocatable :: message
  contains
    procedure :: clear => sn_error_clear
    procedure :: set => sn_error_set
    procedure :: failed => sn_error_failed
  end type sn_error

contains

  subroutine sn_error_clear(self)
    class(sn_error), intent(inout) :: self
    self%code = sn_ok
    if (allocated(self%message)) deallocate(self%message)
  end subroutine sn_error_clear

  subroutine sn_error_set(self, code, message)
    class(sn_error), intent(inout) :: self
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    self%code = code
    self%message = trim(message)
  end subroutine sn_error_set

  logical function sn_error_failed(self) result(failed)
    class(sn_error), intent(in) :: self
    failed = self%code /= sn_ok
  end function sn_error_failed

end module sn_status
