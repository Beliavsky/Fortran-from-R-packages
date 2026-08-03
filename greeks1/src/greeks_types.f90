! SPDX-License-Identifier: MIT
module greeks_types
  use greeks_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: greeks_ok = 0
  integer, parameter, public :: greeks_invalid_argument = 1
  integer, parameter, public :: greeks_unknown_payoff = 2
  integer, parameter, public :: greeks_unknown_greek = 3
  integer, parameter, public :: greeks_no_convergence = 4
  integer, parameter, public :: greeks_numerical_error = 5

  type, public :: greek_result
    integer :: status = greeks_ok
    integer :: iterations = 0
    character(len=256) :: message = ''
    character(len=24), allocatable :: names(:)
    real(dp), allocatable :: values(:)
    real(dp), allocatable :: standard_errors(:)
  end type greek_result

  abstract interface
    pure function payoff_callback(x, strike) result(value)
      import dp
      real(dp), intent(in) :: x, strike
      real(dp) :: value
    end function payoff_callback
  end interface
  public :: payoff_callback
  public :: initialize_result, set_error, find_name

contains

  subroutine initialize_result(result, names)
    type(greek_result), intent(out) :: result
    character(len=*), intent(in) :: names(:)
    integer :: i

    allocate(result%names(size(names)), result%values(size(names)))
    allocate(result%standard_errors(size(names)))
    result%values = 0.0_dp
    result%standard_errors = 0.0_dp
    do i = 1, size(names)
      result%names(i) = trim(names(i))
    end do
    result%status = greeks_ok
    result%iterations = 0
    result%message = ''
  end subroutine initialize_result

  subroutine set_error(result, status, message)
    type(greek_result), intent(inout) :: result
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    result%status = status
    result%message = message
  end subroutine set_error

  integer function find_name(names, target) result(index_value)
    character(len=*), intent(in) :: names(:)
    character(len=*), intent(in) :: target
    integer :: i

    index_value = 0
    do i = 1, size(names)
      if (trim(names(i)) == trim(target)) then
        index_value = i
        return
      end if
    end do
  end function find_name

end module greeks_types
