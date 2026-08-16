module polynom_status
  implicit none
  private

  integer, parameter, public :: poly_ok = 0
  integer, parameter, public :: poly_invalid_argument = 1
  integer, parameter, public :: poly_divide_by_zero = 2
  integer, parameter, public :: poly_duplicate_abscissa = 3
  integer, parameter, public :: poly_root_failure = 4

  type, public :: poly_status_t
    integer :: code = poly_ok
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded
  end type poly_status_t

  public :: set_status

contains

  logical function succeeded(self)
    class(poly_status_t), intent(in) :: self
    succeeded = self%code == poly_ok
  end function succeeded

  subroutine set_status(status, code, message)
    type(poly_status_t), intent(out), optional :: status
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    if (present(status)) then
      status%code = code
      status%message = message
    end if
  end subroutine set_status

end module polynom_status
