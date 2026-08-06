module vamc_status
  implicit none
  private
  integer, parameter, public :: vamc_success = 0
  integer, parameter, public :: vamc_invalid_argument = 1
  integer, parameter, public :: vamc_dimension_error = 2
  integer, parameter, public :: vamc_numerical_error = 3
  integer, parameter, public :: vamc_date_error = 4
  integer, parameter, public :: vamc_not_supported = 5

  type, public :: status_type
    integer :: code = vamc_success
    character(len=:), allocatable :: message
  contains
    procedure :: ok => status_ok
    procedure :: clear => status_clear
    procedure :: set => status_set
  end type status_type
contains
  logical function status_ok(self)
    class(status_type), intent(in) :: self
    status_ok = self%code == vamc_success
  end function status_ok

  subroutine status_clear(self)
    class(status_type), intent(inout) :: self
    self%code = vamc_success
    if (allocated(self%message)) deallocate(self%message)
  end subroutine status_clear

  subroutine status_set(self, code, message)
    class(status_type), intent(inout) :: self
    integer, intent(in) :: code
    character(len=*), intent(in) :: message
    self%code = code
    self%message = trim(message)
  end subroutine status_set
end module vamc_status
