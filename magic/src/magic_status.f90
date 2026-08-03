module magic_status
   implicit none
   private

   integer, parameter, public :: MAGIC_SUCCESS = 0
   integer, parameter, public :: MAGIC_INVALID_ARGUMENT = 1
   integer, parameter, public :: MAGIC_ALLOCATION_ERROR = 2
   integer, parameter, public :: MAGIC_NOT_SUPPORTED = 3

   type, public :: magic_error
      integer :: code = MAGIC_SUCCESS
      character(len=:), allocatable :: message
   contains
      procedure :: clear => clear_error
      procedure :: failed => error_failed
   end type magic_error

   public :: set_error

contains

   subroutine clear_error(self)
      class(magic_error), intent(inout) :: self
      self%code = MAGIC_SUCCESS
      if (allocated(self%message)) deallocate(self%message)
   end subroutine clear_error

   logical function error_failed(self) result(failed)
      class(magic_error), intent(in) :: self
      failed = self%code /= MAGIC_SUCCESS
   end function error_failed

   subroutine set_error(err, code, message)
      type(magic_error), intent(inout), optional :: err
      integer, intent(in) :: code
      character(len=*), intent(in) :: message
      if (.not. present(err)) return
      err%code = code
      err%message = message
   end subroutine set_error

end module magic_status
