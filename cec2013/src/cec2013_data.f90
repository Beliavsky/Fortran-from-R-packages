! SPDX-License-Identifier: GPL-3.0-or-later
module cec2013_data
   use cec2013_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: CEC2013_OK = 0
   integer, parameter, public :: CEC2013_BAD_DIMENSION = 1
   integer, parameter, public :: CEC2013_IO_ERROR = 2
   integer, parameter, public :: CEC2013_BAD_PROBLEM = 3
   integer, parameter, public :: CEC2013_BAD_SHAPE = 4

   type, public :: cec2013_context
      integer :: n = 0
      real(dp), allocatable :: shift(:)
      real(dp), allocatable :: rotation(:)
   contains
      procedure :: init => cec2013_context_init
      procedure :: clear => cec2013_context_clear
      procedure :: initialized => cec2013_context_initialized
   end type cec2013_context

   public :: cec2013_dimension_supported

contains

   pure logical function cec2013_dimension_supported(n)
      integer, intent(in) :: n
      select case (n)
      case (2, 5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
         cec2013_dimension_supported = .true.
      case default
         cec2013_dimension_supported = .false.
      end select
   end function cec2013_dimension_supported

   logical function cec2013_context_initialized(self)
      class(cec2013_context), intent(in) :: self
      cec2013_context_initialized = self%n > 0 .and. allocated(self%shift) .and. allocated(self%rotation)
   end function cec2013_context_initialized

   subroutine cec2013_context_clear(self)
      class(cec2013_context), intent(inout) :: self
      if (allocated(self%shift)) deallocate(self%shift)
      if (allocated(self%rotation)) deallocate(self%rotation)
      self%n = 0
   end subroutine cec2013_context_clear

   subroutine cec2013_context_init(self, n, data_dir, status, message)
      class(cec2013_context), intent(inout) :: self
      integer, intent(in) :: n
      character(len=*), intent(in) :: data_dir
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      character(len=:), allocatable :: mfile, sfile
      character(len=16) :: ndigits
      integer :: unit, ios

      call set_status(status, message, CEC2013_OK, '')
      call self%clear()

      if (.not. cec2013_dimension_supported(n)) then
         call set_status(status, message, CEC2013_BAD_DIMENSION, 'unsupported CEC2013 dimension')
         return
      end if

      write(ndigits, '(i0)') n
      mfile = trim(data_dir) // '/M_D' // trim(ndigits) // '.txt'
      sfile = trim(data_dir) // '/shift_data.txt'

      allocate(self%rotation(10*n*n), self%shift(10*n), stat=ios)
      if (ios /= 0) then
         call self%clear()
         call set_status(status, message, CEC2013_IO_ERROR, 'allocation failure')
         return
      end if

      open(newunit=unit, file=mfile, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call self%clear()
         call set_status(status, message, CEC2013_IO_ERROR, 'cannot open rotation data: ' // mfile)
         return
      end if
      read(unit, *, iostat=ios) self%rotation
      close(unit)
      if (ios /= 0) then
         call self%clear()
         call set_status(status, message, CEC2013_IO_ERROR, 'cannot read rotation data: ' // mfile)
         return
      end if

      open(newunit=unit, file=sfile, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         call self%clear()
         call set_status(status, message, CEC2013_IO_ERROR, 'cannot open shift data: ' // sfile)
         return
      end if
      read(unit, *, iostat=ios) self%shift
      close(unit)
      if (ios /= 0) then
         call self%clear()
         call set_status(status, message, CEC2013_IO_ERROR, 'cannot read shift data: ' // sfile)
         return
      end if

      self%n = n
   end subroutine cec2013_context_init

   subroutine set_status(status, message, code, text)
      integer, intent(out), optional :: status
      character(len=*), intent(out), optional :: message
      integer, intent(in) :: code
      character(len=*), intent(in) :: text
      if (present(status)) status = code
      if (present(message)) then
         message = ''
         if (len_trim(text) > 0) message = text
      end if
   end subroutine set_status

end module cec2013_data
