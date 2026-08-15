module countdm_data
   implicit none
   private
   public :: get_data_criminal, get_data_sbirth, data_criminal, data_sbirth
contains
   subroutine get_data_criminal(x)
      integer, allocatable, intent(out) :: x(:)
      integer :: pos
      allocate(x(4301))
      pos = 1
      call fill_value(x, pos, 4037, 0)
      call fill_value(x, pos, 219, 1)
      call fill_value(x, pos, 29, 2)
      call fill_value(x, pos, 9, 3)
      call fill_value(x, pos, 5, 4)
      call fill_value(x, pos, 2, 5)
   end subroutine get_data_criminal

   subroutine get_data_sbirth(x)
      integer, allocatable, intent(out) :: x(:)
      integer :: pos
      allocate(x(402))
      pos = 1
      call fill_value(x, pos, 314, 0)
      call fill_value(x, pos, 48, 1)
      call fill_value(x, pos, 20, 2)
      call fill_value(x, pos, 7, 3)
      call fill_value(x, pos, 5, 4)
      call fill_value(x, pos, 2, 5)
      call fill_value(x, pos, 6, 6)
   end subroutine get_data_sbirth


   function data_criminal() result(x)
      integer, allocatable :: x(:)
      call get_data_criminal(x)
   end function data_criminal

   function data_sbirth() result(x)
      integer, allocatable :: x(:)
      call get_data_sbirth(x)
   end function data_sbirth

   subroutine fill_value(x, pos, n, value)
      integer, intent(inout) :: x(:), pos
      integer, intent(in) :: n, value
      x(pos:pos+n-1) = value
      pos = pos + n
   end subroutine fill_value
end module countdm_data
