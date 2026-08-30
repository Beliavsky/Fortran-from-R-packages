module polyclip_types
   use polyclip_kinds, only: dp
   implicit none
   private

   integer, parameter, public :: clip_intersection = 1
   integer, parameter, public :: clip_union        = 2
   integer, parameter, public :: clip_difference   = 3
   integer, parameter, public :: clip_xor          = 4

   integer, parameter, public :: fill_evenodd  = 1
   integer, parameter, public :: fill_nonzero  = 2
   integer, parameter, public :: fill_positive = 3
   integer, parameter, public :: fill_negative = 4

   integer, parameter, public :: join_square = 1
   integer, parameter, public :: join_round  = 2
   integer, parameter, public :: join_miter  = 3

   integer, parameter, public :: end_closed_polygon = 1
   integer, parameter, public :: end_closed_line    = 2
   integer, parameter, public :: end_open_butt      = 3
   integer, parameter, public :: end_open_square    = 4
   integer, parameter, public :: end_open_round     = 5

   type, public :: poly_path
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
   contains
      procedure :: size => poly_path_size
   end type poly_path

   type, public :: poly_set
      type(poly_path), allocatable :: path(:)
   contains
      procedure :: size => poly_set_size
   end type poly_set

   public :: make_path, make_set, clear_set, append_path

contains

   pure integer function poly_path_size(self) result(n)
      class(poly_path), intent(in) :: self
      if (allocated(self%x)) then
         n = size(self%x)
      else
         n = 0
      end if
   end function poly_path_size

   pure integer function poly_set_size(self) result(n)
      class(poly_set), intent(in) :: self
      if (allocated(self%path)) then
         n = size(self%path)
      else
         n = 0
      end if
   end function poly_set_size

   function make_path(x, y) result(p)
      real(dp), intent(in) :: x(:), y(:)
      type(poly_path) :: p
      if (size(x) /= size(y)) error stop "polyclip: x and y lengths differ"
      allocate(p%x(size(x)), p%y(size(y)))
      p%x = x
      p%y = y
   end function make_path

   function make_set(paths) result(s)
      type(poly_path), intent(in) :: paths(:)
      type(poly_set) :: s
      allocate(s%path(size(paths)))
      s%path = paths
   end function make_set

   subroutine clear_set(s)
      type(poly_set), intent(inout) :: s
      if (allocated(s%path)) deallocate(s%path)
   end subroutine clear_set

   subroutine append_path(s, p)
      type(poly_set), intent(inout) :: s
      type(poly_path), intent(in) :: p
      type(poly_path), allocatable :: tmp(:)
      integer :: n
      if (.not. allocated(s%path)) then
         allocate(s%path(1))
         s%path(1) = p
         return
      end if
      n = size(s%path)
      allocate(tmp(n + 1))
      if (n > 0) tmp(1:n) = s%path
      tmp(n + 1) = p
      call move_alloc(tmp, s%path)
   end subroutine append_path

end module polyclip_types
