module polyclip_integer
   use polyclip_kinds, only: dp, i8
   use polyclip_types, only: poly_path, poly_set, append_path, make_path
   implicit none
   private

   type, public :: int_path
      integer(i8), allocatable :: x(:)
      integer(i8), allocatable :: y(:)
   contains
      procedure :: size => int_path_size
   end type int_path

   type, public :: int_set
      type(int_path), allocatable :: path(:)
   contains
      procedure :: size => int_set_size
   end type int_set

   public :: scale_parameters_one, scale_parameters_two
   public :: scale_to_int, scale_from_int, append_int_path
   public :: round_away_i8, clean_int_path

contains

   pure integer function int_path_size(self) result(n)
      class(int_path), intent(in) :: self
      if (allocated(self%x)) then
         n = size(self%x)
      else
         n = 0
      end if
   end function int_path_size

   pure integer function int_set_size(self) result(n)
      class(int_set), intent(in) :: self
      if (allocated(self%path)) then
         n = size(self%path)
      else
         n = 0
      end if
   end function int_set_size

   subroutine scale_parameters_one(a, eps, x0, y0)
      type(poly_set), intent(in) :: a
      real(dp), intent(out) :: eps, x0, y0
      real(dp) :: xmin, xmax, ymin, ymax, span
      call bounds_set(a, xmin, xmax, ymin, ymax)
      span = max(xmax - xmin, ymax - ymin)
      if (span <= 0.0_dp) span = 1.0_dp
      eps = span / 1.0e9_dp
      x0 = 0.5_dp * (xmin + xmax)
      y0 = 0.5_dp * (ymin + ymax)
   end subroutine scale_parameters_one

   subroutine scale_parameters_two(a, b, eps, x0, y0)
      type(poly_set), intent(in) :: a, b
      real(dp), intent(out) :: eps, x0, y0
      real(dp) :: aminx, amaxx, aminy, amaxy
      real(dp) :: bminx, bmaxx, bminy, bmaxy
      real(dp) :: xmin, xmax, ymin, ymax, span
      call bounds_set(a, aminx, amaxx, aminy, amaxy)
      call bounds_set(b, bminx, bmaxx, bminy, bmaxy)
      xmin = min(aminx, bminx)
      xmax = max(amaxx, bmaxx)
      ymin = min(aminy, bminy)
      ymax = max(amaxy, bmaxy)
      span = max(xmax - xmin, ymax - ymin)
      if (span <= 0.0_dp) span = 1.0_dp
      eps = span / 1.0e9_dp
      x0 = 0.5_dp * (xmin + xmax)
      y0 = 0.5_dp * (ymin + ymax)
   end subroutine scale_parameters_two

   subroutine bounds_set(a, xmin, xmax, ymin, ymax)
      type(poly_set), intent(in) :: a
      real(dp), intent(out) :: xmin, xmax, ymin, ymax
      integer :: i
      logical :: first
      first = .true.
      xmin = 0.0_dp; xmax = 0.0_dp; ymin = 0.0_dp; ymax = 0.0_dp
      if (.not. allocated(a%path)) return
      do i = 1, size(a%path)
         if (.not. allocated(a%path(i)%x)) cycle
         if (size(a%path(i)%x) == 0) cycle
         if (first) then
            xmin = minval(a%path(i)%x); xmax = maxval(a%path(i)%x)
            ymin = minval(a%path(i)%y); ymax = maxval(a%path(i)%y)
            first = .false.
         else
            xmin = min(xmin, minval(a%path(i)%x)); xmax = max(xmax, maxval(a%path(i)%x))
            ymin = min(ymin, minval(a%path(i)%y)); ymax = max(ymax, maxval(a%path(i)%y))
         end if
      end do
      if (first) then
         xmin = 0.0_dp; xmax = 0.0_dp; ymin = 0.0_dp; ymax = 0.0_dp
      end if
   end subroutine bounds_set

   subroutine scale_to_int(a, x0, y0, eps, out, ierr)
      type(poly_set), intent(in) :: a
      real(dp), intent(in) :: x0, y0, eps
      type(int_set), intent(out) :: out
      integer, intent(out), optional :: ierr
      integer :: i, j, n
      real(dp) :: sx, sy, lim
      if (present(ierr)) ierr = 0
      if (eps <= 0.0_dp) then
         if (present(ierr)) then
            ierr = 1
            return
         else
            error stop "polyclip: eps must be positive"
         end if
      end if
      if (.not. allocated(a%path)) then
         allocate(out%path(0))
         return
      end if
      allocate(out%path(size(a%path)))
      lim = 0.9_dp * real(huge(0_i8), dp)
      do i = 1, size(a%path)
         n = a%path(i)%size()
         if (n /= size(a%path(i)%y)) then
            if (present(ierr)) then
               ierr = 2
               return
            else
               error stop "polyclip: x and y lengths differ"
            end if
         end if
         allocate(out%path(i)%x(n), out%path(i)%y(n))
         do j = 1, n
            sx = (a%path(i)%x(j) - x0) / eps
            sy = (a%path(i)%y(j) - y0) / eps
            if (abs(sx) > lim .or. abs(sy) > lim) then
               if (present(ierr)) then
                  ierr = 3
                  return
               else
                  error stop "polyclip: scaled coordinate exceeds int64 range"
               end if
            end if
            ! Match the package C++ interface: cInt(double) truncates toward zero.
            out%path(i)%x(j) = int(sx, i8)
            out%path(i)%y(j) = int(sy, i8)
         end do
         call clean_int_path(out%path(i), .false.)
      end do
   end subroutine scale_to_int

   subroutine scale_from_int(a, x0, y0, eps, out)
      type(int_set), intent(in) :: a
      real(dp), intent(in) :: x0, y0, eps
      type(poly_set), intent(out) :: out
      integer :: i
      if (.not. allocated(a%path)) then
         allocate(out%path(0))
         return
      end if
      allocate(out%path(size(a%path)))
      do i = 1, size(a%path)
         allocate(out%path(i)%x(a%path(i)%size()), out%path(i)%y(a%path(i)%size()))
         out%path(i)%x = x0 + eps * real(a%path(i)%x, dp)
         out%path(i)%y = y0 + eps * real(a%path(i)%y, dp)
      end do
   end subroutine scale_from_int

   pure integer(i8) function round_away_i8(x) result(v)
      real(dp), intent(in) :: x
      if (x < 0.0_dp) then
         v = int(x - 0.5_dp, i8)
      else
         v = int(x + 0.5_dp, i8)
      end if
   end function round_away_i8

   subroutine append_int_path(s, p)
      type(int_set), intent(inout) :: s
      type(int_path), intent(in) :: p
      type(int_path), allocatable :: tmp(:)
      integer :: n
      if (.not. allocated(s%path)) then
         allocate(s%path(1)); s%path(1) = p; return
      end if
      n = size(s%path)
      allocate(tmp(n + 1))
      if (n > 0) tmp(1:n) = s%path
      tmp(n + 1) = p
      call move_alloc(tmp, s%path)
   end subroutine append_int_path

   subroutine clean_int_path(p, closed)
      type(int_path), intent(inout) :: p
      logical, intent(in) :: closed
      integer(i8), allocatable :: xx(:), yy(:)
      integer :: i, k, n
      if (.not. allocated(p%x)) return
      n = size(p%x)
      if (n <= 1) return
      allocate(xx(n), yy(n))
      k = 0
      do i = 1, n
         if (k == 0) then
            k = 1
            xx(k) = p%x(i); yy(k) = p%y(i)
         else if (p%x(i) /= xx(k) .or. p%y(i) /= yy(k)) then
            k = k + 1
            xx(k) = p%x(i); yy(k) = p%y(i)
         end if
      end do
      if (closed .and. k > 1) then
         if (xx(k) == xx(1) .and. yy(k) == yy(1)) k = k - 1
      end if
      deallocate(p%x, p%y)
      allocate(p%x(k), p%y(k))
      if (k > 0) then
         p%x = xx(1:k); p%y = yy(1:k)
      end if
   end subroutine clean_int_path

end module polyclip_integer
