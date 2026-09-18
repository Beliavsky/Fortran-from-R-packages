module stinepack_api
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use stinepack_kinds, only : dp
   implicit none
   private

   public :: dp
   public :: parabola_slopes
   public :: stineman_slopes
   public :: stinterp
   public :: na_stinterp
   public :: na_stinterp_vector
   public :: na_stinterp_matrix
   public :: stinterp_result
   public :: vector_interp_result
   public :: matrix_interp_result

   type :: stinterp_result
      real(dp), allocatable :: x(:)       !! Requested interpolation coordinates.
      real(dp), allocatable :: y(:)       !! Interpolated ordinates; out-of-range values are NaN.
      integer :: status = 0               !! Zero on success; nonzero for invalid inputs.
      character(len=:), allocatable :: message !! Diagnostic text when status is nonzero.
   end type stinterp_result

   type :: vector_interp_result
      real(dp), allocatable :: values(:)  !! Gap-filled vector, optionally with unresolved edge NaNs removed.
      integer :: status = 0               !! Zero on success; nonzero for invalid inputs.
      character(len=:), allocatable :: message !! Diagnostic text when status is nonzero.
   end type vector_interp_result

   type :: matrix_interp_result
      real(dp), allocatable :: values(:, :) !! Columnwise gap-filled matrix, optionally with NaN rows removed.
      integer :: status = 0                 !! Zero on success; nonzero for invalid inputs.
      character(len=:), allocatable :: message !! Diagnostic text when status is nonzero.
   end type matrix_interp_result

   interface na_stinterp
      module procedure na_stinterp_vector
      module procedure na_stinterp_matrix
   end interface na_stinterp

contains

   pure function parabola_slopes(x, y) result(yp)
      real(dp), intent(in) :: x(:) !! Strictly increasing knot coordinates; size must equal size(y) and be at least two.
      real(dp), intent(in) :: y(:) !! Ordinates at x; same shape as x.
      real(dp), allocatable :: yp(:)

      integer :: i
      integer :: m
      real(dp), allocatable :: dx(:)
      real(dp), allocatable :: dydx(:)

      m = size(x)
      if (m /= size(y) .or. m < 2) then
         allocate (yp(0))
         return
      end if

      allocate (yp(m))
      dx = x(2:m) - x(1:m - 1)
      dydx = (y(2:m) - y(1:m - 1)) / dx

      if (m == 2) then
         yp = dydx(1)
         return
      end if

      yp(1) = (dydx(1) * (2.0_dp * dx(1) + dx(2)) - dydx(2) * dx(1)) / (dx(1) + dx(2))
      do i = 2, m - 1
         yp(i) = (dydx(i - 1) * dx(i) + dydx(i) * dx(i - 1)) / (dx(i) + dx(i - 1))
      end do
      yp(m) = (dydx(m - 1) * (2.0_dp * dx(m - 1) + dx(m - 2)) - &
               dydx(m - 2) * dx(m - 1)) / (dx(m - 1) + dx(m - 2))
   end function parabola_slopes

   pure function stineman_slopes(x, y, scale) result(yp)
      real(dp), intent(in) :: x(:) !! Strictly increasing knot coordinates; size must equal size(y) and be at least two.
      real(dp), intent(in) :: y(:) !! Ordinates at x; same shape as x.
      logical, intent(in), optional :: scale !! If true, normalize x and y ranges before estimating the slopes.
      real(dp), allocatable :: yp(:)

      logical :: do_scale
      integer :: i
      integer :: m
      real(dp) :: denominator
      real(dp) :: dx2dy2m
      real(dp) :: dx2dy2p
      real(dp) :: s
      real(dp) :: sx
      real(dp) :: sy
      real(dp), allocatable :: dx(:)
      real(dp), allocatable :: dy(:)
      real(dp), allocatable :: xs(:)
      real(dp), allocatable :: ys(:)

      m = size(x)
      if (m /= size(y) .or. m < 2) then
         allocate (yp(0))
         return
      end if

      allocate (yp(m))
      if (m == 2) then
         yp = (y(2) - y(1)) / (x(2) - x(1))
         return
      end if

      do_scale = .false.
      if (present(scale)) do_scale = scale
      allocate (xs(m), ys(m))

      if (do_scale) then
         sx = maxval(x) - minval(x)
         sy = maxval(y) - minval(y)
         if (sy <= 0.0_dp) sy = 1.0_dp
         xs = x / sx
         ys = y / sy
      else
         sx = 1.0_dp
         sy = 1.0_dp
         xs = x
         ys = y
      end if

      dx = xs(2:m) - xs(1:m - 1)
      dy = ys(2:m) - ys(1:m - 1)

      do i = 2, m - 1
         dx2dy2m = dx(i - 1)**2 + dy(i - 1)**2
         dx2dy2p = dx(i)**2 + dy(i)**2
         denominator = dx(i - 1) * dx2dy2p + dx(i) * dx2dy2m
         yp(i) = (dy(i - 1) * dx2dy2p + dy(i) * dx2dy2m) / denominator
      end do

      s = dy(1) / dx(1)
      if ((s >= 0.0_dp .and. s >= yp(2)) .or. (s <= 0.0_dp .and. s <= yp(2))) then
         yp(1) = 2.0_dp * s - yp(2)
      else
         yp(1) = s + abs(s) * (s - yp(2)) / (abs(s) + abs(s - yp(2)))
      end if

      s = dy(m - 1) / dx(m - 1)
      if ((s >= 0.0_dp .and. s >= yp(m - 1)) .or. (s <= 0.0_dp .and. s <= yp(m - 1))) then
         yp(m) = 2.0_dp * s - yp(m - 1)
      else
         yp(m) = s + abs(s) * (s - yp(m - 1)) / (abs(s) + abs(s - yp(m - 1)))
      end if

      if (do_scale) yp = yp * sy / sx
   end function stineman_slopes

   pure function stinterp(x, y, xout, yp, method) result(res)
      real(dp), intent(in) :: x(:) !! Strictly increasing knot coordinates; at least two finite/non-NaN values.
      real(dp), intent(in) :: y(:) !! Ordinates at x; same length as x and containing no NaNs.
      real(dp), intent(in) :: xout(:) !! Coordinates at which to evaluate the interpolant; NaNs are invalid.
      real(dp), intent(in), optional :: yp(:) !! Optional slopes at x; when present, method must be omitted.
      character(len=*), intent(in), optional :: method !! Slope method or unique prefix: scaledstineman, stineman, or parabola.
      type(stinterp_result) :: res

      integer :: code
      integer :: i
      integer :: ix
      integer :: m
      integer :: k
      real(dp) :: dxo1
      real(dp) :: dxo2
      real(dp) :: dyo1
      real(dp) :: dyo2
      real(dp) :: product
      real(dp) :: epx
      real(dp) :: y0o
      real(dp), allocatable :: dx(:)
      real(dp), allocatable :: dy(:)
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: slopes(:)

      m = size(x)
      k = size(xout)
      allocate (res%x(k))
      res%x = xout
      allocate (res%y(k))
      res%y = quiet_nan()
      res%status = 0
      res%message = ''

      if (m < 2) then
         call set_stinterp_error(res, 1, 'x must have 2 or more elements')
         return
      end if
      if (size(y) /= m) then
         call set_stinterp_error(res, 2, 'x must have the same number of elements as y')
         return
      end if
      if (any(ieee_is_nan(x)) .or. any(ieee_is_nan(y)) .or. any(ieee_is_nan(xout))) then
         call set_stinterp_error(res, 3, 'NaNs in x, y, or xout are not allowed')
         return
      end if
      if (any(x(2:m) - x(1:m - 1) <= 0.0_dp)) then
         call set_stinterp_error(res, 4, 'The values of x must be strictly increasing')
         return
      end if

      allocate (slopes(m))
      if (present(yp)) then
         if (present(method)) then
            call set_stinterp_error(res, 5, 'method must not be specified when yp is given')
            return
         end if
         if (size(yp) /= m) then
            call set_stinterp_error(res, 6, 'yp must have the same number of elements as y')
            return
         end if
         if (any(ieee_is_nan(yp))) then
            call set_stinterp_error(res, 7, 'NaNs in yp are not allowed')
            return
         end if
         slopes = yp
      else
         if (present(method)) then
            code = match_method(method)
         else
            code = 1
         end if
         select case (code)
         case (1)
            slopes = stineman_slopes(x, y, scale=.true.)
         case (2)
            slopes = stineman_slopes(x, y, scale=.false.)
         case (3)
            slopes = parabola_slopes(x, y)
         case default
            call set_stinterp_error(res, 8, 'method must uniquely match scaledstineman, stineman, or parabola')
            return
         end select
      end if

      dx = x(2:m) - x(1:m - 1)
      dy = y(2:m) - y(1:m - 1)
      s = dy / dx
      epx = 5.0_dp * epsilon(1.0_dp) * (x(m) - x(1))

      do i = 1, k
         ix = interval_index(xout(i), x, epx)
         if (ix < 1 .or. ix > m - 1) cycle

         dxo1 = xout(i) - x(ix)
         dxo2 = xout(i) - x(ix + 1)
         y0o = y(ix) + s(ix) * dxo1
         dyo1 = (slopes(ix) - s(ix)) * dxo1
         dyo2 = (slopes(ix + 1) - s(ix)) * dxo2
         product = dyo1 * dyo2
         res%y(i) = y0o

         if (m > 2 .or. present(yp)) then
            if (product > 0.0_dp) then
               res%y(i) = y0o + product / (dyo1 + dyo2)
            else if (product < 0.0_dp) then
               res%y(i) = y0o + product * (dxo1 + dxo2) / (dyo1 - dyo2) / dx(ix)
            end if
         end if
      end do
   end function stinterp

   pure function na_stinterp_vector(object, along, na_rm, method) result(res)
      real(dp), intent(in) :: object(:) !! Vector whose NaNs are to be replaced by Stineman interpolation.
      real(dp), intent(in), optional :: along(:) !! Coordinates for rows; defaults to 1, 2, ..., size(object).
      logical, intent(in), optional :: na_rm !! If true (default), remove unresolved leading/trailing NaNs.
      character(len=*), intent(in), optional :: method !! Slope method or unique prefix passed to stinterp.
      type(vector_interp_result) :: res

      logical :: remove_na
      integer :: i
      integer :: n
      integer :: n_good
      integer :: n_miss
      integer :: n_keep
      integer :: j
      logical, allocatable :: missing(:)
      real(dp), allocatable :: coords(:)
      real(dp), allocatable :: filled(:)
      real(dp), allocatable :: xgood(:)
      real(dp), allocatable :: xmiss(:)
      real(dp), allocatable :: ygood(:)
      type(stinterp_result) :: interp

      n = size(object)
      res%status = 0
      res%message = ''

      if (present(along)) then
         if (size(along) /= n) then
            call set_vector_error(res, 1, 'along must have the same number of elements as object')
            return
         end if
         coords = along
      else
         allocate (coords(n))
         do i = 1, n
            coords(i) = real(i, dp)
         end do
      end if

      if (any(ieee_is_nan(coords))) then
         call set_vector_error(res, 2, 'NaNs in along are not allowed')
         return
      end if

      remove_na = .true.
      if (present(na_rm)) remove_na = na_rm

      filled = object
      missing = ieee_is_nan(object)
      n_miss = count(missing)

      if (n_miss > 0) then
         n_good = n - n_miss
         if (n_good < 2) then
            call set_vector_error(res, 3, 'at least two non-NaN values are required for interpolation')
            return
         end if

         allocate (xgood(n_good), ygood(n_good), xmiss(n_miss))
         xgood = pack(coords, .not. missing)
         ygood = pack(object, .not. missing)
         xmiss = pack(coords, missing)
         if (present(method)) then
            interp = stinterp(xgood, ygood, xmiss, method=method)
         else
            interp = stinterp(xgood, ygood, xmiss)
         end if
         if (interp%status /= 0) then
            call set_vector_error(res, 10 + interp%status, interp%message)
            return
         end if

         j = 0
         do i = 1, n
            if (missing(i)) then
               j = j + 1
               filled(i) = interp%y(j)
            end if
         end do
      end if

      if (remove_na) then
         n_keep = count(.not. ieee_is_nan(filled))
         allocate (res%values(n_keep))
         res%values = pack(filled, .not. ieee_is_nan(filled))
      else
         res%values = filled
      end if
   end function na_stinterp_vector

   pure function na_stinterp_matrix(object, along, na_rm, method) result(res)
      real(dp), intent(in) :: object(:, :) !! Matrix whose columns are gap-filled independently; rows share along coordinates.
      real(dp), intent(in), optional :: along(:) !! Coordinates for rows; defaults to 1, 2, ..., size(object,1).
      logical, intent(in), optional :: na_rm !! If true (default), remove rows retaining a NaN after interpolation.
      character(len=*), intent(in), optional :: method !! Slope method or unique prefix passed to each column interpolation.
      type(matrix_interp_result) :: res

      logical :: remove_na
      integer :: i
      integer :: j
      integer :: ncol
      integer :: nrow
      integer :: nkeep
      logical, allocatable :: keep(:)
      real(dp), allocatable :: coords(:)
      real(dp), allocatable :: work(:, :)
      type(vector_interp_result) :: column

      nrow = size(object, 1)
      ncol = size(object, 2)
      res%status = 0
      res%message = ''

      if (present(along)) then
         if (size(along) /= nrow) then
            call set_matrix_error(res, 1, 'along must have one element per matrix row')
            return
         end if
         coords = along
      else
         allocate (coords(nrow))
         do i = 1, nrow
            coords(i) = real(i, dp)
         end do
      end if

      if (any(ieee_is_nan(coords))) then
         call set_matrix_error(res, 2, 'NaNs in along are not allowed')
         return
      end if

      remove_na = .true.
      if (present(na_rm)) remove_na = na_rm
      work = object

      do j = 1, ncol
         if (present(method)) then
            column = na_stinterp_vector(object(:, j), along=coords, na_rm=.false., method=method)
         else
            column = na_stinterp_vector(object(:, j), along=coords, na_rm=.false.)
         end if
         if (column%status /= 0) then
            call set_matrix_error(res, 10 + column%status, column%message)
            return
         end if
         work(:, j) = column%values
      end do

      if (.not. remove_na) then
         res%values = work
         return
      end if

      allocate (keep(nrow))
      keep = .true.
      do i = 1, nrow
         keep(i) = .not. any(ieee_is_nan(work(i, :)))
      end do
      nkeep = count(keep)
      allocate (res%values(nkeep, ncol))
      j = 0
      do i = 1, nrow
         if (keep(i)) then
            j = j + 1
            res%values(j, :) = work(i, :)
         end if
      end do
   end function na_stinterp_matrix

   pure integer function interval_index(xout, x, epx) result(ix)
      real(dp), intent(in) :: xout !! Single requested interpolation coordinate.
      real(dp), intent(in) :: x(:) !! Strictly increasing interpolation knots.
      real(dp), intent(in) :: epx !! Tiny endpoint extrapolation tolerance based on machine precision.

      integer :: hi
      integer :: lo
      integer :: mid
      integer :: m

      m = size(x)
      if (xout >= x(1) - epx .and. xout <= x(1)) then
         ix = 1
         return
      end if
      if (xout >= x(m) .and. xout <= x(m) + epx) then
         ix = m - 1
         return
      end if
      if (xout < x(1)) then
         ix = 0
         return
      end if
      if (xout > x(m)) then
         ix = m
         return
      end if

      lo = 1
      hi = m
      do while (lo + 1 < hi)
         mid = (lo + hi) / 2
         if (xout < x(mid)) then
            hi = mid
         else
            lo = mid
         end if
      end do
      ix = lo
   end function interval_index

   pure integer function match_method(method) result(code)
      character(len=*), intent(in) :: method !! Candidate method name, interpreted with R-style unique-prefix matching.

      character(len=*), parameter :: names(3) = [character(len=15) :: &
                                                  'scaledstineman', 'stineman', 'parabola']
      integer :: i
      integer :: nmatch
      integer :: n

      n = len_trim(method)
      code = 0
      if (n == 0) return

      nmatch = 0
      do i = 1, size(names)
         if (n <= len_trim(names(i))) then
            if (method(1:n) == names(i)(1:n)) then
               code = i
               nmatch = nmatch + 1
            end if
         end if
      end do
      if (nmatch /= 1) code = 0
   end function match_method

   pure function quiet_nan() result(value)
      real(dp) :: value

      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function quiet_nan

   pure subroutine set_stinterp_error(res, status, message)
      type(stinterp_result), intent(inout) :: res !! Interpolation result whose diagnostic state is to be set.
      integer, intent(in) :: status !! Nonzero package-defined error code.
      character(len=*), intent(in) :: message !! Human-readable explanation of the invalid input.

      res%status = status
      res%message = message
   end subroutine set_stinterp_error

   pure subroutine set_vector_error(res, status, message)
      type(vector_interp_result), intent(inout) :: res !! Vector gap-fill result whose diagnostic state is to be set.
      integer, intent(in) :: status !! Nonzero package-defined error code.
      character(len=*), intent(in) :: message !! Human-readable explanation of the invalid input.

      res%status = status
      res%message = message
      allocate (res%values(0))
   end subroutine set_vector_error

   pure subroutine set_matrix_error(res, status, message)
      type(matrix_interp_result), intent(inout) :: res !! Matrix gap-fill result whose diagnostic state is to be set.
      integer, intent(in) :: status !! Nonzero package-defined error code.
      character(len=*), intent(in) :: message !! Human-readable explanation of the invalid input.

      res%status = status
      res%message = message
      allocate (res%values(0, 0))
   end subroutine set_matrix_error

end module stinepack_api
