! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_interpolation
   use fmbasics_kinds, only : dp, FM_OK, FM_INVALID_ARGUMENT
   implicit none
   private

   integer, parameter, public :: INTERP_CONSTANT = 1
   integer, parameter, public :: INTERP_LOG_DF = 2
   integer, parameter, public :: INTERP_LINEAR = 3
   integer, parameter, public :: INTERP_CUBIC = 4
   integer, parameter, public :: INTERP_LINEAR_CUBIC_TIME_VAR = 5

   type, public :: interpolation_t
      integer :: method = INTERP_LOG_DF
   end type interpolation_t

   public :: constant_interpolation, logdf_interpolation
   public :: linear_interpolation, cubic_interpolation
   public :: linear_cubic_time_var_interpolation
   public :: interpolate_1d, natural_spline_values

contains

   pure function constant_interpolation() result(value)
      type(interpolation_t) :: value
      value%method = INTERP_CONSTANT
   end function constant_interpolation

   pure function logdf_interpolation() result(value)
      type(interpolation_t) :: value
      value%method = INTERP_LOG_DF
   end function logdf_interpolation

   pure function linear_interpolation() result(value)
      type(interpolation_t) :: value
      value%method = INTERP_LINEAR
   end function linear_interpolation

   pure function cubic_interpolation() result(value)
      type(interpolation_t) :: value
      value%method = INTERP_CUBIC
   end function cubic_interpolation

   pure function linear_cubic_time_var_interpolation() result(value)
      type(interpolation_t) :: value
      value%method = INTERP_LINEAR_CUBIC_TIME_VAR
   end function linear_cubic_time_var_interpolation

   function interpolate_1d(method, x, y, xq, status) result(yq)
      integer, intent(in) :: method
      real(dp), intent(in) :: x(:), y(:), xq(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: yq(:)
      integer :: i, stat_i
      if (.not. valid_grid(x, y)) then
         allocate(yq(size(xq)))
         yq = 0.0_dp
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      select case (method)
      case (INTERP_CONSTANT)
         allocate(yq(size(xq)))
         do i = 1, size(xq)
            yq(i) = constant_value(x, y, xq(i))
         end do
         stat_i = FM_OK
      case (INTERP_LINEAR)
         allocate(yq(size(xq)))
         do i = 1, size(xq)
            yq(i) = linear_value(x, y, xq(i))
         end do
         stat_i = FM_OK
      case (INTERP_CUBIC)
         yq = natural_spline_values(x, y, xq, stat_i)
      case default
         allocate(yq(size(xq)))
         yq = 0.0_dp
         stat_i = FM_INVALID_ARGUMENT
      end select
      if (present(status)) status = stat_i
   end function interpolate_1d

   function natural_spline_values(x, y, xq, status) result(yq)
      real(dp), intent(in) :: x(:), y(:), xq(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: yq(:)
      real(dp), allocatable :: second(:), lower(:), diag(:), upper(:), rhs(:)
      real(dp) :: h0, h1, a, b, h
      integer :: n, i, k
      n = size(x)
      allocate(yq(size(xq)))
      if (.not. valid_grid(x, y)) then
         yq = 0.0_dp
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      if (n == 2) then
         do i = 1, size(xq)
            yq(i) = linear_value(x, y, xq(i))
         end do
         if (present(status)) status = FM_OK
         return
      end if
      allocate(second(n), lower(n), diag(n), upper(n), rhs(n))
      second = 0.0_dp
      lower = 0.0_dp
      diag = 0.0_dp
      upper = 0.0_dp
      rhs = 0.0_dp
      diag(1) = 1.0_dp
      diag(n) = 1.0_dp
      do i = 2, n - 1
         h0 = x(i) - x(i-1)
         h1 = x(i+1) - x(i)
         lower(i) = h0
         diag(i) = 2.0_dp * (h0 + h1)
         upper(i) = h1
         rhs(i) = 6.0_dp * ((y(i+1)-y(i))/h1 - (y(i)-y(i-1))/h0)
      end do
      do i = 2, n
         a = lower(i) / diag(i-1)
         diag(i) = diag(i) - a * upper(i-1)
         rhs(i) = rhs(i) - a * rhs(i-1)
      end do
      second(n) = rhs(n) / diag(n)
      do i = n - 1, 1, -1
         second(i) = (rhs(i) - upper(i) * second(i+1)) / diag(i)
      end do
      do i = 1, size(xq)
         if (xq(i) <= x(1)) then
            yq(i) = y(1)
         else if (xq(i) >= x(n)) then
            yq(i) = y(n)
         else
            k = locate_interval(x, xq(i))
            h = x(k+1) - x(k)
            a = (x(k+1) - xq(i)) / h
            b = (xq(i) - x(k)) / h
            yq(i) = a * y(k) + b * y(k+1) + &
               ((a**3-a)*second(k) + (b**3-b)*second(k+1)) * h*h / 6.0_dp
         end if
      end do
      if (present(status)) status = FM_OK
   end function natural_spline_values

   pure logical function valid_grid(x, y) result(value)
      real(dp), intent(in) :: x(:), y(:)
      value = size(x) == size(y) .and. size(x) >= 2
      if (value) value = all(x(2:) > x(:size(x)-1))
   end function valid_grid

   pure real(dp) function constant_value(x, y, q) result(value)
      real(dp), intent(in) :: x(:), y(:), q
      integer :: k
      if (q <= x(1)) then
         value = y(1)
      else if (q >= x(size(x))) then
         value = y(size(y))
      else
         k = locate_interval(x, q)
         value = y(k)
      end if
   end function constant_value

   pure real(dp) function linear_value(x, y, q) result(value)
      real(dp), intent(in) :: x(:), y(:), q
      integer :: k
      real(dp) :: w
      if (q <= x(1)) then
         value = y(1)
      else if (q >= x(size(x))) then
         value = y(size(y))
      else
         k = locate_interval(x, q)
         w = (q - x(k)) / (x(k+1) - x(k))
         value = (1.0_dp - w) * y(k) + w * y(k+1)
      end if
   end function linear_value

   pure integer function locate_interval(x, q) result(k)
      real(dp), intent(in) :: x(:), q
      integer :: lo, hi, mid
      lo = 1
      hi = size(x) - 1
      do while (lo < hi)
         mid = (lo + hi + 1) / 2
         if (x(mid) <= q) then
            lo = mid
         else
            hi = mid - 1
         end if
      end do
      k = lo
   end function locate_interval

end module fmbasics_interpolation
