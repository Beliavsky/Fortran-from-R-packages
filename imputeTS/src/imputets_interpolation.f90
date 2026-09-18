module imputets_interpolation
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use imputets_kinds, only : dp
   use imputets_basic, only : na_locf
   use imputets_utils, only : apply_maxgap, count_observed, lower_string
   use stinepack_api, only : stinterp, stinterp_result
   implicit none
   private

   public :: na_interpolation

contains

   pure function na_interpolation(x, option, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series whose NaN entries are filled by interpolation.
      character(len=*), intent(in), optional :: option !! Interpolator: linear, spline, or stine; defaults to linear.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: knots_x(:)
      real(dp), allocatable :: knots_y(:)
      real(dp), allocatable :: interp(:)
      real(dp), allocatable :: xout(:)
      character(len=:), allocatable :: method
      type(stinterp_result) :: stine_result
      integer :: i
      integer :: j
      integer :: n
      integer :: nobs

      n = size(x)
      if (count_observed(x) < 2) error stop 'na_interpolation: at least two observations are required'
      method = 'linear'
      if (present(option)) method = trim(lower_string(option))
      nobs = count_observed(x)
      allocate(knots_x(nobs), knots_y(nobs), xout(n))
      j = 0
      do i = 1, n
         xout(i) = real(i, dp)
         if (.not. ieee_is_nan(x(i))) then
            j = j + 1
            knots_x(j) = real(i, dp)
            knots_y(j) = x(i)
         end if
      end do

      select case (method)
      case ('linear')
         interp = linear_interp_constant(knots_x, knots_y, xout)
      case ('spline')
         if (n == 1) then
            interp = x
         else
            do i = 1, n
               xout(i) = knots_x(1) + real(i - 1, dp) * &
                         (knots_x(nobs) - knots_x(1)) / real(n - 1, dp)
            end do
            interp = natural_cubic_spline(knots_x, knots_y, xout)
         end if
      case ('stine')
         stine_result = stinterp(knots_x, knots_y, xout)
         if (stine_result%status /= 0) error stop 'na_interpolation: stinepack interpolation failed'
         interp = stine_result%y
         if (any_nan(interp)) interp = na_locf(interp, na_remaining = 'rev')
      case default
         error stop 'na_interpolation: option must be linear, spline, or stine'
      end select

      out = x
      do i = 1, n
         if (ieee_is_nan(out(i))) out(i) = interp(i)
      end do
      call apply_maxgap(x, out, maxgap)
   end function na_interpolation

   pure function linear_interp_constant(x, y, xout) result(values)
      real(dp), intent(in) :: x(:) !! Strictly increasing interpolation knots.
      real(dp), intent(in) :: y(:) !! Ordinates at the interpolation knots; same length as x.
      real(dp), intent(in) :: xout(:) !! Coordinates to evaluate with constant endpoint extrapolation.
      real(dp), allocatable :: values(:)
      integer :: i
      integer :: j
      real(dp) :: fraction

      allocate(values(size(xout)))
      do i = 1, size(xout)
         if (xout(i) <= x(1)) then
            values(i) = y(1)
         else if (xout(i) >= x(size(x))) then
            values(i) = y(size(y))
         else
            j = 1
            do while (xout(i) > x(j + 1))
               j = j + 1
            end do
            fraction = (xout(i) - x(j)) / (x(j + 1) - x(j))
            values(i) = y(j) + fraction * (y(j + 1) - y(j))
         end if
      end do
   end function linear_interp_constant

   pure function natural_cubic_spline(x, y, xout) result(values)
      real(dp), intent(in) :: x(:) !! Strictly increasing spline knots.
      real(dp), intent(in) :: y(:) !! Ordinates at spline knots; same length as x and at least two.
      real(dp), intent(in) :: xout(:) !! Coordinates at which to evaluate the natural cubic spline.
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: second(:)
      real(dp), allocatable :: u(:)
      real(dp) :: a
      real(dp) :: b
      real(dp) :: h
      real(dp) :: p
      real(dp) :: sig
      integer :: i
      integer :: j
      integer :: k
      integer :: n

      n = size(x)
      if (n < 2 .or. size(y) /= n) error stop 'natural_cubic_spline: invalid knots'
      allocate(values(size(xout)), second(n), u(n))
      second = 0.0_dp
      u = 0.0_dp

      do i = 2, n - 1
         sig = (x(i) - x(i - 1)) / (x(i + 1) - x(i - 1))
         p = sig * second(i - 1) + 2.0_dp
         second(i) = (sig - 1.0_dp) / p
         u(i) = (6.0_dp * ((y(i + 1) - y(i)) / (x(i + 1) - x(i)) - &
                 (y(i) - y(i - 1)) / (x(i) - x(i - 1))) / (x(i + 1) - x(i - 1)) - &
                 sig * u(i - 1)) / p
      end do
      do k = n - 1, 1, -1
         second(k) = second(k) * second(k + 1) + u(k)
      end do

      do i = 1, size(xout)
         if (xout(i) <= x(1)) then
            values(i) = y(1)
            cycle
         end if
         if (xout(i) >= x(n)) then
            values(i) = y(n)
            cycle
         end if
         j = 1
         do while (xout(i) > x(j + 1))
            j = j + 1
         end do
         h = x(j + 1) - x(j)
         a = (x(j + 1) - xout(i)) / h
         b = (xout(i) - x(j)) / h
         values(i) = a * y(j) + b * y(j + 1) + &
                     ((a ** 3 - a) * second(j) + (b ** 3 - b) * second(j + 1)) * h ** 2 / 6.0_dp
      end do
   end function natural_cubic_spline

   pure logical function any_nan(x) result(found)
      real(dp), intent(in) :: x(:) !! Numeric vector tested for at least one NaN.
      integer :: i

      found = .false.
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            found = .true.
            return
         end if
      end do
   end function any_nan

end module imputets_interpolation
