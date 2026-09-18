! SPDX-License-Identifier: MIT
! SPDX-FileComment: FMM and natural cubic interpolation splines corresponding to R stats.
module r_stats_spline_interpolation
   use r_kinds, only: dp
   use r_linalg, only: solve_system
   use r_optional, only: optval
   use r_stats_interpolation, only: approxfun
   use r_stats_types, only: cubic_spline_t, interpolation_t, smooth_xy_t
   implicit none
   private

   integer, parameter, public :: spline_method_fmm = 1, spline_method_natural = 2

   public :: predict_splinefun, spline, splinefun

contains

   pure function spline(x, y, xout, method, derivative) result(value)
      !! Interpolates paired observations with an FMM or natural cubic spline.
      real(dp), intent(in) :: x(:) !! Predictor values with size `n`.
      real(dp), intent(in) :: y(:) !! Responses with size `n`.
      real(dp), intent(in) :: xout(:) !! Requested evaluation locations.
      integer, intent(in), optional :: method !! FMM or natural method; defaults to FMM.
      integer, intent(in), optional :: derivative !! Derivative order from zero through three.
      type(smooth_xy_t) :: value !! Evaluation locations and spline values or derivatives.
      type(cubic_spline_t) :: interpolator

      interpolator = splinefun(x, y, method)
      value%x = xout
      value%y = predict_splinefun(interpolator, xout, derivative)
   end function spline

   pure function splinefun(x, y, method) result(value)
      !! Constructs a reusable FMM/not-a-knot or natural cubic interpolation spline.
      real(dp), intent(in) :: x(:) !! Predictor values with size `n`; ties are averaged.
      real(dp), intent(in) :: y(:) !! Responses with size `n`.
      integer, intent(in), optional :: method !! FMM or natural method; defaults to FMM.
      type(cubic_spline_t) :: value !! Piecewise-cubic spline coefficients.
      type(interpolation_t) :: points
      real(dp), allocatable :: h(:), matrix(:, :), rhs(:), second(:)
      real(dp) :: divided, slope
      integer :: i, info, n

      value%method = optval(method, spline_method_fmm)
      if (value%method /= spline_method_fmm .and. value%method /= spline_method_natural) then
         error stop "splinefun: unsupported spline method"
      end if
      points = approxfun(x, y)
      value%x = points%x
      n = size(value%x)
      h = value%x(2:) - value%x(:n - 1)
      allocate (value%a(n - 1), value%b(n - 1), value%c(n - 1), value%d(n - 1))
      value%a = points%y(:n - 1)
      if (n == 2) then
         value%b(1) = (points%y(2) - points%y(1))/h(1)
         value%c = 0.0_dp
         value%d = 0.0_dp
      else if (n == 3 .and. value%method == spline_method_fmm) then
         slope = (points%y(2) - points%y(1))/h(1)
         divided = ((points%y(3) - points%y(2))/h(2) - slope)/(value%x(3) - value%x(1))
         value%b(1) = slope - divided*h(1)
         value%b(2) = slope + divided*h(1)
         value%c = divided
         value%d = 0.0_dp
      else
         allocate (matrix(n, n), source=0.0_dp)
         allocate (rhs(n), source=0.0_dp)
         allocate (second(n))
         if (value%method == spline_method_natural) then
            matrix(1, 1) = 1.0_dp
            matrix(n, n) = 1.0_dp
         else
            matrix(1, 1) = -h(1)
            matrix(1, 2) = h(1)
            matrix(n, n - 1) = h(n - 1)
            matrix(n, n) = -h(n - 1)
            rhs(1) = 6.0_dp*h(1)**2/(value%x(4) - value%x(1))*(&
                     ((points%y(4) - points%y(3))/h(3) - &
                       (points%y(3) - points%y(2))/h(2))/(value%x(4) - value%x(2)) - &
                     ((points%y(3) - points%y(2))/h(2) - &
                       (points%y(2) - points%y(1))/h(1))/(value%x(3) - value%x(1)))
            rhs(n) = -6.0_dp*h(n - 1)**2/(value%x(n) - value%x(n - 3))*(&
                     ((points%y(n) - points%y(n - 1))/h(n - 1) - &
                       (points%y(n - 1) - points%y(n - 2))/h(n - 2))/(value%x(n) - value%x(n - 2)) - &
                     ((points%y(n - 1) - points%y(n - 2))/h(n - 2) - &
                       (points%y(n - 2) - points%y(n - 3))/h(n - 3))/(value%x(n - 1) - value%x(n - 3)))
         end if
         do i = 2, n - 1
            matrix(i, i - 1) = h(i - 1)
            matrix(i, i) = 2.0_dp*(h(i - 1) + h(i))
            matrix(i, i + 1) = h(i)
            rhs(i) = 6.0_dp*((points%y(i + 1) - points%y(i))/h(i) - &
                             (points%y(i) - points%y(i - 1))/h(i - 1))
         end do
         call solve_system(matrix, rhs, second, info)
         if (info /= 0) error stop "splinefun: spline coefficient solution failed"
         do i = 1, n - 1
            value%b(i) = (points%y(i + 1) - points%y(i))/h(i) - &
                         h(i)*(2.0_dp*second(i) + second(i + 1))/6.0_dp
            value%c(i) = 0.5_dp*second(i)
            value%d(i) = (second(i + 1) - second(i))/(6.0_dp*h(i))
         end do
      end if
      value%left_slope = value%b(1)
      value%right_slope = value%b(n - 1) + 2.0_dp*value%c(n - 1)*h(n - 1) + &
                          3.0_dp*value%d(n - 1)*h(n - 1)**2
   end function splinefun

   pure function predict_splinefun(interpolator, xout, derivative) result(yout)
      !! Evaluates a cubic interpolation spline or its first three derivatives.
      type(cubic_spline_t), intent(in) :: interpolator !! Spline returned by `splinefun`.
      real(dp), intent(in) :: xout(:) !! Requested evaluation locations.
      integer, intent(in), optional :: derivative !! Derivative order from zero through three.
      real(dp), allocatable :: yout(:) !! Spline values or derivatives.
      real(dp) :: distance
      integer :: high, i, interval, low, middle, order

      order = optval(derivative, 0)
      if (order < 0 .or. order > 3) error stop "predict_splinefun: derivative must be in [0,3]"
      if (size(interpolator%x) < 2) error stop "predict_splinefun: invalid spline object"
      allocate (yout(size(xout)))
      do i = 1, size(xout)
         if (interpolator%method == spline_method_natural .and. xout(i) < interpolator%x(1)) then
            distance = xout(i) - interpolator%x(1)
            select case (order)
            case (0)
               yout(i) = interpolator%a(1) + interpolator%left_slope*distance
            case (1)
               yout(i) = interpolator%left_slope
            case default
               yout(i) = 0.0_dp
            end select
            cycle
         end if
         if (interpolator%method == spline_method_natural .and. &
             xout(i) > interpolator%x(size(interpolator%x))) then
            distance = xout(i) - interpolator%x(size(interpolator%x))
            select case (order)
            case (0)
               yout(i) = interval_endpoint(interpolator) + interpolator%right_slope*distance
            case (1)
               yout(i) = interpolator%right_slope
            case default
               yout(i) = 0.0_dp
            end select
            cycle
         end if
         low = 1
         high = size(interpolator%x) - 1
         if (xout(i) <= interpolator%x(1)) then
            interval = 1
         else if (xout(i) >= interpolator%x(size(interpolator%x))) then
            interval = size(interpolator%x) - 1
         else
            do while (low <= high)
               middle = (low + high)/2
               if (xout(i) < interpolator%x(middle)) then
                  high = middle - 1
               else if (xout(i) >= interpolator%x(middle + 1)) then
                  low = middle + 1
               else
                  low = middle
                  exit
               end if
            end do
            interval = low
         end if
         distance = xout(i) - interpolator%x(interval)
         select case (order)
         case (0)
            yout(i) = interpolator%a(interval) + distance*(interpolator%b(interval) + &
                      distance*(interpolator%c(interval) + distance*interpolator%d(interval)))
         case (1)
            yout(i) = interpolator%b(interval) + 2.0_dp*interpolator%c(interval)*distance + &
                      3.0_dp*interpolator%d(interval)*distance**2
         case (2)
            yout(i) = 2.0_dp*interpolator%c(interval) + 6.0_dp*interpolator%d(interval)*distance
         case (3)
            yout(i) = 6.0_dp*interpolator%d(interval)
         end select
      end do
   end function predict_splinefun

   pure function interval_endpoint(interpolator) result(value)
      !! Evaluates the final interval at its right endpoint.
      type(cubic_spline_t), intent(in) :: interpolator !! Valid cubic spline object.
      real(dp) :: value !! Response at the final knot.
      real(dp) :: distance
      integer :: interval

      interval = size(interpolator%x) - 1
      distance = interpolator%x(interval + 1) - interpolator%x(interval)
      value = interpolator%a(interval) + distance*(interpolator%b(interval) + &
              distance*(interpolator%c(interval) + distance*interpolator%d(interval)))
   end function interval_endpoint

end module r_stats_spline_interpolation
