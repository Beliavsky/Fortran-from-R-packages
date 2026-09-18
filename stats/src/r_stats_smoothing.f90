! SPDX-License-Identifier: MIT
! SPDX-FileComment: Smoothing functions corresponding to selected R stats functions.
module r_stats_smoothing
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_ordering, only: r_order
   use r_quantiles, only: r_median
   use r_status, only: r_invalid_input, r_ok
   use r_stats_types, only: loess_fit_t, smooth_xy_t
   implicit none
   private

   public :: ksmooth, loess_fit, lowess, predict_loess, runmed

contains

   pure function runmed(x, k, endrule) result(smoothed)
      !! Computes a running median with R-compatible odd windows and endpoint rules.
      real(dp), intent(in) :: x(:) !! Input observations in their original order.
      integer, intent(in) :: k !! Odd window width; values outside the valid range return an empty vector.
      character(len=*), intent(in), optional :: endrule !! Endpoint rule: `median`, `keep`, or `constant`.
      real(dp), allocatable :: smoothed(:)
      real(dp), allocatable :: initial(:)
      character(len=:), allocatable :: rule
      integer :: half_width, i, n, window_width

      n = size(x)
      if (k < 1 .or. mod(k, 2) == 0 .or. k > n) then
         allocate (smoothed(0))
         return
      end if
      allocate (smoothed(n))
      if (n == 0 .or. k == 1) then
         smoothed = x
         return
      end if

      rule = "median"
      if (present(endrule)) rule = trim(adjustl(endrule))
      half_width = k / 2
      smoothed = x
      do i = half_width + 1, n - half_width
         smoothed(i) = r_median(x(i - half_width:i + half_width))
      end do

      select case (rule)
      case ("keep")
         return
      case ("constant")
         smoothed(1:half_width) = smoothed(half_width + 1)
         smoothed(n - half_width + 1:n) = smoothed(n - half_width)
      case ("median")
         initial = smoothed
         do i = 2, half_width
            window_width = 2*i - 1
            smoothed(i) = r_median(initial(1:window_width))
            smoothed(n - i + 1) = r_median(initial(n - window_width + 1:n))
         end do
         smoothed(1) = r_median([initial(1), smoothed(2), &
                                 smoothed(2) - 2.0_dp*(smoothed(3) - smoothed(2))])
         smoothed(n) = r_median([initial(n), smoothed(n - 1), &
                                 smoothed(n - 1) - 2.0_dp*(smoothed(n - 2) - smoothed(n - 1))])
      case default
         deallocate (smoothed)
         allocate (smoothed(0))
      end select
   end function runmed

   pure function ksmooth(x, y, kernel, bandwidth, x_points) result(smoothed)
      !! Computes Nadaraya-Watson kernel estimates using R's bandwidth conventions.
      real(dp), intent(in) :: x(:) !! Predictor values.
      real(dp), intent(in) :: y(:) !! Responses; must have the same size as `x`.
      character(len=*), intent(in), optional :: kernel !! Kernel name, `normal` or `box`; defaults to normal.
      real(dp), intent(in), optional :: bandwidth !! Positive kernel bandwidth; defaults to 0.5.
      real(dp), intent(in), optional :: x_points(:) !! Evaluation locations; defaults to 100 equally spaced points.
      type(smooth_xy_t) :: smoothed
      character(len=:), allocatable :: kernel_name
      integer, allocatable :: order(:)
      real(dp) :: bandwidth_value, maximum_x, minimum_x
      integer :: i, n_points

      if (size(x) /= size(y) .or. size(x) == 0) return
      kernel_name = "normal"
      if (present(kernel)) kernel_name = trim(adjustl(kernel))
      if (kernel_name /= "normal" .and. kernel_name /= "box") return
      bandwidth_value = 0.5_dp
      if (present(bandwidth)) bandwidth_value = bandwidth
      if (bandwidth_value <= 0.0_dp) return

      if (present(x_points)) then
         call r_order(x_points, order)
         allocate (smoothed%x(size(x_points)))
         smoothed%x = x_points(order)
      else
         n_points = max(100, size(x))
         allocate (smoothed%x(n_points))
         minimum_x = minval(x)
         maximum_x = maxval(x)
         do i = 1, n_points
            smoothed%x(i) = minimum_x + (maximum_x - minimum_x)*real(i - 1, dp)/real(n_points - 1, dp)
         end do
      end if
      allocate (smoothed%y(size(smoothed%x)))
      do i = 1, size(smoothed%x)
         smoothed%y(i) = kernel_estimate(x, y, smoothed%x(i), bandwidth_value, kernel_name)
      end do
   end function ksmooth

   pure function kernel_estimate(x, y, x_point, bandwidth, kernel) result(estimate)
      !! Evaluates one kernel-weighted local mean using R's truncated kernels.
      real(dp), intent(in) :: x(:) !! Predictor values.
      real(dp), intent(in) :: y(:) !! Responses with the same size as `x`.
      real(dp), intent(in) :: x_point !! Evaluation location.
      real(dp), intent(in) :: bandwidth !! Positive R-compatible bandwidth.
      character(len=*), intent(in) :: kernel !! Kernel name, `normal` or `box`.
      real(dp) :: estimate
      real(dp) :: scale, sum_weights, weighted_sum, weight, z
      integer :: i

      if (kernel == "box") then
         scale = 0.5_dp*bandwidth
      else
         ! R defines the normal-kernel bandwidth as its interquartile width.
         scale = bandwidth*0.3706506_dp
      end if
      sum_weights = 0.0_dp
      weighted_sum = 0.0_dp
      do i = 1, size(x)
         z = (x_point - x(i))/scale
         if (kernel == "box") then
            weight = merge(1.0_dp, 0.0_dp, abs(z) <= 1.0_dp)
         else if (abs(z) <= 4.0_dp) then
            weight = exp(-0.5_dp*z*z)
         else
            weight = 0.0_dp
         end if
         sum_weights = sum_weights + weight
         weighted_sum = weighted_sum + weight*y(i)
      end do
      if (sum_weights > 0.0_dp) then
         estimate = weighted_sum/sum_weights
      else
         estimate = ieee_value(estimate, ieee_quiet_nan)
      end if
   end function kernel_estimate

   pure function lowess(x, y, f, iter, delta) result(smoothed)
      !! Computes Cleveland's robust locally linear LOWESS smoother.
      real(dp), intent(in) :: x(:) !! Predictor values.
      real(dp), intent(in) :: y(:) !! Responses; must have the same size as `x`.
      real(dp), intent(in), optional :: f !! Neighborhood fraction; defaults to two thirds.
      integer, intent(in), optional :: iter !! Number of robustifying iterations; defaults to three.
      real(dp), intent(in), optional :: delta !! Interpolation cutoff; defaults to one percent of the x range.
      type(smooth_xy_t) :: smoothed
      real(dp), allocatable :: absolute_residuals(:), fitted(:), robust_weights(:), sorted_x(:), sorted_y(:)
      integer, allocatable :: order(:)
      real(dp) :: alpha, beta, cutoff, denominator, delta_value, fraction, h, h1, h9
      real(dp) :: median_absolute_residual, radius, sum_w, sum_wx, sum_wxx, sum_wxy, sum_wy, u, weight
      integer :: i, j, last, left, n, neighborhood_size, next_i, pass, right, robust_iterations

      n = size(x)
      if (n /= size(y) .or. n == 0) return
      call r_order(x, order)
      allocate (sorted_x(n), sorted_y(n))
      sorted_x = x(order)
      sorted_y = y(order)
      smoothed%x = sorted_x

      fraction = 2.0_dp/3.0_dp
      if (present(f)) fraction = max(0.01_dp, min(1.0_dp, f))
      robust_iterations = 3
      if (present(iter)) robust_iterations = max(0, iter)
      delta_value = 0.01_dp*(maxval(sorted_x) - minval(sorted_x))
      if (present(delta)) delta_value = max(0.0_dp, delta)
      neighborhood_size = max(2, min(n, int(fraction*real(n, dp))))
      allocate (absolute_residuals(n), fitted(n), robust_weights(n))
      robust_weights = 1.0_dp
      fitted = sorted_y

      do pass = 0, robust_iterations
         left = 1
         right = neighborhood_size
         last = 0
         i = 1
         do
            if (i > n) exit
            do while (right < n)
               if (sorted_x(i) - sorted_x(left) <= sorted_x(right + 1) - sorted_x(i)) exit
               left = left + 1
               right = right + 1
            end do
            h = max(sorted_x(i) - sorted_x(left), sorted_x(right) - sorted_x(i))
            h1 = 0.001_dp*h
            h9 = 0.999_dp*h
            sum_w = 0.0_dp
            sum_wx = 0.0_dp
            sum_wy = 0.0_dp
            sum_wxx = 0.0_dp
            sum_wxy = 0.0_dp
            do j = left, right
               radius = abs(sorted_x(j) - sorted_x(i))
               if (h <= sqrt(tiny(1.0_dp))) then
                  weight = merge(robust_weights(j), 0.0_dp, radius <= sqrt(tiny(1.0_dp)))
               else if (radius <= h9) then
                  if (radius <= h1) then
                     weight = robust_weights(j)
                  else
                     radius = radius/h
                     weight = (1.0_dp - radius**3)**3*robust_weights(j)
                  end if
               else
                  weight = 0.0_dp
               end if
               sum_w = sum_w + weight
               sum_wx = sum_wx + weight*sorted_x(j)
               sum_wy = sum_wy + weight*sorted_y(j)
               sum_wxx = sum_wxx + weight*sorted_x(j)*sorted_x(j)
               sum_wxy = sum_wxy + weight*sorted_x(j)*sorted_y(j)
            end do
            if (sum_w <= sqrt(tiny(1.0_dp))) then
               fitted(i) = sorted_y(i)
            else
               denominator = sum_w*sum_wxx - sum_wx*sum_wx
               if (abs(denominator) <= 100.0_dp*epsilon(1.0_dp)* &
                   max(1.0_dp, abs(sum_w*sum_wxx), abs(sum_wx*sum_wx))) then
                  fitted(i) = sum_wy/sum_w
               else
                  beta = (sum_w*sum_wxy - sum_wx*sum_wy)/denominator
                  alpha = (sum_wy - beta*sum_wx)/sum_w
                  fitted(i) = alpha + beta*sorted_x(i)
               end if
            end if
            if (last > 0 .and. i > last + 1) then
               do j = last + 1, i - 1
                  if (sorted_x(i) == sorted_x(last)) then
                     fitted(j) = fitted(i)
                  else
                     fitted(j) = fitted(last) + (fitted(i) - fitted(last))* &
                                 (sorted_x(j) - sorted_x(last))/(sorted_x(i) - sorted_x(last))
                  end if
               end do
            end if
            last = i
            cutoff = sorted_x(last) + delta_value
            next_i = last + 1
            do while (next_i <= n)
               if (sorted_x(next_i) > cutoff) exit
               if (sorted_x(next_i) > sorted_x(last)) fitted(next_i) = fitted(last)
               next_i = next_i + 1
            end do
            if (next_i > n) then
               i = n
            else
               i = max(last + 1, next_i - 1)
            end if
            if (i <= last) exit
         end do
         if (last < n) fitted(last + 1:n) = fitted(last)
         if (pass >= robust_iterations) exit
         absolute_residuals = abs(sorted_y - fitted)
         median_absolute_residual = r_median(absolute_residuals)
         if (median_absolute_residual <= 100.0_dp*epsilon(1.0_dp)* &
             max(1.0_dp, maxval(abs(sorted_y)))) exit
         do j = 1, n
            u = absolute_residuals(j)/(6.0_dp*median_absolute_residual)
            if (u >= 1.0_dp) then
               robust_weights(j) = 0.0_dp
            else
               robust_weights(j) = (1.0_dp - u*u)**2
            end if
         end do
      end do
      smoothed%y = fitted
   end function lowess

   pure function loess_fit(x, y, span, degree, family, weights, iterations) result(fit)
      !! Fits univariate LOESS by direct tricube-weighted local polynomial regression.
      real(dp), intent(in) :: x(:) !! Finite training predictor values.
      real(dp), intent(in) :: y(:) !! Finite training responses conformable with `x`.
      real(dp), intent(in), optional :: span !! Positive neighborhood fraction; defaults to 0.75.
      integer, intent(in), optional :: degree !! Local-polynomial degree zero, one, or two.
      character(len=*), intent(in), optional :: family !! `gaussian` or robust `symmetric` family.
      real(dp), intent(in), optional :: weights(:) !! Nonnegative prior weights conformable with `x`.
      integer, intent(in), optional :: iterations !! Symmetric-family fitting passes; defaults to four.
      type(loess_fit_t) :: fit
      real(dp), allocatable :: absolute_residuals(:)
      real(dp) :: median_absolute_residual, robust_argument
      integer :: i, pass, requested_iterations

      fit%span = optval(span, 0.75_dp)
      fit%degree = optval(degree, 2)
      fit%family = "gaussian"
      if (present(family)) fit%family = trim(adjustl(family))
      requested_iterations = optval(iterations, 4)
      if (size(x) /= size(y) .or. size(x) == 0 .or. &
          .not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y)) .or. &
          fit%span <= 0.0_dp .or. fit%degree < 0 .or. fit%degree > 2 .or. &
          (fit%family /= "gaussian" .and. fit%family /= "symmetric") .or. &
          requested_iterations < 1) then
         fit%status = r_invalid_input
         return
      end if
      if (present(weights)) then
         if (size(weights) /= size(x) .or. any(weights < 0.0_dp) .or. &
             .not. all(ieee_is_finite(weights)) .or. sum(weights) <= 0.0_dp) then
            fit%status = r_invalid_input
            return
         end if
         fit%weights = weights
      else
         allocate (fit%weights(size(x)), source=1.0_dp)
      end if
      fit%x = x
      fit%y = y
      allocate (fit%robustness(size(x)), source=1.0_dp)
      allocate (fit%fitted_values(size(x)), fit%residuals(size(x)))
      fit%iterations = merge(requested_iterations, 1, fit%family == "symmetric")
      do pass = 1, fit%iterations
         call evaluate_loess(fit%x, fit%y, fit%weights*fit%robustness, fit%span, &
                             fit%degree, fit%x, fit%fitted_values)
         fit%residuals = fit%y - fit%fitted_values
         if (pass == fit%iterations) exit
         absolute_residuals = abs(fit%residuals)
         median_absolute_residual = r_median(absolute_residuals)
         if (median_absolute_residual <= 100.0_dp*epsilon(1.0_dp)* &
             max(1.0_dp, maxval(abs(fit%y)))) exit
         do i = 1, size(x)
            robust_argument = absolute_residuals(i)/(6.0_dp*median_absolute_residual)
            if (robust_argument >= 1.0_dp) then
               fit%robustness(i) = 0.0_dp
            else
               fit%robustness(i) = (1.0_dp - robust_argument**2)**2
            end if
         end do
      end do
      fit%iterations = pass
      fit%status = r_ok
   end function loess_fit

   pure function predict_loess(fit, x_new) result(prediction)
      !! Predicts a fitted direct-surface univariate LOESS model at new locations.
      type(loess_fit_t), intent(in) :: fit !! Previously fitted LOESS model.
      real(dp), intent(in) :: x_new(:) !! Finite prediction locations.
      real(dp), allocatable :: prediction(:)

      allocate (prediction(size(x_new)))
      if (fit%status /= r_ok .or. .not. allocated(fit%x) .or. &
          .not. allocated(fit%y) .or. .not. allocated(fit%weights) .or. &
          .not. allocated(fit%robustness) .or. &
          .not. all(ieee_is_finite(x_new))) then
         prediction = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      if (size(fit%y) /= size(fit%x) .or. size(fit%weights) /= size(fit%x) .or. &
          size(fit%robustness) /= size(fit%x)) then
         prediction = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      call evaluate_loess(fit%x, fit%y, fit%weights*fit%robustness, fit%span, &
                          fit%degree, x_new, prediction)
   end function predict_loess

   pure subroutine evaluate_loess(x, y, observation_weights, span, degree, x_new, prediction)
      !! Evaluates independent direct local-polynomial regressions at requested locations.
      real(dp), intent(in) :: x(:) !! Training predictor values.
      real(dp), intent(in) :: y(:) !! Training responses.
      real(dp), intent(in) :: observation_weights(:) !! Prior times robustness weights.
      real(dp), intent(in) :: span !! Positive neighborhood fraction.
      integer, intent(in) :: degree !! Requested polynomial degree from zero through two.
      real(dp), intent(in) :: x_new(:) !! Prediction locations.
      real(dp), intent(out) :: prediction(size(x_new)) !! Local fitted values.
      integer :: i

      do i = 1, size(x_new)
         prediction(i) = local_polynomial_estimate(x, y, observation_weights, x_new(i), &
                                                   span, degree)
      end do
   end subroutine evaluate_loess

   pure function local_polynomial_estimate(x, y, observation_weights, x_new, span, degree) &
      result(estimate)
      !! Evaluates one tricube-weighted local polynomial with singular-degree fallback.
      real(dp), intent(in) :: x(:) !! Training predictor values.
      real(dp), intent(in) :: y(:) !! Training responses.
      real(dp), intent(in) :: observation_weights(:) !! Prior times robustness weights.
      real(dp), intent(in) :: x_new !! Prediction location.
      real(dp), intent(in) :: span !! Positive neighborhood fraction.
      integer, intent(in) :: degree !! Requested polynomial degree from zero through two.
      real(dp) :: estimate
      integer, allocatable :: distance_order(:)
      real(dp), allocatable :: distance(:), local_weight(:), offset(:)
      real(dp) :: determinant, denominator, radius
      real(dp) :: s0, s1, s2, s3, s4, t0, t1, t2
      integer :: i, neighborhood_size

      allocate (distance(size(x)), local_weight(size(x)), offset(size(x)))
      offset = x - x_new
      distance = abs(offset)
      call r_order(distance, distance_order)
      neighborhood_size = min(size(x), max(degree + 1, int(span*real(size(x), dp))))
      radius = distance(distance_order(neighborhood_size))
      if (span > 1.0_dp) radius = radius*span
      if (radius <= sqrt(tiny(1.0_dp))) radius = maxval(distance)
      if (radius <= sqrt(tiny(1.0_dp))) then
         local_weight = observation_weights
      else
         do i = 1, size(x)
            if (distance(i) < 0.999_dp*radius) then
               local_weight(i) = observation_weights(i)*(1.0_dp - (distance(i)/radius)**3)**3
            else
               local_weight(i) = 0.0_dp
            end if
         end do
      end if
      s0 = sum(local_weight)
      if (s0 <= sqrt(tiny(1.0_dp))) then
         do i = 1, size(x)
            if (observation_weights(distance_order(i)) > 0.0_dp) then
               estimate = y(distance_order(i))
               return
            end if
         end do
         estimate = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      t0 = dot_product(local_weight, y)
      if (degree == 0) then
         estimate = t0/s0
         return
      end if
      s1 = dot_product(local_weight, offset)
      s2 = dot_product(local_weight, offset**2)
      t1 = dot_product(local_weight, offset*y)
      denominator = s0*s2 - s1*s1
      if (degree == 1 .or. abs(denominator) <= 100.0_dp*epsilon(1.0_dp)* &
          max(1.0_dp, abs(s0*s2), abs(s1*s1))) then
         if (abs(denominator) <= 100.0_dp*epsilon(1.0_dp)* &
             max(1.0_dp, abs(s0*s2), abs(s1*s1))) then
            estimate = t0/s0
         else
            estimate = (t0*s2 - t1*s1)/denominator
         end if
         return
      end if
      s3 = dot_product(local_weight, offset**3)
      s4 = dot_product(local_weight, offset**4)
      t2 = dot_product(local_weight, offset**2*y)
      determinant = s0*(s2*s4 - s3*s3) - s1*(s1*s4 - s3*s2) + &
                    s2*(s1*s3 - s2*s2)
      if (abs(determinant) <= 100.0_dp*epsilon(1.0_dp)* &
          max(1.0_dp, abs(s0*s2*s4))) then
         estimate = (t0*s2 - t1*s1)/denominator
      else
         estimate = (t0*(s2*s4 - s3*s3) - s1*(t1*s4 - s3*t2) + &
                     s2*(t1*s3 - s2*t2))/determinant
      end if
   end function local_polynomial_estimate

end module r_stats_smoothing
