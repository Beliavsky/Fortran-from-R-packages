! SPDX-License-Identifier: MIT
! SPDX-FileComment: Independent modern Fortran implementation of STL decomposition.
module r_stats_stl
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use r_kinds, only: dp
   use r_stats_types, only: stl_result_t
   implicit none
   private

   public :: stl

contains

   pure function stl(x, period, s_window, s_degree, t_window, t_degree, l_window, &
                     l_degree, s_jump, t_jump, l_jump, robust, inner, outer, &
                     periodic) result(result)
      !! Decomposes a finite univariate series using seasonal-trend LOESS iterations.
      real(dp), intent(in) :: x(:) !! Finite observations spanning more than two cycles.
      integer, intent(in) :: period !! Number of observations in one seasonal cycle.
      integer, intent(in), optional :: s_window !! Seasonal LOESS window; defaults to seven.
      integer, intent(in), optional :: s_degree !! Seasonal local degree, zero or one.
      integer, intent(in), optional :: t_window !! Trend LOESS window.
      integer, intent(in), optional :: t_degree !! Trend local degree, zero or one.
      integer, intent(in), optional :: l_window !! Low-pass LOESS window.
      integer, intent(in), optional :: l_degree !! Low-pass local degree, zero or one.
      integer, intent(in), optional :: s_jump !! Seasonal smoother jump.
      integer, intent(in), optional :: t_jump !! Trend smoother jump.
      integer, intent(in), optional :: l_jump !! Low-pass smoother jump.
      logical, intent(in), optional :: robust !! Enable robustness iterations.
      integer, intent(in), optional :: inner !! Inner iterations per pass.
      integer, intent(in), optional :: outer !! Robustness iterations.
      logical, intent(in), optional :: periodic !! Make each seasonal phase constant.
      type(stl_result_t) :: result
      integer :: i, phase
      logical :: use_robust

      result%period = period
      result%window(1) = odd_window(7)
      if (present(s_window)) result%window(1) = odd_window(s_window)
      result%degree = [0, 1, 1]
      if (present(s_degree)) result%degree(1) = s_degree
      if (present(t_degree)) result%degree(2) = t_degree
      if (present(l_degree)) result%degree(3) = l_degree
      result%periodic = .false.
      if (present(periodic)) result%periodic = periodic
      if (result%periodic) then
         result%window(1) = 10*size(x) + 1
         result%degree(1) = 0
      end if
      result%window(2) = odd_window(ceiling(1.5_dp*real(period, dp)/ &
                                           (1.0_dp - 1.5_dp/real(result%window(1), dp))))
      if (present(t_window)) result%window(2) = odd_window(t_window)
      result%window(3) = odd_window(period)
      if (present(l_window)) result%window(3) = odd_window(l_window)
      result%jump = max(1, (result%window + 9)/10)
      if (present(s_jump)) result%jump(1) = max(1, s_jump)
      if (present(t_jump)) result%jump(2) = max(1, t_jump)
      if (present(l_jump)) result%jump(3) = max(1, l_jump)
      use_robust = .false.
      if (present(robust)) use_robust = robust
      result%inner = merge(1, 2, use_robust)
      if (present(inner)) result%inner = inner
      result%outer = merge(15, 0, use_robust)
      if (present(outer)) result%outer = outer

      if (period < 2 .or. size(x) <= 2*period .or. .not. all(ieee_is_finite(x))) then
         result%status = 1
         return
      end if
      if (any(result%degree < 0) .or. any(result%degree > 1) .or. &
          result%inner < 1 .or. result%outer < 0) then
         result%status = 2
         return
      end if

      allocate (result%seasonal(size(x)), result%trend(size(x)), &
                result%remainder(size(x)), result%weights(size(x)))
      call stl_iterations(x, period, result%window, result%degree, result%jump, &
                          result%inner, result%outer, result%seasonal, result%trend, &
                          result%weights)
      if (result%periodic) then
         do phase = 1, period
            result%seasonal(phase:size(x):period) = &
               sum(result%seasonal(phase:size(x):period))/ &
               real(size(result%seasonal(phase:size(x):period)), dp)
         end do
      end if
      do i = 1, size(x)
         result%remainder(i) = x(i) - result%seasonal(i) - result%trend(i)
      end do
   end function stl

   pure subroutine stl_iterations(x, period, window, degree, jump, inner, outer, &
                                  seasonal, trend, robustness_weights)
      !! Executes STL inner and outer iteration loops.
      real(dp), intent(in) :: x(:) !! Finite observations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: window(3) !! Seasonal, trend, and low-pass windows.
      integer, intent(in) :: degree(3) !! Local polynomial degrees.
      integer, intent(in) :: jump(3) !! Smoother jumps.
      integer, intent(in) :: inner !! Inner iterations per pass.
      integer, intent(in) :: outer !! Robustness iterations.
      real(dp), intent(out) :: seasonal(:) !! Estimated seasonal component.
      real(dp), intent(out) :: trend(:) !! Estimated trend component.
      real(dp), intent(out) :: robustness_weights(:) !! Final robustness weights.
      integer :: pass
      logical :: use_weights

      trend = 0.0_dp
      seasonal = 0.0_dp
      robustness_weights = 1.0_dp
      use_weights = .false.
      do pass = 0, outer
         call stl_inner_step(x, period, window, degree, jump, inner, use_weights, &
                             robustness_weights, seasonal, trend)
         if (pass == outer) exit
         call calculate_robustness_weights(x - seasonal - trend, robustness_weights)
         use_weights = .true.
      end do
   end subroutine stl_iterations

   pure subroutine stl_inner_step(x, period, window, degree, jump, iterations, &
                                  use_weights, robustness_weights, seasonal, trend)
      !! Performs seasonal smoothing, low-pass removal, and trend smoothing.
      real(dp), intent(in) :: x(:) !! Finite observations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: window(3) !! Seasonal, trend, and low-pass windows.
      integer, intent(in) :: degree(3) !! Local polynomial degrees.
      integer, intent(in) :: jump(3) !! Smoother jumps.
      integer, intent(in) :: iterations !! Number of inner iterations.
      logical, intent(in) :: use_weights !! Whether robustness weights are active.
      real(dp), intent(in) :: robustness_weights(:) !! Robustness weights.
      real(dp), intent(inout) :: seasonal(:) !! Seasonal component.
      real(dp), intent(inout) :: trend(:) !! Trend component.
      real(dp), allocatable :: detrended(:), deseasonalized(:), extended(:)
      real(dp), allocatable :: low_pass(:), ma_one(:), ma_two(:), ma_three(:)
      integer :: iteration, n

      n = size(x)
      allocate (detrended(n), deseasonalized(n), extended(n + 2*period))
      allocate (ma_one(n + period + 1), ma_two(n + 2), ma_three(n), low_pass(n))
      do iteration = 1, iterations
         detrended = x - trend
         call smooth_seasonal_subseries(detrended, period, window(1), degree(1), &
                                        jump(1), use_weights, robustness_weights, extended)
         call moving_average(extended, period, ma_one)
         call moving_average(ma_one, period, ma_two)
         call moving_average(ma_two, 3, ma_three)
         call loess_smooth(ma_three, window(3), degree(3), jump(3), .false., &
                           robustness_weights, low_pass)
         seasonal = extended(period + 1:period + n) - low_pass
         deseasonalized = x - seasonal
         call loess_smooth(deseasonalized, window(2), degree(2), jump(2), &
                           use_weights, robustness_weights, trend)
      end do
   end subroutine stl_inner_step

   pure subroutine smooth_seasonal_subseries(x, period, window, degree, jump, &
                                             use_weights, robustness_weights, extended)
      !! Smooths and extrapolates each seasonal subseries independently.
      real(dp), intent(in) :: x(:) !! Detrended observations.
      integer, intent(in) :: period !! Seasonal period.
      integer, intent(in) :: window !! Seasonal LOESS window.
      integer, intent(in) :: degree !! Local polynomial degree.
      integer, intent(in) :: jump !! Smoother jump.
      logical, intent(in) :: use_weights !! Whether robustness weights are active.
      real(dp), intent(in) :: robustness_weights(:) !! Robustness weights.
      real(dp), intent(out) :: extended(:) !! Extended seasonal sequence `(n + 2*period)`.
      real(dp), allocatable :: fit(:), subseries(:), subweights(:), work(:)
      real(dp) :: estimate
      integer :: i, n_subseries, phase
      logical :: succeeded

      do phase = 1, period
         n_subseries = (size(x) - phase)/period + 1
         allocate (subseries(n_subseries), subweights(n_subseries), &
                   fit(n_subseries), work(n_subseries))
         do i = 1, n_subseries
            subseries(i) = x(phase + (i - 1)*period)
            subweights(i) = robustness_weights(phase + (i - 1)*period)
         end do
         call loess_smooth(subseries, window, degree, jump, use_weights, subweights, fit)
         call loess_estimate(subseries, window, degree, 0.0_dp, 1, min(window, n_subseries), &
                             use_weights, subweights, work, estimate, succeeded)
         if (.not. succeeded) estimate = fit(1)
         extended(phase) = estimate
         do i = 1, n_subseries
            extended((i)*period + phase) = fit(i)
         end do
         call loess_estimate(subseries, window, degree, real(n_subseries + 1, dp), &
                             max(1, n_subseries - window + 1), n_subseries, use_weights, &
                             subweights, work, estimate, succeeded)
         if (.not. succeeded) estimate = fit(n_subseries)
         extended((n_subseries + 1)*period + phase) = estimate
         deallocate (subseries, subweights, fit, work)
      end do
   end subroutine smooth_seasonal_subseries

   pure subroutine loess_smooth(x, window, degree, jump, use_weights, &
                                robustness_weights, smoothed)
      !! Smooths equally spaced observations using STL's local regression rules.
      real(dp), intent(in) :: x(:) !! Values at consecutive integer locations.
      integer, intent(in) :: window !! Local regression window length.
      integer, intent(in) :: degree !! Local polynomial degree, zero or one.
      integer, intent(in) :: jump !! Distance between directly evaluated locations.
      logical, intent(in) :: use_weights !! Whether robustness weights are active.
      real(dp), intent(in) :: robustness_weights(:) !! Robustness weights.
      real(dp), intent(out) :: smoothed(:) !! Smoothed values with shape `(size(x))`.
      real(dp), allocatable :: work(:)
      real(dp) :: delta, estimate
      integer :: anchor, i, k, left, new_jump, next_anchor, right, shift
      logical :: succeeded

      if (size(x) < 2) then
         smoothed = x
         return
      end if
      allocate (work(size(x)))
      new_jump = min(jump, size(x) - 1)
      if (window >= size(x)) then
         left = 1
         right = size(x)
         do anchor = 1, size(x), new_jump
            call loess_estimate(x, window, degree, real(anchor, dp), left, right, &
                                use_weights, robustness_weights, work, estimate, succeeded)
            smoothed(anchor) = merge(estimate, x(anchor), succeeded)
         end do
      else if (new_jump == 1) then
         left = 1
         right = window
         shift = (window + 1)/2
         do i = 1, size(x)
            if (i > shift .and. right /= size(x)) then
               left = left + 1
               right = right + 1
            end if
            call loess_estimate(x, window, degree, real(i, dp), left, right, &
                                use_weights, robustness_weights, work, estimate, succeeded)
            smoothed(i) = merge(estimate, x(i), succeeded)
         end do
         return
      else
         shift = (window + 1)/2
         do anchor = 1, size(x), new_jump
            if (anchor < shift) then
               left = 1
               right = window
            else if (anchor >= size(x) - shift + 1) then
               left = size(x) - window + 1
               right = size(x)
            else
               left = anchor - shift + 1
               right = window + anchor - shift
            end if
            call loess_estimate(x, window, degree, real(anchor, dp), left, right, &
                                use_weights, robustness_weights, work, estimate, succeeded)
            smoothed(anchor) = merge(estimate, x(anchor), succeeded)
         end do
      end if

      if (new_jump /= 1) then
         do anchor = 1, size(x) - new_jump, new_jump
            next_anchor = anchor + new_jump
            delta = (smoothed(next_anchor) - smoothed(anchor))/real(new_jump, dp)
            do i = anchor + 1, next_anchor - 1
               smoothed(i) = smoothed(anchor) + delta*real(i - anchor, dp)
            end do
         end do
         k = ((size(x) - 1)/new_jump)*new_jump + 1
         if (k /= size(x)) then
            call loess_estimate(x, window, degree, real(size(x), dp), left, right, &
                                use_weights, robustness_weights, work, estimate, succeeded)
            smoothed(size(x)) = merge(estimate, x(size(x)), succeeded)
            delta = (smoothed(size(x)) - smoothed(k))/real(size(x) - k, dp)
            do i = k + 1, size(x) - 1
               smoothed(i) = smoothed(k) + delta*real(i - k, dp)
            end do
         end if
      end if
   end subroutine loess_smooth

   pure subroutine loess_estimate(x, window, degree, target, left, right, use_weights, &
                                  robustness_weights, weights, estimate, succeeded)
      !! Evaluates one STL local constant or local linear regression.
      real(dp), intent(in) :: x(:) !! Values at consecutive integer locations.
      integer, intent(in) :: window !! Requested local regression window.
      integer, intent(in) :: degree !! Local polynomial degree, zero or one.
      real(dp), intent(in) :: target !! Integer or extrapolated evaluation location.
      integer, intent(in) :: left !! First observation in the local window.
      integer, intent(in) :: right !! Last observation in the local window.
      logical, intent(in) :: use_weights !! Whether robustness weights are active.
      real(dp), intent(in) :: robustness_weights(:) !! Robustness weights.
      real(dp), intent(inout) :: weights(:) !! Work array for local weights.
      real(dp), intent(out) :: estimate !! Smoothed value at `target`.
      logical, intent(out) :: succeeded !! Whether the local weights have positive mass.
      real(dp) :: center, correction, distance, h, h_large, h_small
      real(dp) :: range, scale, sum_weights
      integer :: j

      range = real(size(x) - 1, dp)
      h = max(target - real(left, dp), real(right, dp) - target)
      if (window > size(x)) h = h + real((window - size(x))/2, dp)
      h_large = 0.999_dp*h
      h_small = 0.001_dp*h
      sum_weights = 0.0_dp
      weights = 0.0_dp
      do j = left, right
         distance = abs(real(j, dp) - target)
         if (distance <= h_large) then
            if (distance <= h_small) then
               weights(j) = 1.0_dp
            else
               weights(j) = (1.0_dp - (distance/h)**3)**3
            end if
            if (use_weights) weights(j) = weights(j)*robustness_weights(j)
            sum_weights = sum_weights + weights(j)
         end if
      end do
      if (sum_weights <= 0.0_dp) then
         estimate = 0.0_dp
         succeeded = .false.
         return
      end if

      weights(left:right) = weights(left:right)/sum_weights
      if (h > 0.0_dp .and. degree > 0) then
         center = 0.0_dp
         do j = left, right
            center = center + weights(j)*real(j, dp)
         end do
         correction = target - center
         scale = 0.0_dp
         do j = left, right
            scale = scale + weights(j)*(real(j, dp) - center)**2
         end do
         if (sqrt(scale) > range*0.001_dp) then
            correction = correction/scale
            do j = left, right
               weights(j) = weights(j)*(correction*(real(j, dp) - center) + 1.0_dp)
            end do
         end if
      end if
      estimate = 0.0_dp
      do j = left, right
         estimate = estimate + weights(j)*x(j)
      end do
      succeeded = .true.
   end subroutine loess_estimate

   pure subroutine moving_average(x, window, average)
      !! Computes all complete running means for a fixed window.
      real(dp), intent(in) :: x(:) !! Input values.
      integer, intent(in) :: window !! Running-mean width.
      real(dp), intent(out) :: average(:) !! Complete-window means.
      real(dp) :: total
      integer :: i

      total = sum(x(:window))
      average(1) = total/real(window, dp)
      do i = 2, size(average)
         total = total + x(i + window - 1) - x(i - 1)
         average(i) = total/real(window, dp)
      end do
   end subroutine moving_average

   pure subroutine calculate_robustness_weights(residuals, weights)
      !! Computes STL Tukey-bisquare weights from six times the median absolute residual.
      real(dp), intent(in) :: residuals(:) !! Current decomposition residuals.
      real(dp), intent(out) :: weights(:) !! Robustness weights in the closed unit interval.
      real(dp) :: absolute_residuals(size(residuals)), cutoff, lower_cutoff, upper_cutoff
      real(dp) :: ratio, temporary
      integer :: i, j, lower_middle, upper_middle

      absolute_residuals = abs(residuals)
      do i = 2, size(absolute_residuals)
         temporary = absolute_residuals(i)
         j = i - 1
         do while (j >= 1)
            if (absolute_residuals(j) <= temporary) exit
            absolute_residuals(j + 1) = absolute_residuals(j)
            j = j - 1
         end do
         absolute_residuals(j + 1) = temporary
      end do
      upper_middle = size(residuals)/2 + 1
      lower_middle = size(residuals) - upper_middle + 1
      cutoff = 3.0_dp*(absolute_residuals(lower_middle) + &
                       absolute_residuals(upper_middle))
      lower_cutoff = 0.001_dp*cutoff
      upper_cutoff = 0.999_dp*cutoff
      do i = 1, size(residuals)
         if (abs(residuals(i)) <= lower_cutoff) then
            weights(i) = 1.0_dp
         else if (abs(residuals(i)) <= upper_cutoff) then
            ratio = residuals(i)/cutoff
            weights(i) = (1.0_dp - ratio**2)**2
         else
            weights(i) = 0.0_dp
         end if
      end do
   end subroutine calculate_robustness_weights

   pure elemental integer function odd_window(value)
      !! Returns the smallest odd integer not below three and not below `value`.
      integer, intent(in) :: value !! Requested window length.

      odd_window = max(3, value)
      if (modulo(odd_window, 2) == 0) odd_window = odd_window + 1
   end function odd_window

end module r_stats_stl
