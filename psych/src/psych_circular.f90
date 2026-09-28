module psych_circular
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp, r_pi
   use psych_stats, only : finite_mean, finite_sd, finite_correlation, finite_variance
   use psych_types, only : circadian_result, cosinor_result
   implicit none
   private
   public :: circular_mean_value, circular_correlation_matrix, circadian_stats_vector
   public :: circadian_correlation_matrix, circadian_linear_correlation
   public :: cosinor_fit, cosinor_period_fit, circadian_f_test, circ_tests_values
contains
   pure real(dp) function circular_mean_value(angle, hours, period) result(mean_angle)
      real(dp), intent(in) :: angle(:) !! Circular observations in radians or clock units.
      logical, intent(in), optional :: hours !! Interpret input/output as clock units when true.
      real(dp), intent(in), optional :: period !! Clock period when hours is true; defaults to 24.
      real(dp) :: mx, my, scale, use_period
      integer :: i, n
      logical :: use_hours
      mx = 0.0_dp
      my = 0.0_dp
      n = 0
      use_hours = .false.
      if (present(hours)) use_hours = hours
      use_period = 24.0_dp
      if (present(period)) use_period = period
      scale = merge(2.0_dp * r_pi / use_period, 1.0_dp, use_hours)
      do i = 1, size(angle)
         if (ieee_is_nan(angle(i))) cycle
         mx = mx + cos(angle(i) * scale)
         my = my + sin(angle(i) * scale)
         n = n + 1
      end do
      if (n == 0) then
         mean_angle = 0.0_dp
         return
      end if
      mean_angle = atan2(my / real(n, dp), mx / real(n, dp))
      if (use_hours) then
         mean_angle = modulo(mean_angle * use_period / (2.0_dp * r_pi), use_period)
      end if
   end function circular_mean_value

   pure subroutine circadian_stats_vector(angle, result, hours, period)
      real(dp), intent(in) :: angle(:) !! Circular observations in radians or clock units; IEEE NaNs are ignored.
      type(circadian_result), intent(out) :: result !! Circular center, dispersion, resultant, Rayleigh z, p, and sample size.
      logical, intent(in), optional :: hours !! Interpret input as clock units when true.
      real(dp), intent(in), optional :: period !! Clock period when hours is true; defaults to 24.
      real(dp) :: sx, sy, scale, use_period, r
      integer :: i, n
      logical :: use_hours
      use_hours = .false.
      if (present(hours)) use_hours = hours
      use_period = 24.0_dp
      if (present(period)) use_period = period
      scale = merge(2.0_dp * r_pi / use_period, 1.0_dp, use_hours)
      sx = 0.0_dp
      sy = 0.0_dp
      n = 0
      do i = 1, size(angle)
         if (ieee_is_nan(angle(i))) cycle
         sx = sx + cos(angle(i) * scale)
         sy = sy + sin(angle(i) * scale)
         n = n + 1
      end do
      result%n = n
      if (n == 0) then
         result%mean_phase = 0.0_dp
         result%sd_phase = 0.0_dp
         result%resultant = 0.0_dp
         result%z = 0.0_dp
         result%p = 1.0_dp
         return
      end if
      r = sqrt(sx * sx + sy * sy) / real(n, dp)
      result%mean_phase = atan2(sy, sx)
      result%sd_phase = sqrt(max(0.0_dp, -2.0_dp * log(max(r, tiny(1.0_dp)))))
      if (use_hours) then
         result%mean_phase = modulo(result%mean_phase * use_period / (2.0_dp * r_pi), use_period)
         result%sd_phase = result%sd_phase * use_period / (2.0_dp * r_pi)
      end if
      result%resultant = r
      result%z = real(n, dp) * r * r
      result%p = exp(-result%z)
   end subroutine circadian_stats_vector

   pure subroutine circular_correlation_matrix(angle, correlation, variance)
      real(dp), intent(in) :: angle(:, :) !! Observation-by-variable circular angles in radians; IEEE NaNs are ignored pairwise.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Circular sine-centered correlation matrix.
      real(dp), allocatable, intent(out), optional :: variance(:) !! Unnormalized sine-deviation sums of squares.
      real(dp), allocatable :: mean_angle(:), sums(:, :)
      real(dp) :: xi, xj, denom
      integer :: i, j, k, p
      p = size(angle, 2)
      allocate(mean_angle(p), sums(p, p), correlation(p, p))
      do j = 1, p
         mean_angle(j) = circular_mean_value(angle(:, j))
      end do
      sums = 0.0_dp
      do j = 1, p
         do k = 1, p
            do i = 1, size(angle, 1)
               if (ieee_is_nan(angle(i, j)) .or. ieee_is_nan(angle(i, k))) cycle
               xi = sin(angle(i, j) - mean_angle(j))
               xj = sin(angle(i, k) - mean_angle(k))
               sums(j, k) = sums(j, k) + xi * xj
            end do
         end do
      end do
      do j = 1, p
         do k = 1, p
            denom = sqrt(max(sums(j, j), 0.0_dp) * max(sums(k, k), 0.0_dp))
            if (denom > tiny(1.0_dp)) then
               correlation(j, k) = sums(j, k) / denom
            else
               correlation(j, k) = 0.0_dp
            end if
         end do
      end do
      if (present(variance)) then
         allocate(variance(p))
         do j = 1, p
            variance(j) = sums(j, j)
         end do
      end if
   end subroutine circular_correlation_matrix

   pure subroutine circadian_correlation_matrix(angle, correlation, hours, period)
      real(dp), intent(in) :: angle(:, :) !! Observation-by-variable circular values in radians or clock units.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Circular correlation matrix.
      logical, intent(in), optional :: hours !! Interpret values as clock units when true.
      real(dp), intent(in), optional :: period !! Clock period when hours is true; defaults to 24.
      real(dp), allocatable :: converted(:, :)
      real(dp) :: use_period
      logical :: use_hours
      use_hours = .false.
      if (present(hours)) use_hours = hours
      use_period = 24.0_dp
      if (present(period)) use_period = period
      allocate(converted(size(angle, 1), size(angle, 2)))
      if (use_hours) then
         converted = angle * (2.0_dp * r_pi / use_period)
      else
         converted = angle
      end if
      call circular_correlation_matrix(converted, correlation)
   end subroutine circadian_correlation_matrix

   pure real(dp) function circadian_linear_correlation(angle, x, hours, period) result(value)
      real(dp), intent(in) :: angle(:) !! Circular predictor in radians or clock units.
      real(dp), intent(in) :: x(:) !! Linear response matched to angle observations.
      logical, intent(in), optional :: hours !! Interpret angle as clock units when true.
      real(dp), intent(in), optional :: period !! Clock period when hours is true; defaults to 24.
      real(dp), allocatable :: a(:), c(:), s(:)
      real(dp) :: rc, rs, rcs, scale, use_period, radicand
      logical :: use_hours
      use_hours = .false.
      if (present(hours)) use_hours = hours
      use_period = 24.0_dp
      if (present(period)) use_period = period
      scale = merge(2.0_dp * r_pi / use_period, 1.0_dp, use_hours)
      allocate(a(size(angle)), c(size(angle)), s(size(angle)))
      a = angle * scale
      c = cos(a)
      s = sin(a)
      rc = finite_correlation(c, x)
      rs = finite_correlation(s, x)
      rcs = finite_correlation(c, s)
      radicand = (rc * rc + rs * rs - 2.0_dp * rc * rs * rcs) / max(1.0_dp - rcs * rcs, tiny(1.0_dp))
      value = sqrt(max(radicand, 0.0_dp))
      if (rc < 0.0_dp) value = -value
   end function circadian_linear_correlation

   pure subroutine cosinor_fit(angle, y, result, hours, period)
      real(dp), intent(in) :: angle(:) !! Measurement times in radians or clock units.
      real(dp), intent(in) :: y(:) !! Observed response values matched to measurement times.
      type(cosinor_result), intent(out) :: result !! Cosinor phase, fit, amplitude, dispersion, mean, and intercept.
      logical, intent(in), optional :: hours !! Interpret angle as clock units when true.
      real(dp), intent(in), optional :: period !! Clock period when hours is true; defaults to 24.
      real(dp), allocatable :: a(:), c(:), s(:), fitted(:)
      real(dp) :: rcy, rsy, rcs, beta1, beta2, scale, use_period, sy, phase_rad
      logical :: use_hours
      use_hours = .false.
      if (present(hours)) use_hours = hours
      use_period = 24.0_dp
      if (present(period)) use_period = period
      scale = merge(2.0_dp * r_pi / use_period, 1.0_dp, use_hours)
      allocate(a(size(angle)), c(size(angle)), s(size(angle)), fitted(size(angle)))
      a = angle * scale
      c = cos(a)
      s = sin(a)
      rcy = finite_correlation(c, y)
      rsy = finite_correlation(s, y)
      rcs = finite_correlation(c, s)
      beta1 = (rcy - rsy * rcs) / max(1.0_dp - rcs * rcs, tiny(1.0_dp))
      beta2 = (rsy - rcy * rcs) / max(1.0_dp - rcs * rcs, tiny(1.0_dp))
      phase_rad = atan2(beta2, beta1)
      result%phase = phase_rad
      if (use_hours) result%phase = modulo(phase_rad * use_period / (2.0_dp * r_pi), use_period)
      result%amplitude = sqrt(beta1 * beta1 + beta2 * beta2)
      sy = finite_sd(y)
      result%sd = sy
      result%mean = finite_mean(y)
      result%intercept = cos(phase_rad) * result%amplitude * sy + result%mean
      fitted = cos(a - phase_rad)
      result%fit = finite_correlation(fitted, y)
   end subroutine cosinor_fit

   pure subroutine cosinor_period_fit(angle, y, periods, best_period, best_result, hours)
      real(dp), intent(in) :: angle(:) !! Measurement times in radians or clock units shared by all response columns.
      real(dp), intent(in) :: y(:, :) !! Observation-by-response matrix whose columns are fitted independently.
      real(dp), intent(in) :: periods(:) !! Candidate positive periods searched for the largest cosinor fit correlation.
      real(dp), allocatable, intent(out) :: best_period(:) !! Selected period for each response column.
      type(cosinor_result), allocatable, intent(out) :: best_result(:) !! Cosinor statistics at each selected period.
      logical, intent(in), optional :: hours !! Interpret `angle` as clock units when true; default false.
      type(cosinor_result) :: trial
      real(dp) :: best_fit
      integer :: j, k
      logical :: use_hours

      allocate(best_period(size(y, 2)), best_result(size(y, 2)))
      use_hours = .false.
      if (present(hours)) use_hours = hours
      do j = 1, size(y, 2)
         best_fit = -huge(1.0_dp)
         best_period(j) = 0.0_dp
         do k = 1, size(periods)
            call cosinor_fit(angle, y(:, j), trial, hours=use_hours, period=periods(k))
            if (trial%fit > best_fit) then
               best_fit = trial%fit
               best_period(j) = periods(k)
               best_result(j) = trial
            end if
         end do
      end do
   end subroutine cosinor_period_fit

   pure subroutine circadian_f_test(angle, group, f_value, df_between, df_within, hours, period)
      real(dp), intent(in) :: angle(:) !! Circular observations in radians or clock units.
      integer, intent(in) :: group(:) !! Positive integer group labels aligned with angle.
      real(dp), intent(out) :: f_value !! Circadian between/within F-ratio approximation.
      integer, intent(out) :: df_between !! Between-group degrees of freedom.
      integer, intent(out) :: df_within !! Within-group degrees of freedom.
      logical, intent(in), optional :: hours !! Interpret angles as clock units when true.
      real(dp), intent(in), optional :: period !! Clock period when hours is true; defaults to 24.
      type(circadian_result) :: all_stats, group_stats
      real(dp), allocatable :: subset(:)
      real(dp) :: sum_group_nr, ss_between, ss_within
      integer :: g, i, ng, nsub
      ng = 0
      if (size(group) > 0) ng = maxval(group)
      call circadian_stats_vector(angle, all_stats, hours, period)
      sum_group_nr = 0.0_dp
      do g = 1, ng
         nsub = count(group == g)
         allocate(subset(nsub))
         nsub = 0
         do i = 1, size(group)
            if (group(i) == g) then
               nsub = nsub + 1
               subset(nsub) = angle(i)
            end if
         end do
         call circadian_stats_vector(subset, group_stats, hours, period)
         sum_group_nr = sum_group_nr + real(group_stats%n, dp) * group_stats%resultant
         deallocate(subset)
      end do
      ss_between = sum_group_nr - real(all_stats%n, dp) * all_stats%resultant
      ss_within = real(all_stats%n, dp) - sum_group_nr
      df_between = max(ng - 1, 0)
      df_within = max(all_stats%n - ng, 0)
      if (df_between > 0 .and. df_within > 0 .and. ss_within > tiny(1.0_dp)) then
         f_value = (ss_between / real(df_between, dp)) / (ss_within / real(df_within, dp))
      else
         f_value = 0.0_dp
      end if
   end subroutine circadian_f_test

   pure subroutine circ_tests_values(loadings, gap_test, fisher_test, rotation_test, variance_test, sorting)
      real(dp), intent(in) :: loadings(:, :) !! Factor-loading matrix; the first two columns define circumplex locations.
      real(dp), intent(out) :: gap_test !! Variance of adjacent sorted angular gaps.
      real(dp), intent(out) :: fisher_test !! Coefficient of variation of loading radii.
      real(dp), intent(out) :: rotation_test !! Coefficient of variation of row-square-variance criterion across rotations.
      real(dp), intent(out) :: variance_test !! Coefficient of variation of first-axis square proportions across rotations.
      logical, intent(in), optional :: sorting !! Sort angles before gap differences; defaults to true.
      real(dp), allocatable :: theta(:), gaps(:), radius(:), criteria(:), rotated(:, :), rowvars(:)
      real(dp) :: angle, c, s, total
      logical :: do_sort
      integer :: i, k, n
      n = size(loadings, 1)
      allocate(theta(n), radius(n))
      do i = 1, n
         theta(i) = atan2(loadings(i, 2), loadings(i, 1))
         radius(i) = sqrt(loadings(i, 1) ** 2 + loadings(i, 2) ** 2)
      end do
      do_sort = .true.
      if (present(sorting)) do_sort = sorting
      if (do_sort) call sort_real(theta)
      if (n > 1) then
         allocate(gaps(n - 1))
         gaps = theta(2:n) - theta(1:n - 1)
         gap_test = finite_variance(gaps)
      else
         gap_test = 0.0_dp
      end if
      if (abs(finite_mean(radius)) > tiny(1.0_dp)) then
         fisher_test = finite_sd(radius) / abs(finite_mean(radius))
      else
         fisher_test = 0.0_dp
      end if
      allocate(criteria(10), rotated(n, 2), rowvars(n))
      do k = 0, 9
         angle = real(5 * k, dp) * r_pi / 180.0_dp
         c = cos(angle)
         s = sin(angle)
         rotated(:, 1) = c * loadings(:, 1) - s * loadings(:, 2)
         rotated(:, 2) = s * loadings(:, 1) + c * loadings(:, 2)
         do i = 1, n
            rowvars(i) = 0.5_dp * (rotated(i, 1) ** 2 - rotated(i, 2) ** 2) ** 2
         end do
         criteria(k + 1) = sum(rowvars)
      end do
      rotation_test = coefficient_of_variation(criteria)
      do k = 0, 9
         angle = real(5 * k, dp) * r_pi / 180.0_dp
         c = cos(angle)
         s = sin(angle)
         rotated(:, 1) = c * loadings(:, 1) - s * loadings(:, 2)
         rotated(:, 2) = s * loadings(:, 1) + c * loadings(:, 2)
         total = sum(rotated * rotated)
         if (total > tiny(1.0_dp)) then
            criteria(k + 1) = finite_variance(rotated(:, 1) ** 2 / total)
         else
            criteria(k + 1) = 0.0_dp
         end if
      end do
      variance_test = coefficient_of_variation(criteria)
   end subroutine circ_tests_values

   pure real(dp) function coefficient_of_variation(x) result(value)
      real(dp), intent(in) :: x(:) !! Numeric values whose standard-deviation-to-mean ratio is requested.
      real(dp) :: mu
      mu = finite_mean(x)
      if (abs(mu) > tiny(1.0_dp)) then
         value = finite_sd(x) / abs(mu)
      else
         value = 0.0_dp
      end if
   end function coefficient_of_variation

   pure subroutine sort_real(x)
      real(dp), intent(inout) :: x(:) !! Real vector sorted into ascending order in place.
      real(dp) :: key
      integer :: i, j
      do i = 2, size(x)
         key = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = key
      end do
   end subroutine sort_real
end module psych_circular
