module psych_irt_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use psych_distributions, only : normal_quantile
   use psych_irt, only : irt_information
   use psych_factor_extended, only : factor_to_irt_parameters
   implicit none
   private
   public :: irt_thresholds, irt_standard_errors, irt_response_curves
   public :: irt_stats_from_loadings, simulate_item_circumplex, irt_score_poly

contains

   pure subroutine irt_thresholds(items, minimum_category, maximum_category, thresholds)
      integer, intent(in) :: items(:, :) !! Subject-by-item integer responses; values outside the category range are ignored.
      integer, intent(in) :: minimum_category !! Lowest valid response category.
      integer, intent(in) :: maximum_category !! Highest valid response category.
      real(dp), allocatable, intent(out) :: thresholds(:, :) !! Item-by-threshold latent-normal cut points.
      integer :: i, j, k, n_valid, cumulative, n_thresholds
      real(dp) :: probability

      n_thresholds = max(0, maximum_category - minimum_category)
      allocate(thresholds(size(items, 2), n_thresholds))
      thresholds = 0.0_dp
      do j = 1, size(items, 2)
         n_valid = 0
         do i = 1, size(items, 1)
            if (items(i, j) >= minimum_category .and. items(i, j) <= maximum_category) n_valid = n_valid + 1
         end do
         if (n_valid == 0) cycle
         cumulative = 0
         do k = 1, n_thresholds
            do i = 1, size(items, 1)
               if (items(i, j) == minimum_category + k - 1) cumulative = cumulative + 1
            end do
            probability = real(cumulative, dp) / real(n_valid, dp)
            probability = max(1.0e-12_dp, min(1.0_dp - 1.0e-12_dp, probability))
            thresholds(j, k) = normal_quantile(probability)
         end do
      end do
   end subroutine irt_thresholds

   pure subroutine irt_standard_errors(scores, difficulty, discrimination, d_constant, information, standard_error)
      real(dp), intent(in) :: scores(:) !! Latent-trait locations at which information is evaluated.
      real(dp), intent(in) :: difficulty(:, :, :) !! Item-by-threshold-by-factor difficulty array.
      real(dp), intent(in) :: discrimination(:, :) !! Item-by-factor discrimination parameters.
      real(dp), intent(in) :: d_constant !! Link scaling constant, commonly 1.702 for logistic approximations.
      real(dp), allocatable, intent(out) :: information(:, :) !! Score-by-factor test information.
      real(dp), allocatable, intent(out) :: standard_error(:, :) !! Reciprocal square root of information.
      real(dp) :: a, b
      integer :: i, f, s, t

      allocate(information(size(scores), size(discrimination, 2)))
      allocate(standard_error(size(scores), size(discrimination, 2)))
      information = 0.0_dp
      do f = 1, size(discrimination, 2)
         do s = 1, size(scores)
            do i = 1, size(discrimination, 1)
               a = discrimination(i, f) * d_constant
               do t = 1, size(difficulty, 2)
                  b = difficulty(i, t, f)
                  information(s, f) = information(s, f) + irt_information(scores(s), b * sign(1.0_dp, a), a)
               end do
            end do
            if (information(s, f) > tiny(1.0_dp)) then
               standard_error(s, f) = sqrt(1.0_dp / information(s, f))
            else
               standard_error(s, f) = huge(1.0_dp)
            end if
         end do
      end do
   end subroutine irt_standard_errors

   pure subroutine irt_response_curves(theta, items, levels, n_bins, bin_centers, proportions, counts)
      real(dp), intent(in) :: theta(:) !! Person latent-trait estimates.
      integer, intent(in) :: items(:, :) !! Subject-by-item integer responses aligned with `theta`.
      integer, intent(in) :: levels(:) !! Response levels whose binwise proportions are required.
      integer, intent(in) :: n_bins !! Number of equal-width theta bins.
      real(dp), allocatable, intent(out) :: bin_centers(:) !! Midpoint of each theta bin.
      real(dp), allocatable, intent(out) :: proportions(:, :, :) !! Bin-by-item-by-level response proportions.
      integer, allocatable, intent(out) :: counts(:, :) !! Bin-by-item number of responses belonging to requested levels.
      real(dp) :: lo, hi, width
      integer :: b, i, j, lev, n_use

      allocate(bin_centers(max(1, n_bins)))
      allocate(proportions(max(1, n_bins), size(items, 2), size(levels)))
      allocate(counts(max(1, n_bins), size(items, 2)))
      proportions = 0.0_dp
      counts = 0
      if (size(theta) == 0) return
      lo = minval(theta)
      hi = maxval(theta)
      if (hi <= lo) hi = lo + 1.0_dp
      width = (hi - lo) / real(max(1, n_bins), dp)
      do b = 1, max(1, n_bins)
         bin_centers(b) = lo + (real(b, dp) - 0.5_dp) * width
      end do
      do i = 1, min(size(theta), size(items, 1))
         b = min(max(1, int((theta(i) - lo) / width) + 1), max(1, n_bins))
         do j = 1, size(items, 2)
            n_use = 0
            do lev = 1, size(levels)
               if (items(i, j) == levels(lev)) then
                  proportions(b, j, lev) = proportions(b, j, lev) + 1.0_dp
                  n_use = 1
                  exit
               end if
            end do
            counts(b, j) = counts(b, j) + n_use
         end do
      end do
      do b = 1, max(1, n_bins)
         do j = 1, size(items, 2)
            if (counts(b, j) > 0) proportions(b, j, :) = proportions(b, j, :) / real(counts(b, j), dp)
         end do
      end do
   end subroutine irt_response_curves

   pure subroutine irt_stats_from_loadings(items, minimum_category, maximum_category, loadings, thresholds, &
      discrimination, difficulty)
      integer, intent(in) :: items(:, :) !! Subject-by-item integer responses.
      integer, intent(in) :: minimum_category !! Lowest response category.
      integer, intent(in) :: maximum_category !! Highest response category.
      real(dp), intent(in) :: loadings(:, :) !! Item-by-factor loadings.
      real(dp), allocatable, intent(out) :: thresholds(:, :) !! Item-by-threshold latent-normal cut points.
      real(dp), allocatable, intent(out) :: discrimination(:, :) !! Item-by-factor IRT discriminations.
      real(dp), allocatable, intent(out) :: difficulty(:, :, :) !! Item-by-threshold-by-factor difficulties.

      call irt_thresholds(items, minimum_category, maximum_category, thresholds)
      call factor_to_irt_parameters(loadings, thresholds, discrimination, difficulty)
   end subroutine irt_stats_from_loadings

   subroutine simulate_item_circumplex(n_subjects, n_items, circumplex, x_loading, y_loading, g_loading, x_bias, y_bias, &
      categorical, low, high, truncate_below, cutpoint, items)
      integer, intent(in) :: n_subjects !! Number of simulated subjects.
      integer, intent(in) :: n_items !! Number of simulated items.
      logical, intent(in) :: circumplex !! Spread item directions uniformly around the circle when true.
      real(dp), intent(in) :: x_loading !! Loading scale for the x latent dimension.
      real(dp), intent(in) :: y_loading !! Loading scale for the y latent dimension.
      real(dp), intent(in) :: g_loading !! General-factor loading added to every item.
      real(dp), intent(in) :: x_bias !! Mean shift of the x latent variable.
      real(dp), intent(in) :: y_bias !! Mean shift of the y latent variable.
      logical, intent(in) :: categorical !! Round and bound generated responses when true.
      real(dp), intent(in) :: low !! Lower categorical bound.
      real(dp), intent(in) :: high !! Upper categorical bound.
      logical, intent(in) :: truncate_below !! Replace responses below `cutpoint` by zero when true.
      real(dp), intent(in) :: cutpoint !! Truncation threshold.
      real(dp), allocatable, intent(out) :: items(:, :) !! Simulated subject-by-item response matrix.
      real(dp), allocatable :: gx(:), x(:), y(:), error(:, :), angles(:)
      real(dp) :: av_loading, error_weight, pi_value
      integer :: i, j, quarter

      pi_value = acos(-1.0_dp)
      allocate(gx(n_subjects), x(n_subjects), y(n_subjects), error(n_subjects, n_items), angles(n_items))
      allocate(items(n_subjects, n_items))
      call fill_normal(gx)
      call fill_normal(x)
      call fill_normal(y)
      call fill_normal_matrix(error)
      x = x * x_loading + x_bias
      y = y * y_loading + y_bias
      if (circumplex) then
         do j = 1, n_items
            angles(j) = 2.0_dp * pi_value * real(j - 1, dp) / real(max(1, n_items), dp)
         end do
      else
         quarter = max(1, n_items / 4)
         do j = 1, n_items
            select case (min(3, (j - 1) / quarter))
            case (0)
               angles(j) = 0.0_dp
            case (1)
               angles(j) = 0.5_dp * pi_value
            case (2)
               angles(j) = pi_value
            case default
               angles(j) = 1.5_dp * pi_value
            end select
         end do
      end if
      av_loading = 0.5_dp * (x_loading + y_loading)
      error_weight = sqrt(max(0.0_dp, 1.0_dp - (av_loading * av_loading + g_loading * g_loading)))
      do j = 1, n_items
         do i = 1, n_subjects
            items(i, j) = g_loading * gx(i) + x(i) * cos(angles(j)) + y(i) * sin(angles(j)) + error_weight * error(i, j)
            if (categorical) items(i, j) = max(low, min(high, real(nint(items(i, j)), dp)))
            if (truncate_below .and. items(i, j) < cutpoint) items(i, j) = 0.0_dp
         end do
      end do
   end subroutine simulate_item_circumplex

   pure subroutine irt_score_poly(items, minimum_category, difficulty, discrimination, lower_bound, upper_bound, theta, fit)
      integer, intent(in) :: items(:, :) !! Subject-by-item integer ordinal responses.
      integer, intent(in) :: minimum_category !! Lowest item category represented by the thresholds.
      real(dp), intent(in) :: difficulty(:, :) !! Item-by-threshold difficulty matrix for one factor.
      real(dp), intent(in) :: discrimination(:) !! Item discrimination parameters.
      real(dp), intent(in) :: lower_bound !! Lower search bound for person theta.
      real(dp), intent(in) :: upper_bound !! Upper search bound for person theta.
      real(dp), allocatable, intent(out) :: theta(:) !! Maximum-likelihood person scores.
      real(dp), allocatable, intent(out) :: fit(:) !! Mean negative log likelihood at each optimum.
      real(dp) :: lo, hi, x1, x2, f1, f2
      integer :: i, iter

      allocate(theta(size(items, 1)), fit(size(items, 1)))
      do i = 1, size(items, 1)
         lo = lower_bound
         hi = upper_bound
         do iter = 1, 80
            x1 = lo + (hi - lo) / 3.0_dp
            x2 = hi - (hi - lo) / 3.0_dp
            f1 = ordinal_negloglik(x1, items(i, :), minimum_category, difficulty, discrimination)
            f2 = ordinal_negloglik(x2, items(i, :), minimum_category, difficulty, discrimination)
            if (f1 <= f2) then
               hi = x2
            else
               lo = x1
            end if
         end do
         theta(i) = 0.5_dp * (lo + hi)
         fit(i) = ordinal_negloglik(theta(i), items(i, :), minimum_category, difficulty, discrimination)
      end do
   end subroutine irt_score_poly

   pure real(dp) function ordinal_negloglik(theta, response, minimum_category, difficulty, discrimination) result(value)
      real(dp), intent(in) :: theta !! Candidate person trait value.
      integer, intent(in) :: response(:) !! One subject's item responses.
      integer, intent(in) :: minimum_category !! Lowest category represented by zero thresholds below it.
      real(dp), intent(in) :: difficulty(:, :) !! Item-by-threshold difficulty matrix.
      real(dp), intent(in) :: discrimination(:) !! Item discrimination vector.
      real(dp) :: lower_tail, upper_tail, probability, a
      integer :: i, category, n_used, n_thresholds

      value = 0.0_dp
      n_used = 0
      n_thresholds = size(difficulty, 2)
      do i = 1, min(size(response), size(discrimination))
         category = response(i) - minimum_category
         if (category < 0 .or. category > n_thresholds) cycle
         a = discrimination(i)
         if (category == 0) then
            probability = 1.0_dp - logistic_cdf(theta, difficulty(i, 1), a)
         else if (category == n_thresholds) then
            probability = logistic_cdf(theta, difficulty(i, n_thresholds), a)
         else
            lower_tail = logistic_cdf(theta, difficulty(i, category), a)
            upper_tail = logistic_cdf(theta, difficulty(i, category + 1), a)
            probability = lower_tail - upper_tail
         end if
         probability = max(1.0e-12_dp, probability)
         value = value - log(probability)
         n_used = n_used + 1
      end do
      if (n_used > 0) value = value / real(n_used, dp)
   end function ordinal_negloglik

   pure real(dp) function logistic_cdf(theta, difficulty, discrimination) result(value)
      real(dp), intent(in) :: theta !! Person latent trait.
      real(dp), intent(in) :: difficulty !! Threshold difficulty.
      real(dp), intent(in) :: discrimination !! Item discrimination.
      real(dp) :: z

      z = discrimination * (theta - difficulty)
      if (z >= 0.0_dp) then
         value = 1.0_dp / (1.0_dp + exp(-z))
      else
         value = exp(z) / (1.0_dp + exp(z))
      end if
   end function logistic_cdf

   subroutine fill_normal(x)
      real(dp), intent(out) :: x(:) !! Vector filled with independent standard-normal draws.
      real(dp) :: u1, u2, radius, angle
      integer :: i

      i = 1
      do while (i <= size(x))
         call random_number(u1)
         call random_number(u2)
         u1 = max(u1, tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u1))
         angle = 2.0_dp * acos(-1.0_dp) * u2
         x(i) = radius * cos(angle)
         if (i + 1 <= size(x)) x(i + 1) = radius * sin(angle)
         i = i + 2
      end do
   end subroutine fill_normal

   subroutine fill_normal_matrix(x)
      real(dp), intent(out) :: x(:, :) !! Matrix filled with independent standard-normal draws.
      real(dp), allocatable :: work(:)
      integer :: j

      allocate(work(size(x, 1)))
      do j = 1, size(x, 2)
         call fill_normal(work)
         x(:, j) = work
      end do
   end subroutine fill_normal_matrix

end module psych_irt_extended
