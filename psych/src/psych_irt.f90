module psych_irt
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use psych_stats, only : finite_mean
   implicit none
   private
   public :: irt_probability, irt_information, irt_zero_parameter, rasch_item_difficulty
   public :: rasch_person_scores, two_parameter_person_scores, irt_discrimination
   public :: irt_expected_scores, irt_tau_from_thresholds
contains
   pure elemental real(dp) function irt_probability(theta, difficulty, discrimination) result(probability)
      real(dp), intent(in) :: theta !! Person latent-trait location.
      real(dp), intent(in) :: difficulty !! Item difficulty on the same latent scale.
      real(dp), intent(in) :: discrimination !! Item slope/discrimination parameter.
      real(dp) :: z
      z = discrimination * (theta - difficulty)
      if (z >= 0.0_dp) then
         probability = 1.0_dp / (1.0_dp + exp(-z))
      else
         probability = exp(z) / (1.0_dp + exp(z))
      end if
   end function irt_probability

   pure elemental real(dp) function irt_information(theta, difficulty, discrimination) result(information)
      real(dp), intent(in) :: theta !! Person latent-trait location.
      real(dp), intent(in) :: difficulty !! Item difficulty on the same latent scale.
      real(dp), intent(in) :: discrimination !! Item slope/discrimination parameter.
      real(dp) :: p
      p = irt_probability(theta, difficulty, discrimination)
      information = discrimination * discrimination * p * (1.0_dp - p)
   end function irt_information

   pure subroutine irt_zero_parameter(items, raw, theta, valid)
      real(dp), intent(in) :: items(:, :) !! Dichotomous subject-by-item response matrix; IEEE NaNs denote missing values.
      real(dp), allocatable, intent(out) :: raw(:) !! Mean observed item score for each subject.
      real(dp), allocatable, intent(out) :: theta(:) !! Logit transform after the same endpoint correction as upstream irt.0p.
      integer, allocatable, intent(out) :: valid(:) !! Number of nonmissing responses for each subject.
      real(dp) :: adjusted
      integer :: i, j, n_items
      n_items = size(items, 2)
      allocate(raw(size(items, 1)), theta(size(items, 1)), valid(size(items, 1)))
      do i = 1, size(items, 1)
         raw(i) = finite_mean(items(i, :))
         valid(i) = 0
         do j = 1, n_items
            if (.not. ieee_is_nan(items(i, j))) valid(i) = valid(i) + 1
         end do
         adjusted = raw(i)
         if (.not. ieee_is_nan(adjusted)) then
            if (adjusted < 0.0001_dp .and. n_items > 0) adjusted = 0.5_dp / real(n_items, dp)
            if (adjusted > 0.9999_dp .and. n_items > 0) adjusted = (real(n_items, dp) - 0.5_dp) / real(n_items, dp)
            theta(i) = log(adjusted / (1.0_dp - adjusted))
         else
            theta(i) = adjusted
         end if
      end do
   end subroutine irt_zero_parameter

   pure subroutine rasch_item_difficulty(items, difficulty)
      real(dp), intent(in) :: items(:, :) !! Dichotomous subject-by-item response matrix; IEEE NaNs denote missing values.
      real(dp), allocatable, intent(out) :: difficulty(:) !! Rasch log-odds difficulty estimate for each item.
      real(dp) :: p, lower
      integer :: j
      allocate(difficulty(size(items, 2)))
      if (size(items, 1) > 0) then
         lower = 1.0_dp / real(size(items, 1), dp)
      else
         lower = 0.5_dp
      end if
      do j = 1, size(items, 2)
         p = finite_mean(items(:, j))
         if (.not. ieee_is_nan(p)) then
            p = max(p, lower)
            p = min(p, 1.0_dp - epsilon(1.0_dp))
            difficulty(j) = log((1.0_dp - p) / p)
         else
            difficulty(j) = p
         end if
      end do
   end subroutine rasch_item_difficulty

   pure subroutine rasch_person_scores(difficulty, items, theta, fit, lower, upper)
      real(dp), intent(in) :: difficulty(:) !! Per-item Rasch difficulty parameters.
      real(dp), intent(in) :: items(:, :) !! Dichotomous subject-by-item response matrix.
      real(dp), allocatable, intent(out) :: theta(:) !! Maximum-likelihood trait estimate for each subject.
      real(dp), allocatable, intent(out) :: fit(:) !! Mean negative log likelihood at each optimum.
      real(dp), intent(in), optional :: lower !! Optional lower search bound; defaults to -4.
      real(dp), intent(in), optional :: upper !! Optional upper search bound; defaults to 4.
      call person_score_optimizer(difficulty, ones(size(difficulty)), items, theta, fit, lower, upper)
   end subroutine rasch_person_scores

   pure subroutine two_parameter_person_scores(difficulty, discrimination, items, theta, fit, lower, upper)
      real(dp), intent(in) :: difficulty(:) !! Per-item 2PL difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Per-item 2PL discrimination parameters.
      real(dp), intent(in) :: items(:, :) !! Dichotomous subject-by-item response matrix.
      real(dp), allocatable, intent(out) :: theta(:) !! Maximum-likelihood trait estimate for each subject.
      real(dp), allocatable, intent(out) :: fit(:) !! Mean negative log likelihood at each optimum.
      real(dp), intent(in), optional :: lower !! Optional lower search bound; defaults to -4.
      real(dp), intent(in), optional :: upper !! Optional upper search bound; defaults to 4.
      call person_score_optimizer(difficulty, discrimination, items, theta, fit, lower, upper)
   end subroutine two_parameter_person_scores

   pure subroutine irt_discrimination(item_difficulty, theta, items, discrimination, fit, lower, upper)
      real(dp), intent(in) :: item_difficulty(:) !! Difficulty parameter for each item.
      real(dp), intent(in) :: theta(:) !! Person trait estimates corresponding to response rows.
      real(dp), intent(in) :: items(:, :) !! Dichotomous subject-by-item response matrix.
      real(dp), allocatable, intent(out) :: discrimination(:) !! Maximum-likelihood discrimination estimate for each item.
      real(dp), allocatable, intent(out) :: fit(:) !! Mean negative log likelihood at each optimum.
      real(dp), intent(in), optional :: lower !! Optional lower slope bound; defaults to -5.
      real(dp), intent(in), optional :: upper !! Optional upper slope bound; defaults to 5.
      real(dp) :: lo, hi
      integer :: j
      lo = -5.0_dp
      hi = 5.0_dp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      allocate(discrimination(size(items, 2)), fit(size(items, 2)))
      do j = 1, size(items, 2)
         call golden_minimum_slope(lo, hi, item_difficulty(j), theta, items(:, j), discrimination(j), fit(j))
      end do
   end subroutine irt_discrimination

   pure subroutine irt_expected_scores(theta, difficulty, discrimination, probabilities, totals)
      real(dp), intent(in) :: theta(:) !! Person latent-trait locations.
      real(dp), intent(in) :: difficulty(:) !! Item difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Item discrimination parameters.
      real(dp), allocatable, intent(out) :: probabilities(:, :) !! Person-by-item expected dichotomous probabilities.
      real(dp), allocatable, intent(out) :: totals(:) !! Expected summed score for each person.
      integer :: i, j
      allocate(probabilities(size(theta), size(difficulty)), totals(size(theta)))
      do i = 1, size(theta)
         totals(i) = 0.0_dp
         do j = 1, size(difficulty)
            probabilities(i, j) = irt_probability(theta(i), difficulty(j), discrimination(j))
            totals(i) = totals(i) + probabilities(i, j)
         end do
      end do
   end subroutine irt_expected_scores

   pure subroutine irt_tau_from_thresholds(thresholds, tau)
      real(dp), intent(in) :: thresholds(:, :) !! Ordered item-category threshold matrix.
      real(dp), allocatable, intent(out) :: tau(:, :) !! Negated thresholds in psych's IRT tau convention.
      allocate(tau(size(thresholds, 1), size(thresholds, 2)))
      tau = -thresholds
   end subroutine irt_tau_from_thresholds

   pure subroutine person_score_optimizer(difficulty, discrimination, items, theta, fit, lower, upper)
      real(dp), intent(in) :: difficulty(:) !! Item difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Item discrimination parameters.
      real(dp), intent(in) :: items(:, :) !! Dichotomous subject-by-item response matrix.
      real(dp), allocatable, intent(out) :: theta(:) !! Optimizing trait estimate for each response row.
      real(dp), allocatable, intent(out) :: fit(:) !! Objective value for each response row.
      real(dp), intent(in), optional :: lower !! Optional lower optimization bound.
      real(dp), intent(in), optional :: upper !! Optional upper optimization bound.
      real(dp) :: lo, hi
      integer :: i
      lo = -4.0_dp
      hi = 4.0_dp
      if (present(lower)) lo = lower
      if (present(upper)) hi = upper
      allocate(theta(size(items, 1)), fit(size(items, 1)))
      do i = 1, size(items, 1)
         if (count_finite(items(i, :)) > 0) then
            call golden_minimum_theta(lo, hi, difficulty, discrimination, items(i, :), theta(i), fit(i))
         else
            theta(i) = ieee_value(0.0_dp, ieee_quiet_nan)
            fit(i) = theta(i)
         end if
      end do
   end subroutine person_score_optimizer

   pure subroutine golden_minimum_theta(lo, hi, difficulty, discrimination, response, location, objective)
      real(dp), intent(in) :: lo !! Lower scalar optimization bound.
      real(dp), intent(in) :: hi !! Upper scalar optimization bound.
      real(dp), intent(in) :: difficulty(:) !! Item difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Item discrimination parameters.
      real(dp), intent(in) :: response(:) !! One subject's dichotomous response vector.
      real(dp), intent(out) :: location !! Minimizing trait location.
      real(dp), intent(out) :: objective !! Minimum mean negative log likelihood.
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      real(dp) :: a, b, c, d, fc, fd
      integer :: iter
      a = lo
      b = hi
      c = b - ratio * (b - a)
      d = a + ratio * (b - a)
      fc = person_objective(c, difficulty, discrimination, response)
      fd = person_objective(d, difficulty, discrimination, response)
      do iter = 1, 100
         if (abs(b - a) <= 1.0e-8_dp) exit
         if (fc < fd) then
            b = d
            d = c
            fd = fc
            c = b - ratio * (b - a)
            fc = person_objective(c, difficulty, discrimination, response)
         else
            a = c
            c = d
            fc = fd
            d = a + ratio * (b - a)
            fd = person_objective(d, difficulty, discrimination, response)
         end if
      end do
      location = 0.5_dp * (a + b)
      objective = person_objective(location, difficulty, discrimination, response)
   end subroutine golden_minimum_theta

   pure subroutine golden_minimum_slope(lo, hi, difficulty, theta, response, slope, objective)
      real(dp), intent(in) :: lo !! Lower scalar optimization bound.
      real(dp), intent(in) :: hi !! Upper scalar optimization bound.
      real(dp), intent(in) :: difficulty !! Difficulty of the item being optimized.
      real(dp), intent(in) :: theta(:) !! Person trait estimates.
      real(dp), intent(in) :: response(:) !! Dichotomous responses for one item.
      real(dp), intent(out) :: slope !! Minimizing discrimination parameter.
      real(dp), intent(out) :: objective !! Minimum mean negative log likelihood.
      real(dp), parameter :: ratio = 0.6180339887498948482_dp
      real(dp) :: a, b, c, d, fc, fd
      integer :: iter
      a = lo
      b = hi
      c = b - ratio * (b - a)
      d = a + ratio * (b - a)
      fc = slope_objective(c, difficulty, theta, response)
      fd = slope_objective(d, difficulty, theta, response)
      do iter = 1, 100
         if (abs(b - a) <= 1.0e-8_dp) exit
         if (fc < fd) then
            b = d
            d = c
            fd = fc
            c = b - ratio * (b - a)
            fc = slope_objective(c, difficulty, theta, response)
         else
            a = c
            c = d
            fc = fd
            d = a + ratio * (b - a)
            fd = slope_objective(d, difficulty, theta, response)
         end if
      end do
      slope = 0.5_dp * (a + b)
      objective = slope_objective(slope, difficulty, theta, response)
   end subroutine golden_minimum_slope

   pure real(dp) function person_objective(theta, difficulty, discrimination, response) result(value)
      real(dp), intent(in) :: theta !! Candidate person trait location.
      real(dp), intent(in) :: difficulty(:) !! Item difficulty parameters.
      real(dp), intent(in) :: discrimination(:) !! Item discrimination parameters.
      real(dp), intent(in) :: response(:) !! Dichotomous response vector with IEEE NaNs ignored.
      real(dp) :: p
      integer :: j, n
      value = 0.0_dp
      n = 0
      do j = 1, size(response)
         if (ieee_is_nan(response(j))) cycle
         p = irt_probability(theta, difficulty(j), discrimination(j))
         p = min(max(p, tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
         value = value - response(j) * log(p) - (1.0_dp - response(j)) * log(1.0_dp - p)
         n = n + 1
      end do
      if (n > 0) value = value / real(n, dp)
   end function person_objective

   pure real(dp) function slope_objective(slope, difficulty, theta, response) result(value)
      real(dp), intent(in) :: slope !! Candidate item discrimination parameter.
      real(dp), intent(in) :: difficulty !! Fixed item difficulty.
      real(dp), intent(in) :: theta(:) !! Person trait locations.
      real(dp), intent(in) :: response(:) !! Item response vector with IEEE NaNs ignored.
      real(dp) :: p
      integer :: i, n
      value = 0.0_dp
      n = 0
      do i = 1, size(response)
         if (ieee_is_nan(response(i)) .or. ieee_is_nan(theta(i))) cycle
         p = irt_probability(theta(i), difficulty, slope)
         p = min(max(p, tiny(1.0_dp)), 1.0_dp - epsilon(1.0_dp))
         value = value - response(i) * log(p) - (1.0_dp - response(i)) * log(1.0_dp - p)
         n = n + 1
      end do
      if (n > 0) value = value / real(n, dp)
   end function slope_objective

   pure integer function count_finite(x) result(n)
      real(dp), intent(in) :: x(:) !! Vector whose non-NaN entries are counted.
      integer :: i
      n = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) n = n + 1
      end do
   end function count_finite

   pure function ones(n) result(x)
      integer, intent(in) :: n !! Number of unit discrimination parameters requested.
      real(dp) :: x(n)
      x = 1.0_dp
   end function ones
end module psych_irt
