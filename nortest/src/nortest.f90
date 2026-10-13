! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of nortest 1.0-4 by Juergen Gross and Uwe Ligges.
module nortest
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan
   use r_kinds, only: dp
   use r_descriptive, only: r_mean, r_sd, r_correlation
   use r_distributions, only: r_pnorm, r_qnorm, r_pchisq
   use r_ordering, only: r_sort_values_in_place
   implicit none
   private

   integer, parameter, public :: normality_success = 0
   integer, parameter, public :: normality_invalid_size = 1
   integer, parameter, public :: normality_nonfinite = 2
   integer, parameter, public :: normality_zero_variance = 3
   integer, parameter, public :: normality_invalid_classes = 4

   type, public :: normality_result
      !! Numerical result and validation status, replacing R's htest object.
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 0.0_dp
      integer :: n = 0
      integer :: status = normality_success
      integer :: n_classes = 0
      integer :: df = 0
      logical :: p_value_at_floor = .false.
   end type normality_result

   public :: dp, ad_test, cvm_test, lillie_test, pearson_test, sf_test

contains

   pure function ad_test(x) result(test)
      !! Anderson-Darling test for normality with estimated mean and variance.
      real(dp), intent(in) :: x(:) !! Sample vector, with NaNs omitted and at least eight nonmissing observations.
      type(normality_result) :: test
      real(dp), allocatable :: z(:)
      real(dp) :: nr, total, adjusted, weight
      integer :: i

      call prepare_sample(x, 8, huge(1), z, test)
      if (test%status /= normality_success) return
      nr = real(test%n, dp)
      total = 0.0_dp
      do i = 1, test%n
         weight = 2.0_dp * real(i, dp) - 1.0_dp
         total = total + weight * (r_pnorm(z(i), log_probability = .true.) + &
                                  r_pnorm(-z(test%n - i + 1), log_probability = .true.))
      end do
      test%statistic = -nr - total / nr
      adjusted = (1.0_dp + 0.75_dp / nr + 2.25_dp / nr ** 2) * test%statistic
      if (adjusted < 0.2_dp) then
         test%p_value = 1.0_dp - exp(-13.436_dp + 101.14_dp * adjusted - 223.73_dp * adjusted ** 2)
      else if (adjusted < 0.34_dp) then
         test%p_value = 1.0_dp - exp(-8.318_dp + 42.796_dp * adjusted - 59.938_dp * adjusted ** 2)
      else if (adjusted < 0.6_dp) then
         test%p_value = exp(0.9177_dp - 4.279_dp * adjusted - 1.38_dp * adjusted ** 2)
      else if (adjusted < 10.0_dp) then
         test%p_value = exp(1.2937_dp - 5.709_dp * adjusted + 0.0186_dp * adjusted ** 2)
      else
         test%p_value = 3.7e-24_dp
         test%p_value_at_floor = .true.
      end if
   end function ad_test

   pure function cvm_test(x) result(test)
      !! Cramer-von Mises test with the finite-sample normality calibration.
      real(dp), intent(in) :: x(:) !! Sample vector, with NaNs omitted and at least eight nonmissing observations.
      type(normality_result) :: test
      real(dp), allocatable :: z(:)
      real(dp) :: nr, adjusted, difference
      integer :: i

      call prepare_sample(x, 8, huge(1), z, test)
      if (test%status /= normality_success) return
      nr = real(test%n, dp)
      test%statistic = 1.0_dp / (12.0_dp * nr)
      do i = 1, test%n
         difference = r_pnorm(z(i)) - (2.0_dp * real(i, dp) - 1.0_dp) / (2.0_dp * nr)
         test%statistic = test%statistic + difference ** 2
      end do
      adjusted = (1.0_dp + 0.5_dp / nr) * test%statistic
      if (adjusted < 0.0275_dp) then
         test%p_value = 1.0_dp - exp(-13.953_dp + 775.5_dp * adjusted - 12542.61_dp * adjusted ** 2)
      else if (adjusted < 0.051_dp) then
         test%p_value = 1.0_dp - exp(-5.903_dp + 179.546_dp * adjusted - 1515.29_dp * adjusted ** 2)
      else if (adjusted < 0.092_dp) then
         test%p_value = exp(0.886_dp - 31.62_dp * adjusted + 10.897_dp * adjusted ** 2)
      else if (adjusted < 1.1_dp) then
         test%p_value = exp(1.111_dp - 34.242_dp * adjusted + 12.832_dp * adjusted ** 2)
      else
         test%p_value = 7.37e-10_dp
         test%p_value_at_floor = .true.
      end if
   end function cvm_test

   pure function lillie_test(x) result(test)
      !! Lilliefors normality test with Dallal-Wilkinson and Stephens approximations.
      real(dp), intent(in) :: x(:) !! Sample vector, with NaNs omitted and at least five nonmissing observations.
      type(normality_result) :: test
      real(dp), allocatable :: z(:)
      real(dp) :: nr, probability, kd, nd, kk
      integer :: i

      call prepare_sample(x, 5, huge(1), z, test)
      if (test%status /= normality_success) return
      nr = real(test%n, dp)
      test%statistic = 0.0_dp
      do i = 1, test%n
         probability = r_pnorm(z(i))
         test%statistic = max(test%statistic, real(i, dp) / nr - probability, &
                             probability - real(i - 1, dp) / nr)
      end do
      kd = test%statistic
      nd = nr
      if (test%n > 100) then
         kd = kd * (nr / 100.0_dp) ** 0.49_dp
         nd = 100.0_dp
      end if
      test%p_value = exp(-7.01256_dp * kd ** 2 * (nd + 2.78019_dp) + &
                        2.99587_dp * kd * sqrt(nd + 2.78019_dp) - 0.122119_dp + &
                        0.974598_dp / sqrt(nd) + 1.67997_dp / nd)
      if (test%p_value > 0.1_dp) then
         kk = (sqrt(nr) - 0.01_dp + 0.85_dp / sqrt(nr)) * test%statistic
         if (kk <= 0.302_dp) then
            test%p_value = 1.0_dp
         else if (kk <= 0.5_dp) then
            test%p_value = 2.76773_dp - 19.828315_dp * kk + 80.709644_dp * kk ** 2 - &
                           138.55152_dp * kk ** 3 + 81.218052_dp * kk ** 4
         else if (kk <= 0.9_dp) then
            test%p_value = -4.901232_dp + 40.662806_dp * kk - 97.490286_dp * kk ** 2 + &
                           94.029866_dp * kk ** 3 - 32.355711_dp * kk ** 4
         else if (kk <= 1.31_dp) then
            test%p_value = 6.198765_dp - 19.558097_dp * kk + 23.186922_dp * kk ** 2 - &
                           12.234627_dp * kk ** 3 + 2.423045_dp * kk ** 4
         else
            test%p_value = 0.0_dp
         end if
      end if
   end function lillie_test

   pure function pearson_test(x, n_classes, adjust) result(test)
      !! Pearson equal-probability-class chi-square test for fitted normality.
      real(dp), intent(in) :: x(:) !! Sample vector, with NaNs omitted and at least two nonmissing observations.
      integer, intent(in), optional :: n_classes !! Positive number of classes, default ceiling(2 * n ** 0.4).
      logical, intent(in), optional :: adjust !! Subtract two fitted parameters from degrees of freedom, default true.
      type(normality_result) :: test
      real(dp), allocatable :: z(:)
      real(dp) :: nr, probability, classes, expected, total, current_count, occupied
      integer :: i, bin, previous_bin
      logical :: fitted_adjustment

      call prepare_sample(x, 2, huge(1), z, test)
      if (test%status /= normality_success) return
      nr = real(test%n, dp)
      test%n_classes = ceiling(2.0_dp * nr ** 0.4_dp)
      if (present(n_classes)) test%n_classes = n_classes
      fitted_adjustment = .true.
      if (present(adjust)) fitted_adjustment = adjust
      if (test%n_classes < 1) then
         test%status = normality_invalid_classes
         return
      end if
      test%df = test%n_classes - 1
      if (fitted_adjustment) test%df = test%df - 2
      if (test%df < 0) then
         test%status = normality_invalid_classes
         return
      end if

      ! Count only occupied bins so memory does not grow with n_classes.
      ! As in R tabulate(), a rounded CDF of one produces a discarded bin.
      classes = real(test%n_classes, dp)
      expected = nr / classes
      total = 0.0_dp
      current_count = 0.0_dp
      occupied = 0.0_dp
      previous_bin = 0
      do i = 1, test%n
         probability = r_pnorm(z(i))
         if (probability * classes >= classes) cycle
         bin = floor(1.0_dp + classes * probability)
         if (bin > test%n_classes) cycle
         if (bin /= previous_bin) then
            if (previous_bin > 0) then
               total = total + (current_count - expected) ** 2 / expected
               occupied = occupied + 1.0_dp
            end if
            current_count = 0.0_dp
            previous_bin = bin
         end if
         current_count = current_count + 1.0_dp
      end do
      if (previous_bin > 0) then
         total = total + (current_count - expected) ** 2 / expected
         occupied = occupied + 1.0_dp
      end if
      test%statistic = total + (classes - occupied) * expected
      if (test%df == 0) then
         test%p_value = 0.0_dp
         if (test%statistic == 0.0_dp) test%p_value = 1.0_dp
      else
         test%p_value = r_pchisq(test%statistic, real(test%df, dp), lower_tail = .false.)
      end if
   end function pearson_test

   pure function sf_test(x) result(test)
      !! Shapiro-Francia test using normal order-statistic quantiles and Royston's calibration.
      real(dp), intent(in) :: x(:) !! Sample vector, with NaNs omitted and five to 5000 nonmissing observations.
      type(normality_result) :: test
      real(dp), allocatable :: z(:), y(:)
      real(dp) :: nr, u, v, mu, sigma, transformed
      integer :: i

      call prepare_sample(x, 5, 5000, z, test)
      if (test%status /= normality_success) return
      nr = real(test%n, dp)
      allocate(y(test%n))
      do i = 1, test%n
         y(i) = r_qnorm((real(i, dp) - 0.375_dp) / (nr + 0.25_dp))
      end do
      test%statistic = min(1.0_dp, max(0.0_dp, r_correlation(z, y) ** 2))
      if (test%statistic == 1.0_dp) then
         test%p_value = 1.0_dp
         return
      end if
      u = log(nr)
      v = log(u)
      mu = -1.2725_dp + 1.0521_dp * (v - u)
      sigma = 1.0308_dp - 0.26758_dp * (v + 2.0_dp / u)
      transformed = (log(1.0_dp - test%statistic) - mu) / sigma
      test%p_value = r_pnorm(transformed, lower_tail = .false.)
   end function sf_test

   pure subroutine prepare_sample(x, minimum_n, maximum_n, z, test)
      !! Remove missing values, validate, sort, and standardize with sample standard deviation.
      real(dp), intent(in) :: x(:) !! Original sample, with IEEE NaNs representing R NA or NaN observations.
      integer, intent(in) :: minimum_n !! Smallest accepted number of nonmissing observations for the caller.
      integer, intent(in) :: maximum_n !! Largest accepted number of nonmissing observations for the caller.
      real(dp), allocatable, intent(out) :: z(:) !! Sorted standardized data on success, otherwise undefined contents.
      type(normality_result), intent(out) :: test !! Initialized result with effective sample size and validation status.
      real(dp) :: center, sigma, magnitude

      test%statistic = ieee_value(0.0_dp, ieee_quiet_nan)
      test%p_value = ieee_value(0.0_dp, ieee_quiet_nan)
      z = pack(x, .not. ieee_is_nan(x))
      test%n = size(z)
      if (test%n < minimum_n .or. test%n > maximum_n) then
         test%status = normality_invalid_size
         return
      end if
      if (any(.not. ieee_is_finite(z))) then
         test%status = normality_nonfinite
         return
      end if
      magnitude = maxval(abs(z))
      if (magnitude == 0.0_dp) then
         test%status = normality_zero_variance
         return
      end if
      ! Avoid overflow/underflow while retaining original-scale accuracy normally.
      if (magnitude > sqrt(huge(1.0_dp)) / real(test%n, dp) .or. magnitude < sqrt(tiny(1.0_dp))) then
         z = z / magnitude
      end if
      call r_sort_values_in_place(z)
      if (z(1) == z(test%n)) then
         test%status = normality_zero_variance
         return
      end if
      center = r_mean(z)
      sigma = r_sd(z, center = center)
      if (sigma <= 0.0_dp .or. .not. ieee_is_finite(sigma)) then
         test%status = normality_zero_variance
         return
      end if
      z = (z - center) / sigma
   end subroutine prepare_sample

end module nortest
