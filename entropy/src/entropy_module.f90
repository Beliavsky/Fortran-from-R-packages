! SPDX-License-Identifier: GPL-3.0-or-later
! Derived from entropy 1.3.2, Copyright 2008, 2008-13, 2013, 2021 Korbinian Strimmer.
! Fortran translation Copyright 2026 entropy Fortran contributors.
! This program comes without any warranty. See LICENSE and NOTICE.md.
module entropy_module
   use r_kinds, only: dp
   use r_distributions, only: r_pchisq
   use r_transforms, only: r_log1p, r_expm1
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
                                                                               ieee_quiet_nan, ieee_positive_inf
   implicit none
   private
   public :: dp, entropy, freqs, entropy_plugin, entropy_empirical, entropy_dirichlet, entropy_shrink
   public :: entropy_millermadow, entropy_chaoshen, freqs_empirical, freqs_dirichlet, freqs_shrink
   public :: kl_plugin, kl_empirical, kl_dirichlet, kl_shrink
   public :: chi2_plugin, chi2_empirical, chi2_dirichlet, chi2_shrink
   public :: mi_plugin, mi_empirical, mi_dirichlet, mi_shrink
   public :: chi2indep_plugin, chi2indep_empirical, chi2indep_dirichlet, chi2indep_shrink
   public :: gstat, chi2stat, gstatindep, chi2statindep, test_result
   public :: discretize, discretize2d, shrinkage_intensity

   type :: test_result
      real(dp) :: stat = 0.0_dp !! Statistic in the requested logarithm unit.
      integer :: df = 0 !! Degrees of freedom, including bins whose observed count is zero.
      real(dp) :: pval = 0.0_dp !! Upper-tail chi-squared probability, or NaN for invalid input.
   end type test_result

   interface freqs_empirical
      module procedure freqs_empirical_1d, freqs_empirical_2d
   end interface
   interface freqs_dirichlet
      module procedure freqs_dirichlet_scalar, freqs_dirichlet_vector
      module procedure freqs_dirichlet_matrix, freqs_dirichlet_matrix_prior
   end interface
   interface freqs_shrink
      module procedure freqs_shrink_1d, freqs_shrink_2d
   end interface

contains

   pure function missing() result(value)
      real(dp) :: value
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function missing

   pure elemental function unit_divisor(unit) result(value)
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = 1.0_dp
      if (.not. present(unit)) return
      select case (unit)
      case ('log')
         value = 1.0_dp
      case ('log2')
         value = log(2.0_dp)
      case ('log10')
         value = log(10.0_dp)
      case default
         value = missing()
      end select
   end function unit_divisor

   pure function valid_counts(y) result(valid)
      real(dp), intent(in) :: y(:) !! Nonempty vector of finite, nonnegative weights or counts.
      logical :: valid
      valid = size(y) > 0
      if (.not. valid) return
      valid = all(ieee_is_finite(y))
      if (.not. valid) return
      valid = all(y >= 0.0_dp)
   end function valid_counts

   pure function count_total(y) result(total)
      real(dp), intent(in) :: y(:) !! Finite nonnegative counts; overflowing totals return NaN.
      real(dp) :: total, scale, scaled_total
      total = missing()
      if (.not. valid_counts(y)) return
      scale = maxval(y)
      total = 0.0_dp
      if (scale == 0.0_dp) return
      scaled_total = sum(y / scale)
      if (scale > huge(1.0_dp) / scaled_total) then
         total = missing()
      else
         total = scale * scaled_total
      end if
   end function count_total

   pure function freqs_empirical_1d(y) result(prob)
      real(dp), intent(in) :: y(:) !! Nonnegative bin counts or weights; zero total gives all NaNs.
      real(dp) :: prob(size(y)), scale
      prob = missing()
      if (.not. valid_counts(y)) return
      scale = maxval(y)
      if (scale == 0.0_dp) return
      prob = y / scale
      prob = prob / sum(prob)
   end function freqs_empirical_1d

   pure function freqs_empirical_2d(y) result(prob)
      real(dp), intent(in) :: y(:, :) !! Nonnegative contingency-table counts, normalized over all cells.
      real(dp) :: prob(size(y, 1), size(y, 2))
      prob = reshape(freqs_empirical_1d(reshape(y, [size(y)])), shape(y))
   end function freqs_empirical_2d

   pure function freqs_dirichlet_scalar(y, a) result(prob)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts, including unobserved bins.
      real(dp), intent(in) :: a !! Nonnegative finite pseudocount applied to each bin.
      real(dp) :: prob(size(y))
      prob = freqs_dirichlet_vector(y, spread(a, 1, size(y)))
   end function freqs_dirichlet_scalar

   pure function freqs_dirichlet_vector(y, a) result(prob)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts, including unobserved bins.
      real(dp), intent(in) :: a(:) !! Nonnegative finite pseudocount for every bin; must match y's size.
      real(dp) :: prob(size(y)), scale
      prob = missing()
      if (size(y) /= size(a)) return
      if (.not. valid_counts(y) .or. .not. valid_counts(a)) return
      scale = max(maxval(y), maxval(a))
      if (scale == 0.0_dp) return
      prob = freqs_empirical_1d(y / scale + a / scale)
   end function freqs_dirichlet_vector

   pure function freqs_dirichlet_matrix(y, a) result(prob)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts.
      real(dp), intent(in) :: a !! Nonnegative finite pseudocount applied to each cell.
      real(dp) :: prob(size(y, 1), size(y, 2))
      prob = reshape(freqs_dirichlet_scalar(reshape(y, [size(y)]), a), shape(y))
   end function freqs_dirichlet_matrix

   pure function freqs_dirichlet_matrix_prior(y, a) result(prob)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts.
      real(dp), intent(in) :: a(:, :) !! Cell-specific nonnegative pseudocounts with the same shape as y.
      real(dp) :: prob(size(y, 1), size(y, 2))
      prob = missing()
      if (any(shape(y) /= shape(a))) return
      prob = reshape(freqs_dirichlet_vector(reshape(y, [size(y)]), reshape(a, [size(a)])), shape(y))
   end function freqs_dirichlet_matrix_prior

   pure function shrinkage_intensity(y) result(lambda)
      real(dp), intent(in) :: y(:) !! Nonnegative bin counts; automatic uniform-target shrinkage requires positive total.
      real(dp) :: lambda, n, target, msp, prob(size(y))
      lambda = missing()
      n = count_total(y)
      if (.not. ieee_is_finite(n)) return
      if (n == 0.0_dp .or. n == 1.0_dp) then
         lambda = 1.0_dp
         return
      end if
      prob = freqs_empirical(y)
      target = 1.0_dp / real(size(y), dp)
      msp = sum((prob - target) ** 2)
      lambda = 1.0_dp
      if (msp > 0.0_dp) lambda = max(0.0_dp, min(1.0_dp, sum(prob * (1.0_dp - prob)) / (n - 1.0_dp) / msp))
   end function shrinkage_intensity

   pure function freqs_shrink_1d(y, lambda) result(prob)
      real(dp), intent(in) :: y(:) !! Nonnegative finite counts; an all-zero sample produces NaNs, as in R.
      real(dp), intent(in), optional :: lambda !! Fixed shrinkage intensity in [0,1]; absent means estimate from counts.
      real(dp) :: prob(size(y)), intensity
      prob = missing()
      if (.not. valid_counts(y)) return
      if (present(lambda)) then
         intensity = lambda
      else
         intensity = shrinkage_intensity(y)
      end if
      if (.not. ieee_is_finite(intensity)) return
      if (intensity < 0.0_dp .or. intensity > 1.0_dp) return
      prob = freqs_empirical(y)
      if (any(ieee_is_nan(prob))) return
      prob = intensity / real(size(y), dp) + (1.0_dp - intensity) * prob
   end function freqs_shrink_1d

   pure function freqs_shrink_2d(y, lambda) result(prob)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite table counts; shrinkage target is uniform over all cells.
      real(dp), intent(in), optional :: lambda !! Fixed shrinkage intensity in [0,1]; absent means automatic estimation.
      real(dp) :: prob(size(y, 1), size(y, 2))
      prob = reshape(freqs_shrink_1d(reshape(y, [size(y)]), lambda), shape(y))
   end function freqs_shrink_2d

   pure function entropy_plugin(weights, unit) result(value)
      real(dp), intent(in) :: weights(:) !! Nonnegative finite weights, internally normalized to probabilities.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value, prob(size(weights))
      integer :: i
      value = missing()
      if (.not. valid_counts(weights)) return
      prob = freqs_empirical(weights)
      if (any(ieee_is_nan(prob))) return
      value = 0.0_dp
      do i = 1, size(prob)
         if (prob(i) > 0.0_dp) value = value - prob(i) * log(prob(i))
      end do
      value = value / unit_divisor(unit)
   end function entropy_plugin

   pure function kl_plugin(weights1, weights2, unit) result(value)
      real(dp), intent(in) :: weights1(:) !! Nonnegative numerator-distribution weights with positive total.
      real(dp), intent(in) :: weights2(:) !! Same-length reference weights; zero reference at positive numerator gives infinity.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value, p(size(weights1)), q(size(weights2))
      integer :: i
      value = missing()
      if (size(weights1) /= size(weights2)) return
      if (.not. valid_counts(weights1) .or. .not. valid_counts(weights2)) return
      p = freqs_empirical(weights1)
      q = freqs_empirical(weights2)
      if (any(ieee_is_nan(p)) .or. any(ieee_is_nan(q))) return
      value = 0.0_dp
      do i = 1, size(p)
         if (p(i) == 0.0_dp) cycle
         if (q(i) == 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_positive_inf)
            exit
         end if
         value = value + p(i) * (log(p(i)) - log(q(i)))
      end do
      value = value / unit_divisor(unit)
   end function kl_plugin

   pure function chi2_plugin(weights1, weights2, unit) result(value)
      real(dp), intent(in) :: weights1(:) !! Nonnegative numerator-distribution weights with positive total.
      real(dp), intent(in) :: weights2(:) !! Same-length reference weights; common zero cells give NaN, matching R.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value, p(size(weights1)), q(size(weights2))
      integer :: i
      value = missing()
      if (size(weights1) /= size(weights2)) return
      if (.not. valid_counts(weights1) .or. .not. valid_counts(weights2)) return
      p = freqs_empirical(weights1)
      q = freqs_empirical(weights2)
      if (any(ieee_is_nan(p)) .or. any(ieee_is_nan(q))) return
      if (any(p == 0.0_dp .and. q == 0.0_dp)) return
      value = 0.0_dp
      do i = 1, size(p)
         if (q(i) == 0.0_dp) then
            value = ieee_value(0.0_dp, ieee_positive_inf)
            exit
         end if
         value = value + ((p(i) - q(i)) / sqrt(q(i))) ** 2
      end do
      value = value / unit_divisor(unit)
   end function chi2_plugin

   pure function mi_plugin(weights, unit) result(value)
      real(dp), intent(in) :: weights(:, :) !! Nonnegative finite joint weights, with rows and columns denoting two variables.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value, p(size(weights, 1), size(weights, 2)), null(size(weights, 1), size(weights, 2))
      p = freqs_empirical(weights)
      null = spread(sum(p, dim=2), 2, size(p, 2)) * spread(sum(p, dim=1), 1, size(p, 1))
      value = kl_plugin(reshape(p, [size(p)]), reshape(null, [size(null)]), unit)
   end function mi_plugin

   pure function chi2indep_plugin(weights, unit) result(value)
      real(dp), intent(in) :: weights(:, :) !! Nonnegative finite joint weights; zero row or column margins give NaN, as in R.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value, p(size(weights, 1), size(weights, 2)), null(size(weights, 1), size(weights, 2))
      p = freqs_empirical(weights)
      null = spread(sum(p, dim=2), 2, size(p, 2)) * spread(sum(p, dim=1), 1, size(p, 1))
      value = chi2_plugin(reshape(p, [size(p)]), reshape(null, [size(null)]), unit)
   end function chi2indep_plugin

   pure function entropy_millermadow(y, unit) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts; zero or overflowing total gives NaN.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value, n
      value = missing()
      n = count_total(y)
      if (.not. ieee_is_finite(n)) return
      if (n <= 0.0_dp) return
      value = (entropy_plugin(y) + 0.5_dp * real(count(y > 0.0_dp) - 1, dp) / n) / unit_divisor(unit)
   end function entropy_millermadow

   pure function entropy_chaoshen(y, unit) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts; zeros are omitted and count-one bins are singletons.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value, n, singletons, coverage, adjusted, detection
      integer :: i
      value = missing()
      n = count_total(y)
      if (.not. ieee_is_finite(n)) return
      if (n <= 0.0_dp) return
      singletons = real(count(y == 1.0_dp), dp)
      if (singletons == n) singletons = n - 1.0_dp
      coverage = 1.0_dp - singletons / n
      value = 0.0_dp
      do i = 1, size(y)
         if (y(i) <= 0.0_dp) cycle
         adjusted = coverage * (y(i) / n)
         if (adjusted == 1.0_dp) cycle
         if (adjusted <= 0.0_dp) then
            value = missing()
            return
         end if
         detection = -r_expm1(n * r_log1p(-adjusted))
         value = value - adjusted * log(adjusted) / detection
      end do
      value = value / unit_divisor(unit)
   end function entropy_chaoshen

   pure function entropy_empirical(y, unit) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts with positive total.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = entropy_plugin(freqs_empirical(y), unit)
   end function entropy_empirical

   pure function entropy_dirichlet(y, a, unit) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: a !! Nonnegative finite pseudocount added to every cell of the corresponding sample.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = entropy_plugin(freqs_dirichlet(y, a), unit)
   end function entropy_dirichlet

   pure function entropy_shrink(y, lambda, unit) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in), optional :: lambda !! Shrinkage intensity in [0,1] for the corresponding sample; absent means estimate.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = entropy_plugin(freqs_shrink(y, lambda), unit)
   end function entropy_shrink

   pure function mi_empirical(y, unit) result(value)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = mi_plugin(freqs_empirical(y), unit)
   end function mi_empirical

   pure function mi_dirichlet(y, a, unit) result(value)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      real(dp), intent(in) :: a !! Nonnegative finite pseudocount added to every cell of the corresponding sample.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = mi_plugin(freqs_dirichlet(y, a), unit)
   end function mi_dirichlet

   pure function mi_shrink(y, lambda, unit) result(value)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      real(dp), intent(in), optional :: lambda !! Shrinkage intensity in [0,1] for the corresponding sample; absent means estimate.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = mi_plugin(freqs_shrink(y, lambda), unit)
   end function mi_shrink

   pure function chi2indep_empirical(y, unit) result(value)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = chi2indep_plugin(freqs_empirical(y), unit)
   end function chi2indep_empirical

   pure function chi2indep_dirichlet(y, a, unit) result(value)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      real(dp), intent(in) :: a !! Nonnegative finite pseudocount added to every cell of the corresponding sample.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = chi2indep_plugin(freqs_dirichlet(y, a), unit)
   end function chi2indep_dirichlet

   pure function chi2indep_shrink(y, lambda, unit) result(value)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      real(dp), intent(in), optional :: lambda !! Shrinkage intensity in [0,1] for the corresponding sample; absent means estimate.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = chi2indep_plugin(freqs_shrink(y, lambda), unit)
   end function chi2indep_shrink

   pure function kl_empirical(y1, y2, unit) result(value)
      real(dp), intent(in) :: y1(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: y2(:) !! Nonnegative finite bin counts with positive total.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = kl_plugin(freqs_empirical(y1), freqs_empirical(y2), unit)
   end function kl_empirical

   pure function kl_dirichlet(y1, y2, a1, a2, unit) result(value)
      real(dp), intent(in) :: y1(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: y2(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: a1 !! Nonnegative finite pseudocount added to every cell of the corresponding sample.
      real(dp), intent(in) :: a2 !! Nonnegative finite pseudocount added to every cell of the corresponding sample.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = kl_plugin(freqs_dirichlet(y1, a1), freqs_dirichlet(y2, a2), unit)
   end function kl_dirichlet

   pure function kl_shrink(y1, y2, lambda1, lambda2, unit) result(value)
      real(dp), intent(in) :: y1(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: y2(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in), optional :: lambda1 !! Shrinkage intensity in [0,1] for the corresponding sample; absent means estimate.
      real(dp), intent(in), optional :: lambda2 !! Shrinkage intensity in [0,1] for the corresponding sample; absent means estimate.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = kl_plugin(freqs_shrink(y1, lambda1), freqs_shrink(y2, lambda2), unit)
   end function kl_shrink

   pure function chi2_empirical(y1, y2, unit) result(value)
      real(dp), intent(in) :: y1(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: y2(:) !! Nonnegative finite bin counts with positive total.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = chi2_plugin(freqs_empirical(y1), freqs_empirical(y2), unit)
   end function chi2_empirical

   pure function chi2_dirichlet(y1, y2, a1, a2, unit) result(value)
      real(dp), intent(in) :: y1(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: y2(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: a1 !! Nonnegative finite pseudocount added to every cell of the corresponding sample.
      real(dp), intent(in) :: a2 !! Nonnegative finite pseudocount added to every cell of the corresponding sample.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = chi2_plugin(freqs_dirichlet(y1, a1), freqs_dirichlet(y2, a2), unit)
   end function chi2_dirichlet

   pure function chi2_shrink(y1, y2, lambda1, lambda2, unit) result(value)
      real(dp), intent(in) :: y1(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: y2(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in), optional :: lambda1 !! Shrinkage intensity in [0,1] for the corresponding sample; absent means estimate.
      real(dp), intent(in), optional :: lambda2 !! Shrinkage intensity in [0,1] for the corresponding sample; absent means estimate.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp) :: value
      value = chi2_plugin(freqs_shrink(y1, lambda1), freqs_shrink(y2, lambda2), unit)
   end function chi2_shrink

   pure function freqs(y, method, lambda) result(prob)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts; returned probabilities preserve its length.
      character(len=*), intent(in), optional :: method !! ML (default), MM, Jeffreys, Laplace, SG, minimax, CS, NSB, or shrink.
      real(dp), intent(in), optional :: lambda !! Fixed shrinkage intensity in [0,1], used only by method shrink.
      real(dp) :: prob(size(y)), a, n
      character(len=:), allocatable :: selected
      selected = 'ML'
      if (present(method)) selected = method
      prob = missing()
      if (.not. valid_counts(y)) return
      select case (selected)
      case ('ML')
         prob = freqs_empirical(y)
      case ('Jeffreys')
         prob = freqs_dirichlet(y, 0.5_dp)
      case ('Laplace')
         prob = freqs_dirichlet(y, 1.0_dp)
      case ('SG')
         prob = freqs_dirichlet(y, 1.0_dp / real(size(y), dp))
      case ('minimax')
         n = count_total(y)
         if (.not. ieee_is_finite(n)) return
         a = sqrt(n) / real(size(y), dp)
         prob = freqs_dirichlet(y, a)
      case ('shrink')
         prob = freqs_shrink(y, lambda)
      case default
         ! MM, CS, and NSB have no frequency estimator in upstream R.
         prob = missing()
      end select
   end function freqs

   pure function entropy(y, method, unit, lambda) result(value)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts; invalid data or unsupported method yields NaN.
      character(len=*), intent(in), optional :: method !! ML (default), MM, Jeffreys, Laplace, SG, minimax, CS, or shrink.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      real(dp), intent(in), optional :: lambda !! Fixed shrinkage intensity in [0,1], used only by method shrink.
      real(dp) :: value
      character(len=:), allocatable :: selected
      selected = 'ML'
      if (present(method)) selected = method
      select case (selected)
      case ('MM')
         value = entropy_millermadow(y, unit)
      case ('CS')
         value = entropy_chaoshen(y, unit)
      case default
         value = entropy_plugin(freqs(y, selected, lambda), unit)
      end select
   end function entropy

   pure elemental function make_test(statistic, df, unit) result(test)
      real(dp), intent(in) :: statistic !! Test statistic in the requested logarithm unit, before upper-tail evaluation.
      integer, intent(in) :: df !! Chi-squared degrees of freedom; zero is treated as a point mass at zero.
      character(len=*), intent(in), optional :: unit !! Logarithm unit used by statistic: log (default), log2, or log10.
      type(test_result) :: test
      real(dp) :: natural_stat
      test%stat = statistic
      test%df = df
      test%pval = missing()
      natural_stat = statistic * unit_divisor(unit)
      if (ieee_is_nan(natural_stat) .or. df < 0) return
      if (natural_stat < 0.0_dp) then
         natural_stat = 0.0_dp
      end if
      if (.not. ieee_is_finite(natural_stat)) then
         test%pval = 0.0_dp
      else if (df == 0) then
         test%pval = 0.0_dp
         if (natural_stat == 0.0_dp) test%pval = 1.0_dp
      else
         test%pval = r_pchisq(natural_stat, real(df, dp), lower_tail=.false.)
      end if
   end function make_test

   pure function gstat(y, expected, unit) result(test)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: expected(:) !! Same-length nonnegative reference weights; normalized internally.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      type(test_result) :: test
      real(dp) :: n, statistic
      integer :: df
      n = count_total(y)
      df = size(y) - 1
      statistic = 2.0_dp * (n * kl_empirical(y, expected, unit))
      test = make_test(statistic, df, unit)
   end function gstat

   pure function chi2stat(y, expected, unit) result(test)
      real(dp), intent(in) :: y(:) !! Nonnegative finite bin counts with positive total.
      real(dp), intent(in) :: expected(:) !! Same-length nonnegative reference weights; normalized internally.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      type(test_result) :: test
      real(dp) :: n, statistic
      integer :: df
      n = count_total(y)
      df = size(y) - 1
      statistic = n * chi2_empirical(y, expected, unit)
      test = make_test(statistic, df, unit)
   end function chi2stat

   pure function gstatindep(y, unit) result(test)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      type(test_result) :: test
      real(dp) :: n, statistic
      integer :: df
      n = count_total(reshape(y, [size(y)]))
      df = (size(y, 1) - 1) * (size(y, 2) - 1)
      statistic = 2.0_dp * (n * mi_empirical(y, unit))
      test = make_test(statistic, df, unit)
   end function gstatindep

   pure function chi2statindep(y, unit) result(test)
      real(dp), intent(in) :: y(:, :) !! Nonnegative finite contingency-table counts with positive total.
      character(len=*), intent(in), optional :: unit !! Exact logarithm unit: log (default), log2, or log10.
      type(test_result) :: test
      real(dp) :: n, statistic
      integer :: df
      n = count_total(reshape(y, [size(y)]))
      df = (size(y, 1) - 1) * (size(y, 2) - 1)
      statistic = n * chi2indep_empirical(y, unit)
      test = make_test(statistic, df, unit)
   end function chi2statindep

   pure elemental function equal_width_break(index, lower, upper, num_bins) result(edge)
      integer, intent(in) :: index !! Interior boundary index between one and num_bins minus one.
      real(dp), intent(in) :: lower !! Finite lower endpoint of the strictly increasing bin range.
      real(dp), intent(in) :: upper !! Finite upper endpoint of the strictly increasing bin range.
      integer, intent(in) :: num_bins !! Positive number of intervals in the range.
      real(dp) :: edge, fraction, width
      fraction = real(index, dp) / real(num_bins, dp)
      if (lower < 0.0_dp .and. upper > 0.0_dp) then
         if (upper > huge(1.0_dp) + lower) then
            edge = (1.0_dp - fraction) * lower + fraction * upper
            return
         end if
      end if
      width = upper - lower
      edge = lower + real(index, dp) * (width / real(num_bins, dp))
   end function equal_width_break

   pure function distinct_breaks(lower, upper, num_bins) result(valid)
      real(dp), intent(in) :: lower !! Finite lower endpoint, strictly smaller than upper.
      real(dp), intent(in) :: upper !! Finite upper endpoint of the requested histogram range.
      integer, intent(in) :: num_bins !! Positive bin count; all computed breaks must remain distinct in dp precision.
      logical :: valid
      real(dp) :: previous, edge
      integer :: i
      valid = .false.
      previous = lower
      do i = 1, num_bins - 1
         edge = equal_width_break(i, lower, upper, num_bins)
         if (edge <= previous .or. edge >= upper) return
         previous = edge
      end do
      valid = .true.
   end function distinct_breaks

   pure elemental function bin_index(x, lower, upper, num_bins) result(index)
      real(dp), intent(in) :: x !! Observation; nonfinite values and values outside the closed range receive index zero.
      real(dp), intent(in) :: lower !! Finite lower endpoint, included in the first interval.
      real(dp), intent(in) :: upper !! Finite upper endpoint greater than lower, included in the last interval.
      integer, intent(in) :: num_bins !! Positive number of equal-width, right-closed intervals.
      integer :: index, lo, hi, mid
      real(dp) :: edge
      index = 0
      if (.not. ieee_is_finite(x)) return
      if (x < lower .or. x > upper) return
      lo = 1
      hi = num_bins
      do while (lo < hi)
         mid = lo + (hi - lo) / 2
         edge = equal_width_break(mid, lower, upper, num_bins)
         if (x <= edge) then
            hi = mid
         else
            lo = mid + 1
         end if
      end do
      index = lo
   end function bin_index

   pure subroutine get_range(x, bounds, valid, requested)
      real(dp), intent(in) :: x(:) !! Observations whose finite minimum and maximum define the default range.
      real(dp), intent(out) :: bounds(2) !! Validated finite lower and upper endpoints, or zeros on failure.
      logical, intent(out) :: valid !! True only when a finite strictly increasing range can be established.
      real(dp), intent(in), optional :: requested(2) !! Explicit finite increasing endpoints; bypass default range computation.
      bounds = 0.0_dp
      valid = .false.
      if (present(requested)) then
         bounds = requested
      else
         if (size(x) == 0) return
         if (.not. all(ieee_is_finite(x))) return
         bounds = [minval(x), maxval(x)]
      end if
      if (.not. all(ieee_is_finite(bounds))) return
      valid = bounds(2) > bounds(1)
   end subroutine get_range

   pure subroutine discretize(x, num_bins, counts, range, status)
      real(dp), intent(in) :: x(:) !! Observations; nonfinite or out-of-range entries are omitted when range is explicit.
      integer, intent(in) :: num_bins !! Positive number of equal-width, right-closed bins; lowest endpoint is included.
      integer, allocatable, intent(out) :: counts(:) !! Allocated bin counts; zeros on error, or zero length for invalid bin count.
      real(dp), intent(in), optional :: range(2) !! Finite increasing endpoints; default is [minval(x),maxval(x)].
      integer, intent(out), optional :: status !! Zero on success, one for invalid bin count, data-derived range, or endpoints.
      real(dp) :: bounds(2)
      integer :: i, index
      logical :: valid
      if (present(status)) status = 1
      allocate (counts(max(num_bins, 0)))
      counts = 0
      if (num_bins <= 0) return
      call get_range(x, bounds, valid, range)
      if (.not. valid) return
      if (.not. distinct_breaks(bounds(1), bounds(2), num_bins)) return
      do i = 1, size(x)
         index = bin_index(x(i), bounds(1), bounds(2), num_bins)
         if (index > 0) counts(index) = counts(index) + 1
      end do
      if (present(status)) status = 0
   end subroutine discretize

   pure subroutine discretize2d(x1, x2, num_bins1, num_bins2, counts, range1, range2, status)
      real(dp), intent(in) :: x1(:) !! First coordinate of each paired observation, corresponding to output rows.
      real(dp), intent(in) :: x2(:) !! Second coordinate, corresponding to columns; must have the same length as x1.
      integer, intent(in) :: num_bins1 !! Positive number of equal-width intervals for the first coordinate.
      integer, intent(in) :: num_bins2 !! Positive number of equal-width intervals for the second coordinate.
      integer, allocatable, intent(out) :: counts(:, :) !! Allocated joint counts; pairs missing either interval are omitted.
      real(dp), intent(in), optional :: range1(2) !! First finite increasing range; default is the full range of x1.
      real(dp), intent(in), optional :: range2(2) !! Second finite increasing range; default is the full range of x2.
      integer, intent(out), optional :: status !! Zero on success, one for invalid bins, shape mismatch, or invalid range.
      real(dp) :: bounds1(2), bounds2(2)
      integer :: i, row, col
      logical :: valid1, valid2
      if (present(status)) status = 1
      allocate (counts(max(num_bins1, 0), max(num_bins2, 0)))
      counts = 0
      if (num_bins1 <= 0 .or. num_bins2 <= 0) return
      if (size(x1) /= size(x2)) return
      call get_range(x1, bounds1, valid1, range1)
      call get_range(x2, bounds2, valid2, range2)
      if (.not. valid1 .or. .not. valid2) return
      if (.not. distinct_breaks(bounds1(1), bounds1(2), num_bins1)) return
      if (.not. distinct_breaks(bounds2(1), bounds2(2), num_bins2)) return
      do i = 1, size(x1)
         row = bin_index(x1(i), bounds1(1), bounds1(2), num_bins1)
         col = bin_index(x2(i), bounds2(1), bounds2(2), num_bins2)
         if (row > 0 .and. col > 0) counts(row, col) = counts(row, col) + 1
      end do
      if (present(status)) status = 0
   end subroutine discretize2d

end module entropy_module
