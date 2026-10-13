! SPDX-License-Identifier: GPL-3.0-only
module hdinterval_test_callbacks
   use hdinterval, only: dp
   use r_distributions, only: r_qnorm
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
   implicit none
contains
   pure function normal_quantile(p) result(q)
      real(dp), intent(in) :: p !! Standard normal lower-tail probability in (0,1).
      real(dp) :: q
      q = r_qnorm(p)
   end function normal_quantile

   pure function normal_cdf(x) result(p)
      real(dp), intent(in) :: x !! Standard normal coordinate.
      real(dp) :: p
      p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
   end function normal_cdf

   pure function logistic_quantile(p) result(q)
      real(dp), intent(in) :: p !! Logistic lower-tail probability in (0,1).
      real(dp) :: q
      q = log(p / (1.0_dp - p))
   end function logistic_quantile

   pure function logistic_cdf(x) result(p)
      real(dp), intent(in) :: x !! Standard logistic coordinate; formula avoids overflow.
      real(dp) :: p, e
      e = exp(-abs(x))
      if (x < 0.0_dp) then
         p = e / (1.0_dp + e)
      else
         p = 1.0_dp / (1.0_dp + e)
      end if
   end function logistic_cdf

   pure function shifted_cdf(x) result(p)
      real(dp), intent(in) :: x !! Coordinate for logistic location 1000 and scale 7.
      real(dp) :: p
      p = logistic_cdf((x - 1000.0_dp) / 7.0_dp)
   end function shifted_cdf

   pure function exponential_quantile(p) result(q)
      real(dp), intent(in) :: p !! Unit exponential lower-tail probability in (0,1).
      real(dp) :: q
      q = -log(1.0_dp - p)
   end function exponential_quantile

   pure function beta21_quantile(p) result(q)
      real(dp), intent(in) :: p !! Beta(2,1) lower-tail probability in (0,1).
      real(dp) :: q
      q = sqrt(p)
   end function beta21_quantile

   pure function cauchy_quantile(p) result(q)
      real(dp), intent(in) :: p !! Standard Cauchy lower-tail probability in (0,1).
      real(dp) :: q
      q = tan(acos(-1.0_dp) * (p - 0.5_dp))
   end function cauchy_quantile

   pure function cauchy_cdf(x) result(p)
      real(dp), intent(in) :: x !! Standard Cauchy coordinate.
      real(dp) :: p
      p = 0.5_dp + atan(x) / acos(-1.0_dp)
   end function cauchy_cdf

   pure function uniform_quantile(p) result(q)
      real(dp), intent(in) :: p !! Uniform(0,1) lower-tail probability.
      real(dp) :: q
      q = p
   end function uniform_quantile

   pure function extreme_quantile(p) result(q)
      real(dp), intent(in) :: p !! Probability selecting extreme finite endpoints whose width overflows.
      real(dp) :: q
      q = merge(-huge(1.0_dp), huge(1.0_dp), p < 0.5_dp)
   end function extreme_quantile

   pure function decreasing_quantile(p) result(q)
      real(dp), intent(in) :: p !! Probability passed to an intentionally invalid decreasing quantile.
      real(dp) :: q
      q = 1.0_dp - p
   end function decreasing_quantile

   pure function invalid_callback(x) result(y)
      real(dp), intent(in) :: x !! Coordinate ignored by a deliberately invalid callback except to form its NaN kind.
      real(dp) :: y
      y = ieee_value(x, ieee_quiet_nan)
   end function invalid_callback

   pure function bounded_cdf(x) result(p)
      real(dp), intent(in) :: x !! Uniform(0,1) CDF, deliberately unavailable outside its support.
      real(dp) :: p
      if (x < 0.0_dp .or. x > 1.0_dp) then
         p = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         p = x
      end if
   end function bounded_cdf

   pure function decreasing_cdf(x) result(p)
      real(dp), intent(in) :: x !! Coordinate passed to a deliberately decreasing probability function.
      real(dp) :: p
      p = 1.0_dp - logistic_cdf(x)
   end function decreasing_cdf

   pure function nonmonotone_cdf(x) result(p)
      real(dp), intent(in) :: x !! Coordinate for a probability curve with an invalid internal reversal.
      real(dp) :: p
      if (x < 0.0_dp) then
         p = 0.4_dp
      else if (x == 0.0_dp) then
         p = 0.1_dp
      else
         p = 0.6_dp
      end if
   end function nonmonotone_cdf

   pure function outside_range_cdf(x) result(p)
      real(dp), intent(in) :: x !! Coordinate for an invalid CDF taking values outside [0,1].
      real(dp) :: p
      p = 2.0_dp + logistic_cdf(x)
   end function outside_range_cdf

   pure function constant_cdf(x) result(p)
      real(dp), intent(in) :: x !! Finite coordinate for a deliberately non-bracketing constant probability curve.
      real(dp) :: p
      p = 0.25_dp + 0.0_dp * x
   end function constant_cdf
end module hdinterval_test_callbacks

program test_hdinterval
   use hdinterval
   use r_reference_values
   use hdinterval_test_callbacks
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   real(dp) :: bounds(2), columns(2, 2), matrix(5, 2), mass, height, nan, infinity, tail, q
   real(dp), allocatable :: intervals(:, :), old_intervals(:, :)
   real(dp) :: grid(7), density(7), quantiles(9), qthree(3), qempty(0), bempty(2, 0)
   integer :: status, statuses(2), codes(3), failures, checks, i
   logical :: is_hdi
   procedure(scalar_function_interface), pointer :: callback

   failures = 0
   checks = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   infinity = ieee_value(0.0_dp, ieee_positive_inf)
   call hdi_samples([9.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 8.0_dp, 20.0_dp], bounds, 0.6_dp, status)
   call check(status == hdi_ok, 'sample status')
   call close(bounds, sample_regular, 0.0_dp, 'R sample unsorted')
   call hdi_samples([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp], bounds, 0.6_dp)
   call close(bounds, sample_integer_span, 0.0_dp, 'R floor mass and first tie')
   call hdi_samples([3.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 9.0_dp], bounds, 0.5_dp)
   call close(bounds, sample_ties, 0.0_dp, 'R tied observations')
   call hdi_samples([4.0_dp, 2.0_dp, 7.0_dp], bounds, 0.1_dp)
   call close(bounds, sample_tiny_mass, 0.0_dp, 'R less than one draw mass')
   call hdi_samples([nan, 5.0_dp, nan, 1.0_dp, 2.0_dp, 8.0_dp, 10.0_dp], bounds, 0.6_dp)
   call close(bounds, sample_missing, 0.0_dp, 'R removes missing draws')
   matrix(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp]
   matrix(:, 2) = [2.0_dp, 4.0_dp, 8.0_dp, 16.0_dp, 32.0_dp]
   call hdi_columns(matrix, columns, 0.6_dp, statuses)
   call check(all(statuses == hdi_ok), 'matrix statuses')
   call close(reshape(columns, [4]), sample_matrix, 0.0_dp, 'R column matrix')
   matrix(:, 2) = nan
   call hdi_columns(matrix, columns, statuses = statuses)
   call check(statuses(1) == hdi_ok .and. statuses(2) == hdi_no_interval, 'mixed matrix statuses')
   call hdi_columns(reshape([real(dp) ::], [0, 0]), bempty)
   call check(size(bempty) == 0, 'empty matrix shape')
   call hdi_samples([real(dp) ::], bounds, status = status)
   call check(status == hdi_no_interval .and. all(ieee_is_nan(bounds)), 'empty samples')
   call hdi_samples([nan, nan], bounds, status = status)
   call check(status == hdi_no_interval .and. all(ieee_is_nan(bounds)), 'all missing samples')
   call hdi_samples([42.0_dp], bounds, status = status)
   call close(bounds, [42.0_dp, 42.0_dp], 0.0_dp, 'singleton')
   call hdi_samples([infinity, infinity], bounds, status = status)
   call check(status == hdi_ok .and. all(bounds == infinity), 'identical infinities fallback')
   call hdi_samples([-infinity, infinity], bounds, status = status)
   call check(status == hdi_ok .and. bounds(1) == -infinity .and. bounds(2) == infinity, 'opposite infinities')
   call hdi_samples([-infinity, 1.0_dp, 2.0_dp, 3.0_dp, infinity], bounds, 0.3_dp, status)
   call close(bounds, [1.0_dp, 2.0_dp], 0.0_dp, 'finite narrow interval among infinities')
   call hdi_samples([-huge(1.0_dp), huge(1.0_dp)], bounds, status = status)
   call check(status == hdi_ok .and. bounds(1) == -huge(1.0_dp) .and. &
      bounds(2) == huge(1.0_dp), 'overflow-sized sample interval')
   call hdi_samples([1.0_dp, 2.0_dp], bounds, 0.0_dp, status)
   call check(status == hdi_invalid_mass .and. all(ieee_is_nan(bounds)), 'zero mass')
   call hdi_samples([1.0_dp, 2.0_dp], bounds, 1.0_dp, status)
   call check(status == hdi_invalid_mass, 'unit mass')
   call hdi_samples([1.0_dp, 2.0_dp], bounds, nan, status)
   call check(status == hdi_invalid_mass, 'NaN mass')

   grid = [(real(i, dp), i = -3, 3)]
   density = [3.0_dp, 4.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 4.0_dp, 1.0_dp]
   call hdi_density(grid, density, intervals, height, 0.6_dp, .true., status, is_hdi, mass)
   call check(status == hdi_ok .and. is_hdi .and. size(intervals, 2) == 2, 'split density status and shape')
   call compare_r_density(intervals, density_split, 'R split density')
   call close([height], density_height, 0.0_dp, 'R density threshold')
   call close([mass], [11.0_dp / 14.0_dp], 1e-15_dp, 'threshold mass includes density ties')
   call hdi_density(grid, density, intervals, height, 0.6_dp, status = status, is_hdi = is_hdi)
   call check(status == hdi_ok .and. .not. is_hdi, 'nonsplit fallback flag')
   call close(intervals(:, 1), density_joined, 1e-14_dp, 'R nonsplit interval')
   call hdi_density(mixture_grid_x, mixture_grid_y, intervals, height, 0.8_dp, .true., status, is_hdi, mass)
   call check(status == hdi_ok .and. is_hdi .and. mass >= 0.8_dp, 'mixture threshold mass')
   call compare_r_density(intervals, mixture_split, 'R multimodal Gaussian grid')
   call close([height], mixture_height, 2e-15_dp, 'R mixture cutoff height')
   old_intervals = intervals
   call hdi_density(mixture_grid_x, 1e300_dp * mixture_grid_y, intervals, height, 0.8_dp, .true., status)
   call check(status == hdi_ok, 'large density normalization')
   call close(reshape(intervals, [size(intervals)]), reshape(old_intervals, [size(old_intervals)]), &
      0.0_dp, 'density scale invariant')
   call hdi_density(mixture_grid_x, 1e-300_dp * mixture_grid_y, intervals, height, 0.8_dp, .true., status)
   call compare_r_density(intervals, mixture_split, 'subnormal-scale density')
   call hdi_density(mixture_grid_x, mixture_grid_y, intervals, height, 0.8_dp, status = status, is_hdi = is_hdi)
   call close(intervals(:, 1), mixture_joined, 2e-14_dp, 'R mixture shortest contiguous fallback')
   call check(.not. is_hdi, 'mixture fallback indicated')
   call hdi_density(grid, spread(1.0_dp, 1, 7), intervals, height, 0.2_dp, status = status, achieved_mass = mass)
   call close(intervals(:, 1), [-3.0_dp, 3.0_dp], 0.0_dp, 'constant-density ties retain all nodes')
   call close([mass], [1.0_dp], 0.0_dp, 'constant-density attained mass')
   call hdi_density([0.0_dp, 1.0_dp, 2.0_dp], [1.0_dp, 0.0_dp, 1.0_dp], intervals, height, &
      0.9_dp, status = status)
   call check(status == hdi_no_interval .and. size(intervals, 2) == 0, 'undefined R contiguous fallback')
   call hdi_density([0.0_dp, 1.0_dp, 2.1_dp], [1.0_dp, 2.0_dp, 1.0_dp], intervals, height, status = status)
   call check(status == hdi_invalid_data, 'nonuniform grid rejected')
   call hdi_density([-huge(1.0_dp), huge(1.0_dp)], [1.0_dp, 1.0_dp], intervals, height, status = status)
   call check(status == hdi_invalid_data, 'unrepresentable density spacing rejected')
   call hdi_density([0.0_dp, 0.0_dp], [1.0_dp, 1.0_dp], intervals, height, status = status)
   call check(status == hdi_invalid_data, 'duplicate grid rejected')
   call hdi_density([0.0_dp, nan], [1.0_dp, 1.0_dp], intervals, height, status = status)
   call check(status == hdi_invalid_data, 'NaN grid rejected')
   call hdi_density(grid, -density, intervals, height, status = status)
   call check(status == hdi_invalid_data, 'negative density rejected')
   call hdi_density(grid, density * 0.0_dp, intervals, height, status = status)
   call check(status == hdi_invalid_data, 'zero density rejected')
   call hdi_density(grid, [1.0_dp, nan], intervals, height, status = status)
   call check(status == hdi_invalid_data, 'mismatched density rejected')
   call hdi_density(grid, density, intervals, height, nan, status = status)
   call check(status == hdi_invalid_mass, 'NaN density mass')

   callback => normal_quantile
   call hdi_function(callback, bounds, status = status, lower_tail = tail)
   call check(status == hdi_ok, 'normal function status')
   call close(bounds, function_normal, 3e-7_dp, 'R normal quantile HDI')
   call close([tail], [0.025_dp], 2e-8_dp, 'normal symmetric tail')
   call close(bounds, [-1.959963984540054_dp, 1.959963984540054_dp], 3e-7_dp, 'known normal HDI')
   call hdi_function(callback, bounds, 0.3_dp)
   call close(bounds, function_normal_small, 2e-7_dp, 'R HDI mass under one half')
   callback => logistic_quantile
   call hdi_function(callback, bounds, 0.8_dp)
   call close(bounds, function_logistic, 2e-7_dp, 'R logistic quantile HDI')
   call close(bounds, [-log(9.0_dp), log(9.0_dp)], 2e-7_dp, 'analytic logistic HDI')
   callback => exponential_quantile
   call hdi_function(callback, bounds, 0.95_dp)
   call close(bounds, function_exponential, 3e-7_dp, 'R boundary exponential HDI')
   call close(bounds, [0.0_dp, -log(0.05_dp)], 3e-7_dp, 'analytic exponential HDI')
   callback => beta21_quantile
   call hdi_function(callback, bounds, 0.8_dp)
   call close(bounds, function_beta21, 1e-7_dp, 'R bounded Beta HDI')
   call close(bounds, [sqrt(0.2_dp), 1.0_dp], 1e-7_dp, 'analytic Beta HDI')
   callback => cauchy_quantile
   call hdi_function(callback, bounds, 0.7_dp)
   call close(bounds, function_cauchy, 3e-7_dp, 'R Cauchy quantile HDI')
   callback => uniform_quantile
   call hdi_function(callback, bounds, 0.4_dp, status = status)
   call check(status == hdi_ok .and. bounds(1) >= 0.0_dp .and. bounds(2) <= 1.0_dp, 'flat-width callback bounds')
   call close([bounds(2) - bounds(1)], [0.4_dp], 1e-14_dp, 'nonunique uniform interval width')
   callback => extreme_quantile
   call hdi_function(callback, bounds, status = status)
   call check(status == hdi_invalid_callback, 'overflow-sized quantile width rejected')
   callback => decreasing_quantile
   call hdi_function(callback, bounds, status = status)
   call check(status == hdi_invalid_callback .and. all(ieee_is_nan(bounds)), 'decreasing quantile rejected')
   callback => invalid_callback
   call hdi_function(callback, bounds, status = status)
   call check(status == hdi_invalid_callback, 'NaN quantile rejected')
   callback => logistic_quantile
   call hdi_function(callback, bounds, status = status, max_iterations = 1)
   call check(status == hdi_no_convergence, 'golden iteration budget')
   call hdi_function(callback, bounds, tol = nan, status = status)
   call check(status == hdi_invalid_options, 'NaN optimizer tolerance')
   nullify(callback)
   call hdi_function(callback, bounds, status = status)
   call check(status == hdi_unassociated_callback, 'null quantile pointer')

   callback => logistic_cdf
   call inverse_cdf(probabilities, callback, quantiles, status = status)
   call check(status == hdi_ok, 'inverse vector status')
   call close(quantiles, inverse_logistic, 2e-4_dp, 'R logistic inverse default tolerance')
   call inverse_cdf(probabilities, callback, quantiles, tol = 1e-10_dp, status = status)
   do i = 1, size(probabilities)
      q = max(1e-6_dp, min(1.0_dp - 1e-6_dp, probabilities(i)))
      call close([quantiles(i)], [log(q / (1.0_dp - q))], 2e-9_dp, 'analytic clipped logistic inverse')
   end do
   call check(quantiles(1) == quantiles(2) .and. quantiles(8) == quantiles(9), 'endpoint clipping equivalence')
   callback => normal_cdf
   call inverse_cdf(probabilities, callback, quantiles, status = status)
   call close(quantiles, inverse_normal, 2e-4_dp, 'R normal inverse')
   call inverse_cdf(0.5_dp, callback, q, tol = 1e-12_dp, status = status)
   call close([q], [0.0_dp], 0.0_dp, 'normal median exact root')
   callback => shifted_cdf
   call inverse_cdf([0.1_dp, 0.5_dp, 0.9_dp], callback, qthree, status = status)
   call close(qthree, inverse_shifted, 2e-4_dp, 'R bracket extension to location 1000')
   call inverse_cdf(0.5_dp, callback, q, status = status, max_extensions = 0)
   call check(status == hdi_no_convergence .and. ieee_is_nan(q), 'disabled extension reports failure')
   callback => cauchy_cdf
   call inverse_cdf([0.01_dp, 0.5_dp, 0.99_dp], callback, qthree, status = status)
   call close(qthree, inverse_cauchy, 2e-4_dp, 'R Cauchy tail inversion')
   callback => bounded_cdf
   call inverse_cdf(0.25_dp, callback, q, status = status)
   call check(status == hdi_invalid_callback, 'restricted-domain default bracket rejected')
   call inverse_cdf(0.25_dp, callback, q, lower = 0.0_dp, upper = 1.0_dp, status = status)
   call check(status == hdi_ok, 'restricted domain explicit bracket accepted')
   call close([q], [0.25_dp], 0.0_dp, 'restricted domain exact quantile')
   callback => decreasing_cdf
   call inverse_cdf(0.3_dp, callback, q, status = status)
   call check(status == hdi_invalid_callback, 'decreasing CDF rejected')
   callback => nonmonotone_cdf
   call inverse_cdf(0.5_dp, callback, q, status = status)
   call check(status == hdi_invalid_callback, 'internal CDF monotonicity failure')
   callback => outside_range_cdf
   call inverse_cdf(0.3_dp, callback, q, status = status)
   call check(status == hdi_invalid_callback, 'out of range CDF rejected')
   callback => invalid_callback
   call inverse_cdf(0.3_dp, callback, q, status = status)
   call check(status == hdi_invalid_callback, 'NaN CDF rejected')
   callback => constant_cdf
   call inverse_cdf(0.7_dp, callback, q, status = status, max_extensions = 5)
   call check(status == hdi_no_convergence, 'constant CDF cannot bracket')
   callback => logistic_cdf
   call inverse_cdf([0.25_dp, nan, 2.0_dp], callback, qthree, status = status, element_status = codes)
   call check(status == hdi_invalid_data .and. codes(1) == hdi_ok .and. &
      all(codes(2:) == hdi_invalid_data), 'partial vector failure status')
   call close(qthree(:1), [-log(3.0_dp)], 1e-4_dp, 'partial vector success retained')
   call check(all(ieee_is_nan(qthree(2:))), 'failed vector elements are NaNs')
   call inverse_cdf([real(dp) ::], callback, qempty, status = status)
   call check(status == hdi_ok .and. size(qempty) == 0, 'empty inverse vector')
   call inverse_cdf(0.3_dp, callback, q, tol = nan, status = status)
   call check(status == hdi_invalid_options, 'NaN inverse tolerance')
   call inverse_cdf(0.3_dp, callback, q, lower = 10.0_dp, upper = -10.0_dp, status = status)
   call check(status == hdi_invalid_options, 'reversed bracket')
   call inverse_cdf(0.3_dp, callback, q, max_iterations = 1, status = status)
   call check(status == hdi_no_convergence, 'inverse refinement budget')
   nullify(callback)
   call inverse_cdf(0.3_dp, callback, q, status = status)
   call check(status == hdi_unassociated_callback, 'null CDF pointer')

   if (failures /= 0) then
      print '(a,i0,a,i0)', 'FAIL: ', failures, ' of ', checks
      error stop 1
   end if
   print '(a,i0,a)', 'PASS: ', checks, ' deterministic HDInterval checks'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion that must hold for this deterministic check.
      character(len=*), intent(in) :: label !! Human-readable case name printed if the assertion fails.
      checks = checks + 1
      if (.not. condition) then
         failures = failures + 1
         print '(a)', 'FAIL: ' // label
      end if
   end subroutine check

   subroutine close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Computed finite vector for comparison.
      real(dp), intent(in) :: expected(:) !! Reference vector with the same shape as actual.
      real(dp), intent(in) :: tolerance !! Absolute elementwise error permitted for this check.
      character(len=*), intent(in) :: label !! Name of the numerical reference or identity being checked.
      call check(size(actual) == size(expected), label // ' shape')
      if (size(actual) /= size(expected)) return
      call check(all(abs(actual - expected) <= tolerance), label)
      if (.not. all(abs(actual - expected) <= tolerance)) then
         print *, 'actual: ', actual
         print *, 'expected: ', expected
      end if
   end subroutine close

   subroutine compare_r_density(actual, expected, label)
      real(dp), intent(in) :: actual(:, :) !! Fortran interval endpoints, shape (2,k).
      real(dp), intent(in) :: expected(:) !! Flattened R k-by-2 matrix, with begins followed by ends.
      character(len=*), intent(in) :: label !! Name of the R density fixture being compared.
      integer :: k
      k = size(actual, 2)
      call check(2 * k == size(expected), label // ' shape')
      if (2 * k /= size(expected)) return
      call close(actual(1, :), expected(:k), 1e-14_dp, label // ' lower')
      call close(actual(2, :), expected(k + 1:), 1e-14_dp, label // ' upper')
   end subroutine compare_r_density
end program test_hdinterval
