! SPDX-License-Identifier: GPL-2.0-or-later
program test_properties
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_nan, ieee_is_finite
   use diptest_api
   implicit none
   type(dip_result) :: fit, reference
   type(dip_test_result) :: test, test_again
   real(dp) :: x(127), pair_masses(40), extremes(20), small(8), base(8), uniforms(4, 4), p, nan, infinity, scale_tiny
   real(dp), allocatable :: sample(:)
   integer, allocatable :: seed(:)
   integer :: i, k, n, status, seed_size

   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   infinity = ieee_value(0.0_dp, ieee_positive_inf)
   call dip_compute([real(dp) ::], fit)
   call require(fit%status == dip_empty_sample .and. ieee_is_nan(fit%statistic), 'empty sample status')
   call require(ieee_is_nan(dip([nan])), 'dip rejects NaN')
   call dip_compute([1.0_dp, infinity], fit)
   call require(fit%status == dip_nonfinite, 'dip rejects infinity')
   call dip_compute([2.0_dp, 1.0_dp], fit, sorted=.true.)
   call require(fit%status == dip_unsorted, 'checked presorted contract')
   call dip_compute([1.0_dp, 2.0_dp], fit, sorted=.true.)
   call require(fit%status == dip_ok, 'valid presorted contract')
   call dip_test([1.0_dp, nan, 1.0_dp, 2.0_dp, 2.0_dp], test)
   call require(test%n == 4 .and. test%status == dip_ok, 'NaN omission from dip.test')
   call near(test%statistic, 0.25_dp, 0.0_dp, 'NaN omitted statistic')
   call dip_test([nan, nan], test)
   call require(test%status == dip_empty_sample, 'all missing observations')
   call dip_test([1.0_dp, infinity, 2.0_dp], test)
   call require(test%status == dip_nonfinite, 'infinite test sample')
   do n = 1, 3
      call dip_test([(real(i, dp), i = 1, n)], test)
      call near(test%p_value, 1.0_dp, 0.0_dp, 'small sample p=1')
   end do
   call dip_compute([3.0_dp, 3.0_dp, 3.0_dp, 3.0_dp], fit)
   call near(fit%statistic, 0.125_dp, 0.0_dp, 'constant sample classical floor')
   call require(all(fit%gcm == [4, 1]) .and. all(fit%lcm == [1, 4]), 'defined constant hulls')
   call near(dip([3.0_dp, 3.0_dp], .true.), 0.0_dp, 0.0_dp, 'zero minimum constant sample')
   call dip_compute([3.0_dp], fit)
   call require(all(fit%full_gcm == [1]) .and. all(fit%full_lcm == [1]), 'singleton hulls')

   ! Independent closed-form geometry for every nonconstant two-point sample.
   do n = 2, 40
      do k = 1, n - 1
         pair_masses(:k) = -2.0_dp
         pair_masses(k + 1:n) = 5.0_dp
         call near(dip(pair_masses(:n)), real(min(k, n - k), dp) / (2.0_dp * real(n, dp)), &
            3.0e-15_dp, 'two-mass closed form')
      end do
      call near(dip([(real(i, dp), i = 1, n)]), 0.5_dp / real(n, dp), 3.0e-15_dp, 'uniform-grid closed form')
      call near(dip([(real(i, dp), i = 1, n)], .true.), 0.0_dp, 3.0e-15_dp, 'uniform-grid zero convention')
   end do

   do i = 1, size(x)
      x(i) = sin(real(i, dp) * 0.217_dp) + 0.019_dp * real(i, dp)
   end do
   call dip_compute(x, reference)
   call near(dip(-x), reference%statistic, 4.0e-15_dp, 'reflection invariance')
   call near(dip(x(size(x):1:-1)), reference%statistic, 4.0e-15_dp, 'permutation invariance')
   call near(dip(8.0_dp * x + 64.0_dp), reference%statistic, 4.0e-15_dp, 'affine invariance')
   call near(dip(scale(x, 900)), reference%statistic, 4.0e-15_dp, 'huge-scale invariance')
   call near(dip(scale(x, -900)), reference%statistic, 4.0e-15_dp, 'tiny-scale invariance')
   call check_hull_geometry(reference)

   base = [0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 101.0_dp, 102.0_dp, 103.0_dp, 104.0_dp]
   scale_tiny = nearest(0.0_dp, 1.0_dp)
   small = base * scale_tiny
   call near(dip(small), dip(base), 5.0e-15_dp, 'subnormal geometry')
   call near(dip(2.0_dp ** 40 + base / 1024.0_dp), dip(base), 5.0e-15_dp, 'large-offset geometry')
   call near(dip([-huge(1.0_dp), -huge(1.0_dp), huge(1.0_dp), huge(1.0_dp)]), &
      0.25_dp, 0.0_dp, 'opposite-sign overflowing differences')
   extremes(:9) = [(real(i, dp) * 1.0e-303_dp, i = 0, 8)]
   extremes(10:19) = [(real(i, dp) * 1.0e-303_dp, i = 1000, 1009)]
   extremes(20) = 1.0e300_dp
   call near(dip(extremes), 0.2232_dp, 3.0e-15_dp, 'tiny clusters with huge outlier')
   call near(dip(-extremes), 0.2232_dp, 3.0e-15_dp, 'reflected mixed extreme geometry')

   call dip_p_value(0.125_dp, 4, p, status)
   call near(p, 0.782_dp, 3.0e-15_dp, 'R approx averages repeated minimum knots')
   call dip_p_value(-1.0_dp, 4, p, status)
   call require(status == dip_invalid_argument .and. ieee_is_nan(p), 'invalid dip for calibration')
   call dip_p_value(0.1_dp, 0, p, status)
   call require(status == dip_invalid_argument, 'invalid n for calibration')
   call dip_p_value(nan, 4, p, status)
   call require(status == dip_invalid_argument, 'NaN calibration input')

   uniforms(:, 1) = [0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp]
   uniforms(:, 2) = [0.0_dp, 0.2_dp, 0.7_dp, 1.0_dp]
   uniforms(:, 3) = [0.3_dp, 0.3_dp, 0.3_dp, 0.3_dp]
   uniforms(:, 4) = [0.1_dp, 0.1_dp, 0.9_dp, 0.9_dp]
   call dip_test([0.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], test, uniforms)
   call near(test%p_value, 0.5_dp, 0.0_dp, 'MC inclusive >= and no plus-one adjustment')
   call dip_test([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], test, uniforms(:3, :))
   call require(test%status == dip_invalid_simulation, 'MC wrong row count')
   call dip_test([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], test, uniforms(:, :0))
   call require(test%status == dip_invalid_simulation, 'MC empty replicates')
   uniforms(1, 1) = -0.1_dp
   call dip_test([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], test, uniforms)
   call require(test%status == dip_invalid_simulation, 'MC out-of-range uniform')
   uniforms(1, 1) = nan
   call dip_test([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], test, uniforms)
   call require(test%status == dip_invalid_simulation, 'MC NaN uniform')
   call dip_test([0.0_dp, 1.0_dp], test, uniforms)
   call near(test%p_value, 1.0_dp, 0.0_dp, 'small n ignores Monte Carlo options')

   call random_seed(size=seed_size)
   allocate(seed(seed_size))
   seed = 217
   call random_seed(put=seed)
   call dip_test_simulate(x, test, 37)
   call random_seed(put=seed)
   call dip_test_simulate(x, test_again, 37)
   call require(test%status == dip_ok .and. test%replicates == 37, 'intrinsic simulation metadata')
   call near(test%p_value, test_again%p_value, 0.0_dp, 'intrinsic RNG replay within one compiler')
   call dip_test_simulate(x, test, 0)
   call require(test%status == dip_invalid_simulation, 'invalid intrinsic replicate count')
   call dip_test_simulate([1.0_dp], test, 0)
   call near(test%p_value, 1.0_dp, 0.0_dp, 'small n ignores intrinsic replicate count')
   call dip_test_simulate([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], test)
   call require(test%replicates == 2000 .and. test%status == dip_ok, 'default replicate count')
   allocate(sample(72001))
   sample = [(real(i, dp), i = 1, size(sample))]
   call dip_test(sample(:72000), test)
   call require(test%status == dip_ok .and. .not. test%asymptotic, 'exact 72000 table boundary')
   call near(test%p_value, 1.0_dp, 0.0_dp, '72000 uniform-grid test')
   call dip_test(sample, test)
   call require(test%status == dip_ok .and. test%asymptotic, 'above 72000 asymptotic calibration')
   print '(a)', 'Passed closed forms, hull geometry, invariances, extreme scales, input contracts, and RNG replay.'

contains

   subroutine require(condition, message)
      logical, intent(in) :: condition !! Assertion that must hold for the current mathematical or API property.
      character(*), intent(in) :: message !! Name of the failing property.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // message
         error stop 1
      end if
   end subroutine require

   subroutine near(actual, expected, tolerance, message)
      real(dp), intent(in) :: actual !! Numerical result being checked.
      real(dp), intent(in) :: expected !! Analytic or independently computed reference value.
      real(dp), intent(in) :: tolerance !! Allowed absolute error; zero requests bitwise numerical equality.
      character(*), intent(in) :: message !! Name of the failing numerical property.
      if (.not. (abs(actual - expected) <= tolerance)) then
         print '(a,2es25.16)', 'FAIL: ' // message // ': ', actual, expected
         error stop 1
      end if
   end subroutine near

   subroutine check_hull_geometry(result)
      type(dip_result), intent(in) :: result !! Complete fit to strictly increasing sorted data, used to check envelope geometry.
      integer :: side, segment, a, b, q
      integer, allocatable :: nodes(:)
      real(dp) :: chord, slope, old_slope
      do side = 1, 2
         if (side == 1) then
            nodes = result%full_gcm(size(result%full_gcm):1:-1)
         else
            nodes = result%full_lcm
         end if
         old_slope = 0.0_dp
         do segment = 1, size(nodes) - 1
            a = nodes(segment)
            b = nodes(segment + 1)
            slope = real(b - a, dp) / (result%x(b) - result%x(a))
            if (segment > 1) then
               if (side == 1) call require(slope >= old_slope - 1.0e-11_dp, 'convex full hull slopes')
               if (side == 2) call require(slope <= old_slope + 1.0e-11_dp, 'concave full hull slopes')
            end if
            do q = a, b
               chord = real(a, dp) + (result%x(q) - result%x(a)) * slope
               if (side == 1) call require(chord <= real(q, dp) + 1.0e-11_dp, 'minorant below ECDF index points')
               if (side == 2) call require(chord >= real(q, dp) - 1.0e-11_dp, 'majorant above ECDF index points')
            end do
            old_slope = slope
         end do
      end do
   end subroutine check_hull_geometry
end program test_properties
