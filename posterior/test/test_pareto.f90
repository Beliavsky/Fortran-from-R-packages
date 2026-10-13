! SPDX-License-Identifier: BSD-3-Clause
program test_pareto
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_is_nan, ieee_value, &
      ieee_quiet_nan, ieee_positive_inf, ieee_negative_inf
   use r_kinds, only: dp
   use posterior_pareto
   use posterior_predictive
   implicit none
   integer :: unit, cases, c, n, j, s, status, checks
   character(80) :: label
   character(5), parameter :: sides(3) = ['both ', 'left ', 'right']
   real(dp), allocatable :: x(:), weights(:), expected_x(:), smoothed(:), pointwise(:), expected_test(:)
   real(dp) :: expected(14), actual(14), targets(3), results(3), reference(3), k, sigma, inf, nan, pvalue
   real(dp) :: guard_draws(10), guard_smoothed(10)
   real(dp) :: shapes(5), p, q, pit_ties, invalid(5), original(5)
   type(pareto_diagnostic) :: diags
   checks = 0
   open(newunit = unit, file = 'test/reference_pareto.dat', status = 'old', action = 'read')
   read(unit, *) cases
   do c = 1, cases
      read(unit, '(a)') label
      read(unit, *) n
      allocate(x(n), weights(n), expected_x(n), smoothed(n))
      read(unit, *) x
      read(unit, *) weights
      read(unit, *) expected
      call gpdfit(x - minval(x) + 0.01_dp, actual(1), actual(2), status = status)
      if (status /= 0) error stop 'Expected GPD fit to succeed'
      call gpdfit(x - minval(x) + 0.01_dp, actual(3), actual(4), wip = .false.)
      call gpdfit(x - minval(x) + 0.01_dp, actual(5), actual(6), weights = weights)
      do s = 1, 4
         read(unit, *) expected_x
         if (s <= 3) then
            call pareto_smooth(x, smoothed, actual(s + 6), tail = trim(sides(s)), r_eff = 0.73_dp)
         else
            call pareto_smooth(x, smoothed, actual(s + 6), r_eff = 0.73_dp, are_log_weights = .true.)
         end if
         do j = 1, n
            call check(smoothed(j), expected_x(j), 2.0e-10_dp, trim(label) // ' smoothing', j)
         end do
      end do
      diags = pareto_diags(x, r_eff = 0.73_dp)
      actual(11:14) = [diags%khat, diags%min_ss, diags%khat_threshold, diags%convergence_rate]
      do j = 1, 14
         call check(actual(j), expected(j), 2.0e-9_dp, trim(label) // ' fit', j)
      end do
      read(unit, *) targets
      do s = 1, 4
         read(unit, *) reference
         do j = 1, 3
            select case (s)
            case (1)
               results(j) = pit(x, targets(j))
            case (2)
               results(j) = pit(x, targets(j), weights = weights)
            case (3)
               results(j) = pareto_pit(x, targets(j))
            case (4)
               results(j) = pareto_pit(x, targets(j), weights = log(weights), log_weights = .true.)
            end select
            call check(results(j), reference(j), 2.0e-10_dp, trim(label) // ' PIT', j)
         end do
      end do
      deallocate(x, weights, expected_x, smoothed)
   end do
   close(unit)
   open(newunit = unit, file = 'test/reference_uniformity.dat', status = 'old', action = 'read')
   read(unit, *) cases
   do c = 1, cases
      read(unit, '(a)') label
      read(unit, *) n
      allocate(x(n), pointwise(n), expected_test(n + 1))
      read(unit, *) x
      read(unit, *) expected_test
      call uniformity_test(x, trim(label), pvalue, pointwise, status)
      if (status /= 0) error stop 'Uniformity test unexpectedly unavailable'
      call check(pvalue, expected_test(1), 2.0e-9_dp, label, 0)
      do j = 1, n
         call check(pointwise(j), expected_test(j + 1), 2.0e-9_dp, label, j)
      end do
      deallocate(x, pointwise, expected_test)
   end do
   close(unit)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   shapes = [-0.5_dp, 0.0_dp, 0.00001_dp, 0.5_dp, 1.2_dp]
   do s = 1, size(shapes)
      do j = 0, 10
         p = real(j, dp) / 10.0_dp
         q = qgeneralized_pareto(p, mu = 2.0_dp, sigma = 3.0_dp, k = shapes(s))
         call check(pgeneralized_pareto(q, mu = 2.0_dp, sigma = 3.0_dp, k = shapes(s)), &
            p, 1.0e-11_dp, 'GPD inversion property', j)
      end do
   end do
   call check(qgeneralized_pareto(1.0_dp, k = -0.5_dp), 2.0_dp, 0.0_dp, 'GPD finite endpoint', 1)
   call check(qgeneralized_pareto(1.0_dp), inf, 0.0_dp, 'GPD infinite endpoint', 1)
   call check(pgeneralized_pareto(-inf), 0.0_dp, 0.0_dp, 'GPD left endpoint', 1)
   call check(pgeneralized_pareto(inf), 1.0_dp, 0.0_dp, 'GPD right endpoint', 1)
   call check(ps_min_ss(0.5_dp), 100.0_dp, 1.0e-14_dp, 'minimum sample size', 1)
   call check(ps_min_ss(1.0_dp), inf, 0.0_dp, 'minimum sample size divergent', 1)
   call check(ps_khat_threshold(100), 0.5_dp, 0.0_dp, 'k threshold', 1)
   if (ps_tail_length(1000, 1.0_dp) /= 94) error stop 'Tail length property failed'
   if (ps_tail_length(1000, -1.0_dp) /= -1) error stop 'Invalid efficiency guard failed'
   guard_draws = [(real(j, dp), j = 1, 10)]
   call pareto_smooth(guard_draws, guard_smoothed, k, tail = 'right_extra', r_eff = 1.0_dp, status = status)
   if (status == 0 .or. .not. ieee_is_nan(k)) error stop 'Truncated invalid tail was accepted'
   original = 1.0_dp
   call ps_tail(original, 4, invalid, k, status = status)
   if (status == 0 .or. .not. ieee_is_nan(k)) error stop 'Invalid short tail guard failed'
   if (any(invalid /= original)) error stop 'Invalid tail mutated its data'
   call gpdfit(original, k, sigma, status = status)
   if (status == 0 .or. .not. ieee_is_nan(k)) error stop 'Constant GPD guard failed'
   pit_ties = pit([1.0_dp, 2.0_dp, 2.0_dp, 4.0_dp], 2.0_dp, uniform = 0.25_dp)
   call check(pit_ties, 0.375_dp, 0.0_dp, 'supplied PIT randomizer', 1)
   call check(entropy([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp]), 1.0_dp, 1.0e-14_dp, 'entropy', 1)
   call check(dissent([1.0_dp, 1.0_dp, 3.0_dp, 3.0_dp]), 1.0_dp, 1.0e-14_dp, 'dissent', 1)
   call check(modal_category([3.0_dp, 1.0_dp, 3.0_dp, 1.0_dp]), 3.0_dp, 0.0_dp, 'mode tie', 1)
   print '(a,i0,a)', 'posterior Pareto/predictive: ', checks, ' R-oracle and analytic checks passed'
contains
   subroutine check(actual, expected, tolerance, context, index)
      real(dp), intent(in) :: actual !! Observed numerical result.
      real(dp), intent(in) :: expected !! Reference R result or analytic expectation.
      real(dp), intent(in) :: tolerance !! Relative tolerance with unit absolute floor.
      character(*), intent(in) :: context !! Human-readable test description.
      integer, intent(in) :: index !! Element number in the test case.
      logical :: ok
      if (ieee_is_nan(expected)) then
         ok = ieee_is_nan(actual)
      else if (.not. ieee_is_finite(expected)) then
         ok = actual == expected
      else
         ok = ieee_is_finite(actual)
         if (ok) ok = abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected))
      end if
      if (.not. ok) then
         print *, context, index, ' actual=', actual, ' expected=', expected
         error stop 'Numerical comparison failed'
      end if
      checks = checks + 1
   end subroutine check
end program test_pareto
