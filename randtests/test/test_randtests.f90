! SPDX-License-Identifier: GPL-2.0-or-later
module test_callbacks
   use randtests, only: dp, randtests_aux
   implicit none
contains
   pure function squared_differences(x) result(value)
      real(dp), intent(in) :: x(:) !! One permutation of finite source values.
      real(dp) :: value
      value = sum((x(2:) - x(:size(x) - 1)) ** 2)
   end function squared_differences
end module test_callbacks

program test_randtests
   use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf, &
                                                                               ieee_negative_inf, ieee_is_nan, ieee_is_finite
   use randtests
   use r_reference, only: expected
   use test_callbacks, only: squared_differences
   implicit none
   real(dp), parameter :: inputs(8, 3) = reshape([ &
                                                 4.0_dp, 1.0_dp, 8.0_dp, 2.0_dp, 6.0_dp, 3.0_dp, 7.0_dp, 5.0_dp, &
                                                 1.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 5.0_dp, 4.0_dp, 4.0_dp, &
                                                 8.0_dp, 7.0_dp, 6.0_dp, 5.0_dp, 4.0_dp, 3.0_dp, 2.0_dp, 1.0_dp], [8, 3])
   character(len=11), parameter :: alternatives(3) = [character(len=11) :: &
                                                      'two.sided', 'left.sided', 'right.sided']
   type(randomness_result) :: test
   integer :: d, a, f, n, n1, n2, mask, r, previous, current, i, total, status, checks
   integer, allocatable :: histogram(:), brute(:), seed(:)
   integer :: random_a(100), random_b(100), seed_size
   real(dp) :: probability, cumulative, nan, inf, neginf, q, upper, lower
   real(dp), allocatable :: permutations(:, :), values(:), x(:)
   procedure(permutation_callback), pointer :: callback
   checks = 0
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   neginf = ieee_value(0.0_dp, ieee_negative_inf)

   ! Independently generated expected results from pinned R source.
   do d = 1, 3
      do a = 1, 3
         do f = 1, 8
            select case (f)
            case (1)
               test = difference_sign_test(inputs(:, d), trim(alternatives(a)))
            case (2)
               test = rank_test(inputs(:, d), trim(alternatives(a)))
            case (3)
               test = turning_point_test(inputs(:, d), trim(alternatives(a)))
            case (4)
               test = cox_stuart_test(inputs(:, d), trim(alternatives(a)))
            case (5)
               test = runs_test(inputs(:, d), trim(alternatives(a)))
            case (6)
               test = runs_test(inputs(:, d), trim(alternatives(a)), pvalue='exact')
            case (7)
               test = bartels_rank_test(inputs(:, d), trim(alternatives(a)))
            case (8)
               test = bartels_rank_test(inputs(:, d), trim(alternatives(a)), pvalue='beta')
            end select
            call require(test % status == rt_success, 'R parity status')
            call close (test % statistic, expected(1, f, a, d), 'R parity statistic')
            call close (test % p_value, expected(2, f, a, d), 'R parity p-value')
            call close (real(test % n, dp), expected(3, f, a, d), 'R parity sample size')
         end do
      end do
   end do

   ! Enumerate every binary sequence through length 10 for every class split.
   do n = 2, 10
      allocate (histogram(0:n))
      do n1 = 1, n - 1
         n2 = n - n1
         histogram = 0
         do mask = 0, 2 ** n - 1
            if (popcnt(mask) /= n1) cycle
            r = 1
            previous = ibits(mask, 0, 1)
            do i = 1, n - 1
               current = ibits(mask, i, 1)
               if (current /= previous) r = r + 1
               previous = current
            end do
            histogram(r) = histogram(r) + 1
         end do
         total = sum(histogram)
         cumulative = 0.0_dp
         do r = 0, n + 1
            probability = 0.0_dp
            if (r <= n) probability = real(histogram(r), dp) / real(total, dp)
            cumulative = cumulative + probability
            call close (druns(real(r, dp), n1, n2), probability, 'enumerated runs PMF')
            if (probability > 0.0_dp) then
               if (r >= 2) then
                  call close (qruns(pruns(real(r, dp), n1, n2), n1, n2), real(r, dp), 'self-CDF quantile')
                  call close (qruns(pruns(real(r, dp), n1, n2, .false.), n1, n2, .false.), &
                              real(r, dp), 'self-upper-CDF quantile')
               end if
               call close (druns(real(r, dp), n1, n2, .true.), log(probability), 'runs log PMF')
            end if
            call close (pruns(real(r, dp), n1, n2), cumulative, 'enumerated runs CDF')
            call close (pruns(real(r, dp) + 0.8_dp, n1, n2), cumulative, 'runs CDF floors q')
            call close (pruns(real(r, dp), n1, n2, .false.), 1.0_dp - cumulative, 'runs upper tail')
            if (cumulative > 0.0_dp) call close (pruns(real(r, dp), n1, n2, log_probability=.true.), &
                                                 log(cumulative), 'corrected runs log CDF')
         end do
         do i = 1, 99
            probability = real(i, dp) / 100.0_dp
            q = qruns(probability, n1, n2)
            call require(pruns(q, n1, n2) >= probability - 2.0e-13_dp, 'quantile upper inequality')
            call require(pruns(q - 1.0_dp, n1, n2) < probability + 2.0e-13_dp, 'quantile lower inequality')
            ! Non-boundary p makes equality robust against rounding in equivalent representations.
            if (abs(pruns(q, n1, n2) - probability) > 1.0e-12_dp .and. &
                abs(pruns(q - 1.0_dp, n1, n2) - probability) > 1.0e-12_dp) then
               call close (qruns(log(probability), n1, n2, log_probability=.true.), q, 'log quantile')
               call close (qruns(1.0_dp - probability, n1, n2, .false.), q, 'upper-tail quantile')
            end if
         end do
      end do
      deallocate (histogram)
   end do

   ! Compare streaming Bartels enumeration with the public permutation matrix.
   callback => squared_differences
   do n = 1, 7
      x = [(real(i, dp), i=1, n)]
      call bartels_histogram(n, histogram, status)
      call require(status == 0, 'Bartels histogram status')
      allocate (brute(0:ubound(histogram, 1)))
      brute = 0
      call permut_apply(x, callback, values, status)
      call require(status == 0, 'permutation callback status')
      do i = 1, size(values)
         r = nint(values(i))
         brute(r) = brute(r) + 1
      end do
      call require(all(histogram == brute), 'independent permutation Bartels histogram')
      do i = 1, 4
         q = real(i * ubound(histogram, 1) / 4, dp)
         lower = real(sum(brute(:int(q))), dp) / real(sum(brute), dp)
         call close (pbartelsrank(q, n), lower, 'Bartels exact CDF')
         call close (pbartelsrank(q, n, .false.), 1.0_dp - lower, 'Bartels exact upper tail')
         call close (dbartelsrank(q, n), real(brute(int(q)), dp) / real(sum(brute), dp), 'Bartels exact mass')
      end do
      deallocate (brute)
   end do
   call bartels_histogram(9, histogram, status)
   call require(status == 0 .and. sum(histogram) == 362880, 'Bartels n=9 cap supported')
   call bartels_histogram(10, histogram, status)
   call require(status /= 0 .and. size(histogram) == 0, 'Bartels n=10 rejected')

   ! Exact Bartels bugs: actual NM and inclusive right tail, small-n callbacks.
   test = bartels_rank_test([1.0_dp, 2.0_dp, 3.0_dp], pvalue='exact')
   call close (test % observed, 2.0_dp, 'corrected Bartels NM')
   call close (test % p_value, 2.0_dp / 3.0_dp, 'corrected exact Bartels two-sided')
   test = bartels_rank_test([1.0_dp, 2.0_dp, 3.0_dp], 'right.sided', 'exact')
   call close (test % p_value, 1.0_dp, 'corrected exact Bartels right tail')
   call close (dbartelsrank(1.0_dp, 2), 1.0_dp, 'corrected n=2 Bartels mass')
   call close (dbartelsrank(0.0_dp, 1), 1.0_dp, 'corrected n=1 Bartels mass')
   test = cox_stuart_test([1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp])
   call close (test % p_value, 1.0_dp, 'Cox all tied p-value capped')
   test = runs_test([0.0_dp, 1.0_dp], pvalue='exact')
   call close (test % p_value, 1.0_dp, 'degenerate exact runs has p=1')
   call require(test % status == rt_success, 'degenerate exact runs status')

   ! API and invalid-input coverage, including infinities without integer conversion.
   call close (druns(2.3_dp, 2, 3), 0.0_dp, 'noninteger mass is zero')
   call close (druns(inf, 2, 3), 0.0_dp, 'infinite mass is zero')
   call close (pruns(inf, 2, 3), 1.0_dp, 'positive infinite CDF')
   call close (pruns(neginf, 2, 3), 0.0_dp, 'negative infinite CDF')
   call require(.not. ieee_is_finite(pruns(0.0_dp, 2, 3, log_probability=.true.)), 'zero log CDF is -Inf')
   call close (qruns(0.0_dp, 2, 3), 2.0_dp, 'zero lower quantile')
   call close (qruns(1.0_dp, 2, 3), 5.0_dp, 'unit lower quantile')
   call close (qruns(0.0_dp, 2, 3, .false.), 5.0_dp, 'zero upper quantile')
   call close (qruns(1.0_dp, 2, 3, .false.), 2.0_dp, 'unit upper quantile')
   call close (qruns(neginf, 2, 3, log_probability=.true.), 2.0_dp, 'minus infinite log quantile')
   call require(ieee_is_nan(qruns(2.0_dp, 2, 3)), 'invalid probability rejected')
   call require(ieee_is_nan(qruns(nan, 2, 3)), 'NaN probability rejected')
   call require(ieee_is_nan(druns(2.0_dp, 0, 3)), 'zero count rejected')
   call require(ieee_is_nan(pruns(2.0_dp, huge(1), huge(1))), 'count overflow prevented')
   call require(ieee_is_nan(pbartelsrank(nan, 3)), 'NaN Bartels cutoff rejected')
   call require(ieee_is_nan(dbartelsrank(2.0_dp, 0)), 'invalid Bartels size rejected')
   test = difference_sign_test([1.0_dp, nan, 3.0_dp])
   call require(test % status == 0 .and. test % n == 2, 'NaNs omitted')
   test = turning_point_test([1.0_dp, 1.0_dp, 1.0_dp])
   call require(test % status == rt_insufficient, 'collapsed turning-point sample checked')
   test = rank_test([1.0_dp, inf])
   call require(test % status == rt_invalid, 'infinite data rejected')
   test = cox_stuart_test([1.0_dp, 2.0_dp], 'nonsense')
   call require(test % status == rt_invalid, 'alternative validated')
   test = bartels_rank_test([1.0_dp, 1.0_dp, 1.0_dp])
   call require(test % status == rt_degenerate, 'constant Bartels rejected')
   test = bartels_rank_test([1.0_dp, 1.0_dp, 2.0_dp], pvalue='exact')
   call require(test % status == rt_limit, 'exact Bartels ties rejected')
   test = runs_test([1.0_dp, 2.0_dp, 3.0_dp], threshold=inf)
   call require(test % status == rt_invalid, 'infinite threshold rejected')
   test = runs_test([1.0_dp, 2.0_dp, 3.0_dp], threshold=0.0_dp)
   call require(test % status == rt_degenerate, 'one-class runs rejected')
   test = runs_test([1.0_dp, 2.0_dp, 3.0_dp], pvalue='typo')
   call require(test % status == rt_invalid, 'p-value method validated')
   test = rank_test([real(dp) ::])
   call require(test % status == rt_insufficient .and. ieee_is_nan(test % p_value), 'empty input checked')
   test = difference_sign_test([2.0_dp, 2.0_dp])
   call close (test % p_value, 1.0_dp, 'constant difference-sign sequence')

   call permut([1.0_dp, 2.0_dp, 3.0_dp], permutations, status, 2)
   call require(status == 0 .and. all(shape(permutations) == [6, 2]), 'partial permutations shape')
   call require(all(permutations(1, :) == [1.0_dp, 2.0_dp]), 'permutation order first row')
   call require(all(permutations(6, :) == [3.0_dp, 2.0_dp]), 'permutation order last row')
   call permut([1.0_dp, 2.0_dp], permutations, status, 0)
   call require(status == 0 .and. all(shape(permutations) == [1, 0]), 'empty selection')
   call permut([1.0_dp, 2.0_dp], permutations, status, 3)
   call require(status == 1, 'invalid selection length')
   call permut([(real(i, dp), i=1, 11)], permutations, status)
   call require(status == 2, 'permutation resource bound')
   call permut([nan], permutations, status)
   call require(status == 1, 'permutation NaN rejected')
   nullify (callback)
   call permut_apply([1.0_dp], callback, values, status)
   call require(status == 1, 'null callback rejected')
   call close (randtests_aux([1.0_dp, 3.0_dp, 2.0_dp], 'difference.sign'), 1.0_dp, 'auxiliary rises')
   call close (randtests_aux([1.0_dp, 3.0_dp, 2.0_dp], 'bartels'), 5.0_dp, 'auxiliary differences')
   call close (randtests_aux([real(dp) ::], 'bartels'), 0.0_dp, 'empty auxiliary score')
   call require(ieee_is_nan(randtests_aux([1.0_dp], 'bad')), 'unknown auxiliary method')

   call close(randtests_aux([1.0e12_dp + 1.0_dp, 1.0e12_dp + 2.0_dp], 'bartels'), &
      1.0_dp, 'large offset preserves unit difference')
   call close(randtests_aux([1.0e16_dp + 2.0_dp, 1.0e16_dp + 4.0_dp], 'bartels'), &
      4.0_dp, 'large offset preserves representable difference')
   call require(ieee_is_nan(randtests_aux([huge(1.0_dp), -huge(1.0_dp)], 'bartels')), &
      'unrepresentable auxiliary statistic rejected')
   call close(randtests_aux([real(dp) ::], 'difference.sign'), 0.0_dp, 'empty auxiliary sign count')

   ! Reproducibility for the same Fortran runtime, without claiming R RNG parity.
   call random_seed(size=seed_size)
   allocate (seed(seed_size))
   seed = [(100 + 37 * i, i=1, seed_size)]
   call random_seed(put=seed)
   call rruns(random_a, 2, 3, status)
   call require(status == 0 .and. all(random_a >= 2) .and. all(random_a <= 5), 'random support')
   call random_seed(put=seed)
   call rruns(random_b, 2, 3)
   call require(all(random_a == random_b), 'same-runtime seed reproducibility')
   call rruns(random_b, 0, 3, status)
   call require(status == 1 .and. all(random_b == -1), 'invalid random parameters')
   print '(a,i0,a)', 'PASS: ', checks, ' checks, including 72 R parity cases and exhaustive discrete references.'
contains
   subroutine require(condition, label)
      logical, intent(in) :: condition !! Assertion that must be true for this test to pass.
      character(len=*), intent(in) :: label !! Name printed when the assertion fails.
      checks = checks + 1
      if (.not. condition) then
         print '(a)', 'FAIL: '//label
         error stop 1
      end if
   end subroutine require

   subroutine close (actual, target, label)
      real(dp), intent(in) :: actual !! Computed finite result to compare with its independent reference.
      real(dp), intent(in) :: target !! Expected finite numerical value.
      character(len=*), intent(in) :: label !! Description of the numerical assertion.
      if (.not. ieee_is_finite(actual)) then
         print '(a)', 'FAIL nonfinite: '//label
         error stop 1
      end if
      if (abs(actual - target) > 3.0e-12_dp * max(1.0_dp, abs(target))) then
         print '(a,2es24.15)', 'FAIL: '//label, actual, target
         error stop 1
      end if
      checks = checks + 1
   end subroutine close
end program test_randtests
