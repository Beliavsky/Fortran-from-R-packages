! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_perm
   use, intrinsic :: iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_positive_inf, ieee_is_finite
   use perm_api
   implicit none
   type(perm_result) :: result, other
   type(perm_control_type) :: options
   real(dp) :: scores(6), y(6), expected(30), interval(2), ties(40), levels(4), small(3)
   integer :: group(6), status, unit, i, j, hits(5), large_hits(4), cm_order(25), kgroup(6)
   integer, allocatable :: cm(:, :), indices(:)
   character(len=32) :: name, method
   character(len=16) :: alternatives(4)
   scores = [-2.3_dp, -2.3_dp, 0.0_dp, 1.7_dp, 2.1_dp, 5.2_dp]
   group = [1, 1, 0, 1, 0, 0]
   y = [2.0_dp, 8.0_dp, 4.0_dp, 3.0_dp, 9.0_dp, 5.0_dp]
   kgroup = [1, 1, 2, 2, 3, 3]
   alternatives = [character(len=16) :: 'less', 'greater', 'two.sided', 'two.sidedAbs']
   open(newunit=unit, file='reference/r_oracle.txt', status='old', action='read')
   read(unit, *) name, expected(:5)
   call calc_pvals([-4.0_dp, 0.0_dp, 0.0_dp, 1.0_dp, 7.0_dp], 1.0_dp, result, &
      weights=[1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 1.0_dp])
   call check_tails(result, expected(:5), 'R weighted tails')
   read(unit, *) name, expected(:5)
   call twosample_exact_ce(scores, group, result)
   call check_tails(result, expected(:5), 'R complete enumeration')
   read(unit, *) name, expected(:5)
   call twosample_exact_network(scores, group, result)
   call check_tails(result, expected(:5), 'R network')
   do i = 1, 4
      read(unit, *) name, expected(:7)
      call calc_pvals_mc([1.0_dp, 1.0_dp, 3.0_dp, 6.0_dp, 8.0_dp, 12.0_dp], 5.0_dp, result, &
         alternative=trim(alternatives(i)), confidence=0.95_dp)
      call check_tails(result, expected(:5), 'R simulated tails')
      call near(result%confidence_interval, expected(6:7), 2.0e-11_dp, 'R simulated confidence interval')
   end do
   read(unit, *) name, expected(:4)
   call calc_pvals_mc([1.0_dp, 1.0_dp, 3.0_dp, 6.0_dp, 8.0_dp, 12.0_dp], 5.0_dp, result, &
      two_sided_statistic=.true., confidence=0.95_dp)
   call near([result%p_twosided, result%p_equal, result%confidence_interval], expected(:4), &
      2.0e-11_dp, 'R upper-tail statistic')
   levels = [0.95_dp, 0.99_dp, 0.0_dp, 1.0_dp]
   hits = [0, 1, 5, 9, 10]
   do i = 1, 4
      do j = 1, 5
         read(unit, *) name, expected(:2)
         call p_ci(hits(j), 10, interval, status, levels(i))
         call require(status == perm_ok, 'pCI status')
         call near(interval, expected(:2), 2.0e-11_dp, 'R pCI endpoints')
      end do
   end do
   large_hits = [0, 1, 5000, 9999]
   do i = 1, 4
      read(unit, *) name, expected(:2)
      call p_ci(large_hits(i), 9999, interval, status)
      call near(interval, expected(:2), 2.0e-10_dp, 'R pCI large trials')
   end do
   read(unit, *) name, expected(:5)
   call twosample_pclt(scores, group, result)
   call require(result%status == perm_ok, 'TS CLT status')
   call near([result%p_twosided, result%p_absolute, result%p_less, result%p_greater, result%statistic], &
      expected(:5), 2.0e-12_dp, 'R two-sample CLT')
   read(unit, *) name, expected(:5)
   call trend_pclt(scores, y, result)
   call near([result%p_twosided, result%p_absolute, result%p_less, result%p_greater, result%statistic], &
      expected(:5), 2.0e-12_dp, 'R trend CLT')
   read(unit, *) name, expected(:3)
   call ksample_pclt(scores, kgroup, result)
   call near([result%p_twosided, result%statistic, real(result%df, dp)], expected(:3), 2.0e-12_dp, 'R k-sample CLT')
   ties(:10) = -0.25_dp
   ties(11:20) = 0.5_dp
   ties(21:30) = 1.25_dp
   ties(31:) = 2.75_dp
   read(unit, *) name, expected(:5)
   call twosample_exact_network(ties, [(modulo(i, 2), i=1,40)], result)
   call check_tails(result, expected(:5), 'R tied network')
   read(unit, *) name, cm_order
   call choose_matrix(5, 4, cm, status)
   call require(status == perm_ok, 'chooseMatrix status')
   call require(all(reshape(transpose(cm), [25]) == cm_order), 'R chooseMatrix row order')
   read(unit, *) name, expected(:5)
   call getcnt([1, 2, 1], [1, 3, 4], [2, 1], indices, status)
   call require(status == perm_ok, 'getcnt status')
   call near(real(indices, dp), expected(:5), 0.0_dp, 'R getcnt')
   close(unit)

   ! Wrapper estimates, alternatives and dispatch are separate from their low-level numerical kernels.
   call perm_ts(pack(scores, group == 1), pack(scores, group == 0), result, alternative='less', method='exact.ce')
   call near([result%p_value], [0.1_dp], 1.0e-14_dp, 'two-sample wrapper tail')
   call near([result%estimate], [-3.4_dp], 1.0e-14_dp, 'two-sample wrapper estimate')
   call perm_control(options, status, tsmethod='abs')
   call require(status == perm_ok, 'control construction')
   call perm_ts([0.0_dp], [1.0_dp, 4.0_dp], result, control=options, method='exact.network')
   call require(result%p_value == result%p_absolute, 'absolute wrapper convention')
   call perm_ks(scores, kgroup, result, method='pclt')
   call ksample_pclt(scores, kgroup, other)
   call near([result%p_value], [other%p_value], 0.0_dp, 'k-sample wrapper')
   call perm_trend(scores, y, result, alternative='greater', method='pclt')
   call trend_pclt(scores, y, other)
   call near([result%p_value, result%estimate], [other%p_greater, other%estimate], 0.0_dp, 'trend wrapper')
   call method_rule_ts1(scores, group, method, status)
   call require(status == perm_ok .and. method == 'exact.network', 'TS automatic rule')
   call method_rule_ts1(scores, group, method, status, exact=0)
   call require(status == perm_ok .and. method == 'pclt', 'TS false exact rule')
   call method_rule_ts1(scores, group, method, status, exact=1, bounds=[1])
   call require(status == perm_ok .and. method == 'exact.mc', 'TS exact MC fallback rule')
   call method_rule_ts1([1.0_dp, 1.0_dp, 2.0_dp, 2.0_dp], [0, 0, 1, 1], method, status)
   call require(status == perm_ok .and. method == 'exact.network', 'constant-group rule repair')
   call method_rule_ks1(scores, kgroup, method, status)
   call require(status == perm_ok .and. method == 'exact.mc', 'KS automatic rule')
   call method_rule_ks1(scores, kgroup, method, status, bound=1)
   call require(status == perm_ok .and. method == 'pclt', 'KS threshold rule')
   call method_rule_trend1(scores, y, method, status)
   call require(status == perm_ok .and. method == 'exact.mc', 'trend automatic rule')
   call method_rule_trend1(scores, y, method, status, bound=3)
   call require(status == perm_ok .and. method == 'pclt', 'trend threshold rule')

   ! Both sides of the asymmetric distribution are distinct and tested independently.
   call calc_pvals([0.0_dp, 1.0_dp, 4.0_dp], 0.0_dp, result)
   call near([result%p_twosided, result%p_absolute], [2.0_dp / 3.0_dp, 2.0_dp / 3.0_dp], &
      1.0e-14_dp, 'skew tails')
   call calc_pvals([0.0_dp, 1.0_dp, 4.0_dp], 4.0_dp, result)
   call near([result%p_twosided, result%p_absolute], [2.0_dp / 3.0_dp, 1.0_dp / 3.0_dp], &
      1.0e-14_dp, 'two distinct two-sided conventions')
   small = [1.0e-310_dp, 2.0e-310_dp, 3.0e-310_dp]
   call calc_pvals(small, small(1), result)
   call near([result%p_less, result%p_equal], [1.0_dp / 3.0_dp, 1.0_dp / 3.0_dp], &
      1.0e-14_dp, 'subnormal significant digits')
   block
      real(dp) :: unit_subnormal, atoms(4)
      unit_subnormal = nearest(0.0_dp, 1.0_dp)
      atoms = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp] * unit_subnormal
      call calc_pvals(atoms, atoms(4), result)
      call near([result%p_absolute], [0.5_dp], 0.0_dp, 'minimum-subnormal weighted centering')
      call calc_pvals_mc(atoms, atoms(4), result, alternative='two.sidedAbs')
      call near([result%p_absolute], [0.6_dp], 0.0_dp, 'minimum-subnormal MC centering')
   end block
   call p_ci(1, 10, interval, status, nearest(1.0_dp, -1.0_dp))
   call near([interval(1)], [5.5511151231257830e-18_dp], 1.0e-31_dp, 'near-one confidence tiny beta bound')
   call twosample_exact_ce([0.0_dp, 2.0_dp, -3.0_dp, -2.0_dp, -7.0_dp, 8.0_dp, -5.0_dp], &
      [8, 8, 8, 8, -3, -3, -3], result)
   call near([result%p_absolute], [33.0_dp / 35.0_dp], 1.0e-14_dp, 'decimal-decade carry regression')
   call choose_matrix(6, 3, cm, status)
   call perm_control(options, status, cm=cm)
   call require(status == perm_ok, 'supplied selection matrix control')
   call twosample_exact_ce(scores, group, result, options)
   call twosample_exact_ce(scores, group, other)
   call near([result%p_less, result%p_equal], [other%p_less, other%p_equal], 0.0_dp, 'supplied selection matrix')
   call large_tie_checks()
   call exhaustive_checks()
   call random_input_checks()
   call monte_carlo_checks()
   call invalid_checks()
   print '(a)', 'perm: R fixtures, exhaustive subsets, random-input checks and Monte Carlo calibration passed'
contains
   subroutine require(condition, message)
      logical, intent(in) :: condition !! Assertion that must be true for the test to pass.
      character(len=*), intent(in) :: message !! Human-readable failure label.
      if (.not. condition) then
         print '(a)', 'FAIL: ' // message
         error stop 1
      end if
   end subroutine require

   subroutine near(actual, wanted, tolerance, message)
      real(dp), intent(in) :: actual(:) !! Computed numerical values.
      real(dp), intent(in) :: wanted(:) !! Independently generated expected values of the same shape.
      real(dp), intent(in) :: tolerance !! Absolute numerical comparison tolerance.
      character(len=*), intent(in) :: message !! Human-readable failure label.
      call require(size(actual) == size(wanted), message // ' shape')
      if (any(.not. ieee_is_finite(actual)) .or. any(abs(actual - wanted) > tolerance)) then
         print '(a)', 'FAIL: ' // message
         print *, 'actual: ', actual
         print *, 'wanted: ', wanted
         error stop 1
      end if
   end subroutine near

   subroutine check_tails(r, wanted, message)
      type(perm_result), intent(in) :: r !! Result whose status and five tail probabilities are checked.
      real(dp), intent(in) :: wanted(5) !! Expected central, absolute, lower, upper and equal probabilities.
      character(len=*), intent(in) :: message !! Human-readable check label.
      call require(r%status == perm_ok, message // ' status')
      call near([r%p_twosided, r%p_absolute, r%p_less, r%p_greater, r%p_equal], wanted, 2.0e-11_dp, message)
   end subroutine check_tails

   subroutine large_tie_checks()
      real(dp) :: x(2000)
      integer :: g(2000), k
      type(perm_result) :: r
      x(:1000) = 0.0_dp
      x(1001:) = 1.0_dp
      g = [(modulo(k, 2), k=1,2000)]
      call twosample_exact_network(x, g, r)
      call require(r%status == perm_ok, 'large tied log-weight distribution')
      ! Python arbitrary-precision integers independently evaluate choose(1000,500)^2/choose(2000,1000).
      call near([r%p_equal, r%p_less], [0.035669103903581877_dp, 0.51783455195179096_dp], &
         1.0e-10_dp, 'large binomial multiplicity overflow avoidance')
   end subroutine large_tie_checks

   subroutine exhaustive_checks()
      real(dp), allocatable :: x(:), distribution(:)
      integer, allocatable :: g(:)
      real(dp) :: observed, mean_value, less, greater, equal, absolute
      integer :: n, m, mask, other_mask, count_values, k, trial
      type(perm_result) :: net, ce
      do n = 3, 8
         allocate(x(n), g(n), distribution(2 ** n))
         do trial = 1, 3
            do k = 1, n
               x(k) = real(modulo(k * k + trial * k, 5), dp) - 2.0_dp
            end do
            do mask = 1, 2 ** n - 2
               m = popcnt(mask)
               g = [(merge(1, 0, btest(mask, k - 1)), k=1,n)]
               observed = sum(x * real(g, dp))
               count_values = 0
               do other_mask = 0, 2 ** n - 1
                  if (popcnt(other_mask) /= m) cycle
                  count_values = count_values + 1
                  distribution(count_values) = 0.0_dp
                  do k = 1, n
                     if (btest(other_mask, k - 1)) distribution(count_values) = distribution(count_values) + x(k)
                  end do
               end do
               mean_value = sum(x) * real(m, dp) / real(n, dp)
               less = real(count(distribution(:count_values) <= observed), dp) / real(count_values, dp)
               greater = real(count(distribution(:count_values) >= observed), dp) / real(count_values, dp)
               equal = real(count(distribution(:count_values) == observed), dp) / real(count_values, dp)
               absolute = real(count(abs(distribution(:count_values) - mean_value) >= &
                  abs(observed - mean_value) - 1.0e-12_dp), dp) / real(count_values, dp)
               call twosample_exact_network(x, g, net)
               call check_tails(net, [min(1.0_dp, 2.0_dp * min(less, greater)), absolute, less, greater, equal], &
                  'independent exhaustive network')
               call twosample_exact_ce(x, g, ce)
               call near([ce%p_twosided, ce%p_absolute, ce%p_equal], &
                  [net%p_twosided, net%p_absolute, net%p_equal], 2.0e-11_dp, 'exhaustive CE/network agreement')
            end do
         end do
         deallocate(x, g, distribution)
      end do
   end subroutine exhaustive_checks

   subroutine random_input_checks()
      integer(int64) :: state
      real(dp) :: x(8)
      integer :: g(8), trial, k
      type(perm_result) :: net, ce, shifted, swapped
      state = 271828_int64
      g = [1, 1, 1, 0, 0, 0, 0, 0]
      do trial = 1, 100
         do k = 1, 8
            state = modulo(48271_int64 * state, 2147483647_int64)
            x(k) = real(modulo(state, 2001_int64) - 1000_int64, dp) / 8.0_dp
         end do
         call twosample_exact_ce(x, g, ce)
         call twosample_exact_network(x, g, net)
         call near([net%p_twosided, net%p_absolute, net%p_less, net%p_equal], &
            [ce%p_twosided, ce%p_absolute, ce%p_less, ce%p_equal], 2.0e-12_dp, 'random real-score CE/network')
         call twosample_exact_network(8.0_dp * x + 1000.0_dp, g, shifted)
         call near([shifted%p_twosided, shifted%p_less], [net%p_twosided, net%p_less], &
            2.0e-12_dp, 'positive-scale location invariance')
         call twosample_exact_network(x, 1 - g, swapped)
         call near([swapped%p_twosided, swapped%p_greater], [net%p_twosided, net%p_less], &
            2.0e-12_dp, 'sample swap tail identity')
      end do
   end subroutine random_input_checks

   subroutine monte_carlo_checks()
      type(perm_control_type) :: opt
      type(perm_result) :: first, second, exact_result
      type(perm_rng) :: stream
      real(dp) :: x(6), trend_values(720), permuted(6), expected_mc, total_estimate
      integer :: g(6), status_local, trial, covered, p(6), nvalues
      x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 9.0_dp]
      g = [1, 1, 0, 0, 0, 0]
      call perm_control(opt, status_local, nmc=1999, confidence=0.99_dp)
      call twosample_exact_network(x, g, exact_result)
      total_estimate = 0.0_dp
      covered = 0
      do trial = 1, 100
         opt%seed = 101_int64 + int(trial, int64) * 982451653_int64
         call twosample_exact_mc(x, g, first, opt, 'less')
         call require(first%status == perm_ok .and. first%nmc == 1999, 'fixed-count MC status')
         total_estimate = total_estimate + first%p_less
         if (first%confidence_interval(1) <= exact_result%p_less .and. &
            first%confidence_interval(2) >= exact_result%p_less) covered = covered + 1
      end do
      expected_mc = (1999.0_dp * exact_result%p_less + 1.0_dp) / 2000.0_dp
      call near([total_estimate / 100.0_dp], [expected_mc], 0.003_dp, 'MC expectation calibration')
      call require(covered >= 95, 'Monte Carlo 99 percent interval calibration')
      opt%seed = 98765_int64
      call twosample_exact_mc(x, g, first, opt)
      call twosample_exact_mc(x, g, second, opt)
      call near([first%p_twosided, first%confidence_interval], &
         [second%p_twosided, second%confidence_interval], 0.0_dp, 'deterministic RNG replay')
      stream%state = 456_int64
      opt%set_seed = .false.
      call twosample_exact_mc(x, g, first, opt, rng=stream)
      call require(stream%state /= 456_int64, 'explicit RNG state advanced')
      call twosample_exact_mc(x, g, second, opt)
      call require(second%status == perm_invalid, 'missing continuation stream rejected')
      opt%set_seed = .true.
      opt%nmc = 99999
      call ksample_exact_mc(x, [1, 1, 2, 2, 3, 3], first, opt)
      call require(first%status == perm_ok, 'k-sample MC status')
      call ksample_exact_mc(x * 1.0e-200_dp, [1, 1, 2, 2, 3, 3], second, opt)
      call near([second%p_value], [first%p_value], 0.0_dp, 'subnormal-square MC scale invariance')
      ! Independent six-factorial enumeration for the same quadratic k-sample statistic.
      nvalues = 0
      p = [1, 2, 3, 4, 5, 6]
      call enumerate(1, p, x, nvalues, trend_values)
      permuted = x - sum(x) / 6.0_dp
      expected_mc = sum([sum(permuted(:2)) ** 2, sum(permuted(3:4)) ** 2, sum(permuted(5:)) ** 2]) / 2.0_dp
      expected_mc = real(count(trend_values >= expected_mc - 1.0e-12_dp), dp) / 720.0_dp
      call near([first%p_value], [expected_mc], 0.007_dp, 'k-sample MC exhaustive calibration')
      call trend_exact_mc(x, [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp], first, opt, 'greater')
      call near([first%p_value], [1.0_dp / 720.0_dp], 0.001_dp, 'trend MC exhaustive calibration')
      call trend_exact_mc(x * 1.0e-200_dp, &
         [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp, 6.0_dp] * 1.0e-200_dp, second, opt, 'greater')
      call near([second%p_value], [first%p_value], 0.0_dp, 'underflow-product trend scale invariance')
      call calc_pvals_mc(spread(1.0_dp, 1, 1000), 1.0_dp, second)
      call near(second%confidence_interval, [1.0_dp, 1.0_dp], 0.0_dp, 'bounded all-ties central interval')
   end subroutine monte_carlo_checks

   recursive subroutine enumerate(position, p, x, nvalues, trend_values)
         integer, intent(in) :: position !! Current position in the independent permutation enumeration.
         integer, intent(inout) :: p(6) !! Current permutation of the six observation indices.
         real(dp), intent(in) :: x(6) !! Fixed data values used for the independent enumeration.
         integer, intent(inout) :: nvalues !! Count of complete permutations already evaluated.
         real(dp), intent(inout) :: trend_values(720) !! Storage for all six-factorial statistic values.
         integer :: j, temporary
         real(dp) :: centered(6)
         if (position == 7) then
            nvalues = nvalues + 1
            centered = x(p) - sum(x) / 6.0_dp
            trend_values(nvalues) = (sum(centered(:2)) ** 2 + sum(centered(3:4)) ** 2 + &
               sum(centered(5:)) ** 2) / 2.0_dp
            return
         end if
         do j = position, 6
            temporary = p(j)
            p(j) = p(position)
            p(position) = temporary
            call enumerate(position + 1, p, x, nvalues, trend_values)
            temporary = p(j)
            p(j) = p(position)
            p(position) = temporary
         end do
   end subroutine enumerate

   subroutine invalid_checks()
      type(perm_result) :: r
      type(perm_control_type) :: opt
      integer, allocatable :: matrix(:, :), cnt(:)
      integer :: st
      real(dp) :: ci(2), x(4)
      x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
      call choose_matrix(50, 25, matrix, st)
      call require(st == perm_limit .and. .not. allocated(matrix), 'choose overflow/allocation guard')
      call choose_matrix(4, 5, matrix, st)
      call require(st == perm_invalid, 'invalid combination size')
      call choose_matrix(0, 0, matrix, st)
      call require(st == perm_ok .and. size(matrix, 1) == 1 .and. size(matrix, 2) == 0, 'empty selection identity')
      call getcnt([2], [1, 2], [1], cnt, st)
      call require(st == perm_invalid, 'invalid network node')
      call getcnt([1], [1, 2], [2], cnt, st)
      call require(st == perm_invalid, 'inconsistent network edges')
      call perm_control(opt, st, nmc=0)
      call require(st == perm_invalid, 'invalid simulation count')
      call perm_control(opt, st, digits=16)
      call require(st == perm_invalid, 'unsupported significant digits')
      call perm_control(opt, st, confidence=ieee_value(0.0_dp, ieee_quiet_nan))
      call require(st == perm_invalid, 'NaN confidence')
      call p_ci(-1, 10, ci, st)
      call require(st == perm_invalid, 'negative tail count')
      call p_ci(0, 0, ci, st)
      call require(st == perm_invalid, 'empty Monte Carlo trials')
      call trend_pclt(x, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], r)
      call require(r%status == perm_degenerate, 'constant trend covariate')
      call ksample_pclt(x, [1, 1, 1, 1], r)
      call require(r%status == perm_invalid, 'one-group k-sample rejected')
      call twosample_exact_ce(x, [0, 0, 0, 0], r)
      call require(r%status == perm_invalid, 'one-group exact test rejected')
      call perm_control(opt, st, max_permutations=2_int64)
      call twosample_exact_network(x, [0, 0, 1, 1], r, opt)
      call require(r%status == perm_limit, 'network enumeration cap')
      call twosample_exact_network(x, [0, 0, 0, 1], r, opt)
      call require(r%status == perm_limit, 'singleton-network enumeration cap')
      call perm_control(opt, st, max_work=1_int64)
      call twosample_exact_mc(x, [0, 0, 1, 1], r, opt)
      call require(r%status == perm_limit, 'Monte Carlo operation cap')
      x(1) = ieee_value(0.0_dp, ieee_positive_inf)
      call twosample_exact_network(x, [0, 0, 1, 1], r)
      call require(r%status == perm_invalid, 'infinite score rejected')
      x = huge(1.0_dp) / 2.0_dp
      call twosample_exact_ce(x, [0, 0, 1, 1], r)
      call require(r%status == perm_numeric, 'score overflow guard')
      call calc_pvals([1.0_dp, 2.0_dp], 1.0_dp, r, weights=[0.0_dp, 0.0_dp])
      call require(r%status == perm_invalid, 'zero-weight distribution rejected')
      call calc_pvals_mc([1.0_dp], 1.0_dp, r, alternative='invalid')
      call require(r%status == perm_invalid, 'unknown alternative rejected')
      call perm_ts([1.0_dp], [2.0_dp], r, method='invalid')
      call require(r%status == perm_invalid, 'unknown test method rejected')
   end subroutine invalid_checks
end program test_perm
