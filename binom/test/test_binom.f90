! SPDX-License-Identifier: GPL-2.0-or-later
program test_binom
   use binom
   use binom_numerics, only: beta_quantile, beta_probability, binomial_deviance, binomial_mass
   use r_distributions, only: r_qnorm
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   type(binom_interval) :: ci, all_ci(11)
   type(binom_design_result) :: design
   type(binom_sim_result) :: simulation, from_counts
   character(len=32) :: method, alternative, interval_type
   real(dp) :: x, n, level, mean, lower, upper, p, expected, width, palt, alpha, phi, a, b, sig
   real(dp) :: u(20), q, finite_tail, nan, inf, residual, total
   integer :: unit, io, cases, ib, ia, i, j, captured, counts(20), rows
   cases = 0
   open(newunit=unit, file='test/data/intervals.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) method, x, n, level, mean, lower, upper
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      ci = binom_confint(x, n, trim(method), level)
      call check(ci%status == 0, 'reference interval status')
      call near(ci%mean, mean, 3.0e-13_dp, 'interval mean')
      call near(ci%lower, lower, 3.0e-7_dp, 'interval lower')
      call near(ci%upper, upper, 3.0e-7_dp, 'interval upper')
      cases = cases + 1
   end do
   close(unit)
   call check(rows == 660, 'intervals fixture row count')
   call check(cases >= 600, 'full interval fixture read')
   open(newunit=unit, file='test/data/summary.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) method, p, n, level, expected, width
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      ! Exact Wilson endpoint limits intentionally correct upstream ulp-sensitive endpoint coverage.
      if (trim(method) == 'wilson' .and. (p == 0.0_dp .or. p == 1.0_dp)) expected = 1.0_dp
      call near(binom_coverage(p, int(n), trim(method), level), expected, 3.0e-12_dp, 'R coverage')
      call near(binom_length(p, int(n), trim(method), level), width, 3.0e-7_dp, 'R expected width')
   end do
   close(unit)
   call check(rows == 66, 'summary fixture row count')
   open(newunit=unit, file='test/data/power.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) method, alternative, palt, n, p, alpha, phi, expected
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      q = binom_power(palt, n, p, alpha, phi, trim(alternative), trim(method))
      call near(q, expected, 3.0e-7_dp, 'R power')
   end do
   close(unit)
   call check(rows == 72, 'power fixture row count')
   open(newunit=unit, file='test/data/bayes.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) interval_type, x, n, level, a, b, mean, lower, upper, sig
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      ci = binom_bayes(x, n, level, trim(interval_type), a, b)
      call near(ci%mean, mean, 3.0e-13_dp, 'posterior mean')
      call near(ci%lower, lower, 3.0e-7_dp, 'posterior lower')
      call near(ci%upper, upper, 3.0e-7_dp, 'posterior upper')
      call near(ci%sig, sig, 3.0e-7_dp, 'posterior outside mass')
   end do
   close(unit)
   call check(rows == 8, 'bayes fixture row count')
   open(newunit=unit, file='test/data/lrt.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) x, n, level, ib, ia, lower, upper
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      ci = binom_lrt(x, n, level, bayes=ib == 1, conf_adj=ia == 1)
      call near(ci%lower, lower, 3.0e-8_dp, 'adjusted LRT lower')
      call near(ci%upper, upper, 3.0e-8_dp, 'adjusted LRT upper')
   end do
   close(unit)
   call check(rows == 20, 'lrt fixture row count')
   open(newunit=unit, file='test/data/profile.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) x, n, level, ib, lower, upper
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      ci = binom_profile(x, n, level, 100, 0.15_dp, ib == 1)
      call near(ci%lower, lower, 3.0e-10_dp, 'custom profile lower')
      call near(ci%upper, upper, 3.0e-10_dp, 'custom profile upper')
   end do
   close(unit)
   call check(rows == 10, 'profile fixture row count')
   open(newunit=unit, file='test/data/design.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) alternative, ib, palt, n, expected
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      design = cloglog_sample_size(p_alt=palt, p=0.4_dp, alternative=trim(alternative), &
         exact_n=ib == 1, recompute_power=.true., phi=1.3_dp)
      call check(design%status == 0, 'sample-size status')
      call near(design%n, n, 3.0e-10_dp, 'R sample size')
      call near(design%power, expected, 3.0e-12_dp, 'R attained power')
   end do
   close(unit)
   call check(rows == 6, 'design fixture row count')
   open(newunit=unit, file='test/data/detectable.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) alternative, expected
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      design = cloglog_sample_size(n=100.0_dp, p=0.4_dp, power=0.4_dp, alternative=trim(alternative), phi=1.3_dp)
      call check(design%status == 0, 'detectable-alternative status')
      call near(design%p_alt, expected, 1.0e-4_dp, 'R detectable alternative')
      call near(binom_power(design%p_alt, 100.0_dp, 0.4_dp, phi=1.3_dp, &
         alternative=trim(alternative)), 0.4_dp, 3.0e-12_dp, 'detectable alternative root')
   end do
   close(unit)
   call check(rows == 3, 'detectable fixture row count')
   open(newunit=unit, file='test/data/simulation_counts.csv', status='old', action='read')
   rows = 0
   do i = 1, 20
      read(unit, *) u(i), counts(i)
      rows = rows + 1
   end do
   close(unit)
   call check(rows == 20, 'simulation_counts fixture row count')
   open(newunit=unit, file='test/data/simulation.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) method, captured, mean, lower, upper, width
      call check(io <= 0, 'valid CSV row format')
      if (io < 0) exit
      rows = rows + 1
      simulation = binom_sim(u, 10, 0.4_dp, trim(method))
      from_counts = binom_sim_counts(counts, 10, 0.4_dp, trim(method))
      call check(simulation%status == 0, 'uniform simulation status')
      call check(simulation%captured == captured, 'uniform simulation capture count')
      call check(from_counts%captured == captured, 'supplied-count simulation capture count')
      call near(simulation%capture_interval%mean, mean, 3.0e-12_dp, 'capture mean')
      call near(simulation%capture_interval%lower, lower, 3.0e-7_dp, 'capture lower')
      call near(simulation%capture_interval%upper, upper, 3.0e-7_dp, 'capture upper')
      call near(simulation%width, width, 3.0e-7_dp, 'simulation mean width')
   end do
   close(unit)
   call check(rows == 11, 'simulation fixture row count')
   open(newunit=unit, file='test/data/beta.csv', status='old', action='read')
   rows = 0
   do
      read(unit, *, iostat=io) p, a, b, expected
      call check(io <= 0, 'valid beta CSV row format')
      if (io < 0) exit
      rows = rows + 1
      q = beta_quantile(p, a, b)
      call near(q / max(1.0e-300_dp, expected), expected / max(1.0e-300_dp, expected), &
         1.0e-8_dp, 'R beta quantile tail grid')
   end do
   close(unit)
   call check(rows == 180, 'beta fixture row count')
   ! Independent identities do not depend on R-generated fixture values.
   do i = 1, 9
      p = real(i, dp) / 10.0_dp
      call near(beta_quantile(p, 1.0_dp, 1.0_dp), p, 3.0e-14_dp, 'uniform beta inverse')
      q = beta_quantile(p, 0.5_dp, 0.5_dp)
      call near(q, sin(acos(-1.0_dp) * p / 2.0_dp) ** 2, 3.0e-14_dp, 'arcsine beta inverse')
      q = beta_quantile(p, 3.25_dp, 17.5_dp)
      call near(beta_probability(q, 3.25_dp, 17.5_dp), p, 3.0e-13_dp, 'beta inverse residual')
      total = 0.0_dp
      do j = 0, 100
         total = total + binomial_mass(real(j, dp), 100.0_dp, p)
      end do
      call near(total, 1.0_dp, 3.0e-13_dp, 'binomial mass normalization')
      call check(binom_coverage(p, 30, 'exact') >= 0.95_dp - 1.0e-12_dp, 'exact coverage guarantee')
   end do
   call near(beta_probability(0.0_dp, 0.0_dp, 1.0_dp), 1.0_dp, 0.0_dp, 'zero-shape beta atom at zero')
   call near(beta_probability(0.0_dp, 0.0_dp, 0.0_dp), 0.5_dp, 0.0_dp, 'two zero-shape beta endpoint')
   finite_tail = beta_quantile(1.0e-20_dp, 0.5_dp, 2.0_dp)
   call check(finite_tail > 0.0_dp, 'tiny beta quantile retained')
   call near(beta_probability(finite_tail, 0.5_dp, 2.0_dp) / 1.0e-20_dp, 1.0_dp, 2.0e-13_dp, 'tiny beta tail')
   do i = 1, 99
      ci = binom_lrt(real(i, dp), 100.0_dp, tol=epsilon(1.0_dp))
      residual = r_qnorm(0.025_dp) ** 2
      call near(binomial_deviance(ci%lower, real(i, dp), 100.0_dp), residual, 4.0e-12_dp, 'LRT lower root')
      call near(binomial_deviance(ci%upper, real(i, dp), 100.0_dp), residual, 4.0e-12_dp, 'LRT upper root')
   end do
   call near(binom_coverage(0.0_dp, 5, 'wilson'), 1.0_dp, 0.0_dp, 'Wilson exact zero-boundary coverage')
   call near(binom_coverage(1.0_dp, 20, 'wilson'), 1.0_dp, 0.0_dp, 'Wilson exact unit-boundary coverage')
   ci = binom_bayes(7.0_dp, 30.0_dp)
   call near(ci%sig, 0.05_dp, 3.0e-13_dp, 'HPD posterior mass')
   residual = (ci%shape1 - 1.0_dp) * log(ci%lower / ci%upper) + &
      (ci%shape2 - 1.0_dp) * log((1.0_dp - ci%lower) / (1.0_dp - ci%upper))
   call near(residual, 0.0_dp, 1.0e-10_dp, 'HPD equal endpoint density')
   design = cloglog_sample_size(p_alt=0.7_dp, n=1.0e-320_dp, phi=1.0e308_dp)
   call check(design%status /= binom_ok, 'nonfinite design power status')
   design = cloglog_sample_size(p_alt=0.500001_dp)
   call check(design%status == 0 .and. design%n > 2147483647.0_dp, 'large sample-size support')
   call near(design%n, aint(design%n), 0.0_dp, 'large sample-size upward rounding')
   design = cloglog_sample_size(n=100.0_dp, p=0.4_dp, alternative='greater', phi=1.3_dp)
   call check(design%status == binom_convergence, 'unbracketed sample-size inverse')
   all_ci = binom_confint_all(3.0_dp, 10.0_dp)
   call check(all(all_ci%status == 0), 'all interval methods')
   ci = binom_wilson(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(11)%lower, 1.0e-15_dp, 'Wilson wrapper')
   ci = binom_agresti_coull(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(1)%lower, 1.0e-15_dp, 'Agresti-Coull wrapper')
   ci = binom_asymp(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(2)%lower, 1.0e-15_dp, 'asymptotic wrapper')
   ci = binom_cloglog(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(4)%lower, 1.0e-15_dp, 'cloglog wrapper')
   ci = binom_exact(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(5)%lower, 1.0e-15_dp, 'exact wrapper')
   ci = binom_logit(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(6)%lower, 1.0e-15_dp, 'logit wrapper')
   ci = binom_probit(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(7)%lower, 1.0e-15_dp, 'probit wrapper')
   ci = binom_prop_test(3.0_dp, 10.0_dp)
   call near(ci%lower, all_ci(10)%lower, 1.0e-15_dp, 'prop.test wrapper')
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   ci = binom_confint(nan, 10.0_dp)
   call check(ci%status == binom_domain .and. ieee_is_nan(ci%lower), 'NaN success input')
   ci = binom_confint(1.0_dp, inf)
   call check(ci%status == binom_domain, 'infinite trial input')
   ci = binom_confint(-1.0_dp, 10.0_dp)
   call check(ci%status == binom_domain, 'negative success input')
   ci = binom_confint(1.0_dp, 0.0_dp)
   call check(ci%status == binom_domain, 'zero trials')
   ci = binom_confint(1.0_dp, 10.0_dp, conf_level=1.0_dp)
   call check(ci%status == binom_domain, 'unit confidence level')
   ci = binom_confint(1.0_dp, 10.0_dp, method='no-such-method')
   call check(ci%status == binom_bad_method, 'unknown interval method')
   ci = binom_profile(2.0_dp, 10.0_dp, maxsteps=5000)
   call check(ci%status == binom_resource, 'profile step resource cap')
   ci = binom_profile(2.0_dp, 10.0_dp, maxsteps=1, del=1.0e-6_dp)
   call check(ci%status == binom_convergence, 'profile target not bracketed')
   ci = binom_lrt(1.0_dp, 1.0_dp, bayes=.true., prior_shape1=1.0e308_dp, prior_shape2=1.0e308_dp)
   call check(ci%status == binom_resource, 'overflowing LRT pseudocounts rejected')
   ci = binom_bayes(2.0_dp, 10.0_dp, prior_shape1=-1.0_dp)
   call check(ci%status == binom_domain, 'invalid beta prior')
   ci = binom_bayes(2.0_dp, 10.0_dp, tol=1.0e-30_dp, maxit=1000)
   call check(ci%status == binom_convergence, 'unattainable HPD tolerance terminates')
   ci = binom_bayes(2.0_dp, 10.0_dp, maxit=1)
   call check(ci%status == binom_convergence, 'HPD iteration cap')
   call check(ieee_is_nan(binom_coverage(0.5_dp, binom_enumeration_limit + 1)), 'enumeration cap')
   call check(ieee_is_nan(binom_power(0.6_dp, p=0.0_dp)), 'transformed power endpoint domain')
   call near(binom_power(0.0_dp, method='asymp'), 1.0_dp, 0.0_dp, 'asymptotic zero alternative')
   call near(binom_power(1.0_dp, method='asymp'), 1.0_dp, 0.0_dp, 'asymptotic unit alternative')
   call near(binom_power(0.0_dp, p=0.0_dp, method='asymp'), 0.05_dp, 0.0_dp, 'asymptotic degenerate null')
   call near(binom_power(0.0_dp, method='exact'), 1.0_dp, 1.0e-14_dp, 'exact power limiting alternative')
   simulation = binom_sim([0.5_dp], binom_enumeration_limit + 1, 0.4_dp)
   call check(simulation%status == binom_resource, 'simulation resource cap')
   simulation = binom_sim([1.0_dp], 10, 0.4_dp)
   call check(simulation%status == binom_domain, 'invalid uniform endpoint')
   simulation = binom_sim([0.0_dp, 0.5_dp], 10, 0.0_dp, 'exact')
   call check(simulation%captured == 2, 'degenerate zero-probability simulation')
   simulation = binom_sim([0.0_dp, 0.5_dp], 10, 1.0_dp, 'exact')
   call check(simulation%captured == 2, 'degenerate unit-probability simulation')
   print '(a,i0,a)', 'binom: ', cases, ' R interval rows and design/coverage/simulation/identity tests passed'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Assertion that must be true for the named test to pass.
      character(len=*), intent(in) :: label !! Human-readable name reported on a failed assertion.
      if (.not. condition) then
         print *, 'FAIL: ', label
         error stop 1
      end if
   end subroutine check
   subroutine near(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed finite numerical result.
      real(dp), intent(in) :: expected !! Reference or independent mathematical result.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error.
      character(len=*), intent(in) :: label !! Human-readable test description.
      if (ieee_is_nan(actual) .or. abs(actual - expected) > tolerance) then
         print *, 'FAIL: ', label, ' actual=', actual, ' expected=', expected, ' tolerance=', tolerance
         print *, 'Current fixture context: ', trim(method), x, n, level
         error stop 1
      end if
   end subroutine near
end program test_binom
