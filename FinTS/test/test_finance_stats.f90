! SPDX-License-Identifier: GPL-2.0-or-later
program test_finance_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf
   use fints
   implicit none

   type(summary_result) :: stats
   type(yearmon_result) :: months
   real(dp) :: x(4), gross, net, continuous
   complex(dp), allocatable :: conjugates(:)
   complex(dp) :: roots(6)

   x = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
   call FinTS_stats(x, stats, start=1990.01_dp)
   call check(stats%status == fints_ok, 'summary status')
   call check(stats%size == 4, 'summary size')
   call check_close(stats%mean, 2.5_dp, 1.0e-13_dp, 'summary mean')
   call check_close(stats%standard_deviation, sqrt(5.0_dp / 3.0_dp), 1.0e-13_dp, 'summary sd')
   call check_close(stats%skewness, 0.0_dp, 1.0e-13_dp, 'summary skewness')
   call check_close(stats%excess_kurtosis, -1.36_dp, 1.0e-13_dp, 'summary kurtosis')
   call check_close(stats%minimum, 1.0_dp, 1.0e-13_dp, 'summary minimum')
   call check_close(stats%maximum, 4.0_dp, 1.0e-13_dp, 'summary maximum')

   gross = compoundInterest(0.10_dp, periods=1.0_dp, frequency=2.0_dp)
   net = compoundInterest(0.10_dp, periods=1.0_dp, frequency=2.0_dp, net_value=.true.)
   continuous = compoundInterest(0.0446_dp, frequency=ieee_value(0.0_dp, ieee_positive_inf))
   call check_close(gross, 1.1025_dp, 1.0e-13_dp, 'compound gross')
   call check_close(net, 0.1025_dp, 1.0e-13_dp, 'compound net increase')
   call check_close(continuous, exp(0.0446_dp), 1.0e-13_dp, 'continuous compounding')
   call check_close(simple2logReturns(0.0456_dp), log(1.0456_dp), 1.0e-14_dp, 'simple to log')

   call as_yearmon2([2000.01_dp, 2000.02_dp, 2001.12_dp], months)
   call check(months%status == fints_ok .and. months%converted, 'yearmon conversion')
   call check(all(months%year == [2000, 2000, 2001]), 'yearmon years')
   call check(all(months%month == [1, 2, 12]), 'yearmon months')
   call as_yearmon2([200001, 200001, 200102], months)
   call check(.not. months%converted .and. months%duplicate_count == 1, 'yearmon duplicates')

   roots = [cmplx(1.0_dp, 1.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
      cmplx(1.0_dp, -1.0_dp, dp), cmplx(2.0_dp, -2.0_dp, dp), &
      cmplx(3.0_dp, 0.0_dp, dp), cmplx(2.0_dp, 2.0_dp, dp)]
   call findConjugates(roots, conjugates, 1.0e-12_dp)
   call check(size(conjugates) == 2, 'conjugate count')
   call check(any(abs(conjugates - cmplx(1.0_dp, -1.0_dp, dp)) < 1.0e-12_dp), 'first conjugate')
   call check(any(abs(conjugates - cmplx(2.0_dp, 2.0_dp, dp)) < 1.0e-12_dp), 'second conjugate')

   print '(a)', 'test_finance_stats: PASS'

contains

   subroutine check(condition, label)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: label
      if (.not. condition) then
         write(*, '(a,a)') 'FAIL: ', trim(label)
         error stop 1
      end if
   end subroutine check

   subroutine check_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual, expected, tolerance
      character(len=*), intent(in) :: label
      call check(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected)), label)
   end subroutine check_close

end program test_finance_stats
