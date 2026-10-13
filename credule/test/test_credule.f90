! SPDX-License-Identifier: MIT
program test_credule
   use credule
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   integer :: assertions = 0, reference_prices = 0, reference_bootstraps = 0

   include 'reference_calls.inc'
   call interpolation_tests()
   call analytic_tests()
   call inversion_tests()
   call edge_tests()
   print '(a,i0,a,i0,a,i0)', 'PASS: ', assertions, ' assertions; R pricing cases ', &
                           reference_prices, '; R bootstrap cases ', reference_bootstraps

contains

   subroutine require(condition, label)
      logical, intent(in) :: condition !! Assertion predicate that must be true for the deterministic test to pass.
      character(*), intent(in) :: label !! Diagnostic identifying the numerical behavior or R reference case.
      assertions = assertions + 1
      if (.not. condition) then
         print '(a)', 'FAIL: ' // label
         error stop 1
      end if
   end subroutine require

   subroutine near(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual(:) !! Computed values to compare against an independent or frozen reference.
      real(dp), intent(in) :: expected(:) !! Reference values, with the same shape as actual.
      real(dp), intent(in) :: tolerance !! Mixed absolute/relative tolerance multiplying max(1,abs(expected)).
      character(*), intent(in) :: label !! Name of the reference comparison, printed on failure.
      call require(size(actual) == size(expected), label // ' shape')
      if (.not. all(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected)))) then
         print *, 'actual: ', actual
         print *, 'expected: ', expected
      end if
      call require(all(abs(actual - expected) <= tolerance * max(1.0_dp, abs(expected))), label)
   end subroutine near

   subroutine check_price(yt, yr, ct, cs, mt, expected, recovery, pf, df, accrued, label)
      real(dp), intent(in) :: yt(:) !! Yield times from the frozen R pricing case, in years.
      real(dp), intent(in) :: yr(:) !! R-case continuously compounded annual yield rates.
      real(dp), intent(in) :: ct(:) !! Credit times from the frozen R pricing case, in years.
      real(dp), intent(in) :: cs(:) !! R-case knot survival probabilities.
      real(dp), intent(in) :: mt(:) !! R-case requested CDS maturities in years.
      real(dp), intent(in) :: expected(:) !! Decimal spreads calculated by the authentic credule 0.1.4 package.
      real(dp), intent(in) :: recovery !! Recovery fraction used by the R call.
      integer, intent(in) :: pf !! Premium payments per year used by the R call.
      integer, intent(in) :: df !! Default integration intervals per year used by the R call.
      logical, intent(in) :: accrued !! R-case choice of half-period accrued premiums.
      character(*), intent(in) :: label !! Unique R fixture identifier for diagnostic output.
      real(dp) :: actual(size(mt))
      integer :: status
      call price_cds(yt, yr, ct, cs, mt, recovery, actual, status, pf, df, accrued)
      call require(status == credule_success, label // ' status')
      call near(actual, expected, 2.0e-12_dp, label)
      reference_prices = reference_prices + 1
   end subroutine check_price

   subroutine check_bootstrap(yt, yr, mt, quotes, expected_sp, expected_h, recovery, pf, df, accrued, label)
      real(dp), intent(in) :: yt(:) !! R bootstrap case yield-curve times, in years.
      real(dp), intent(in) :: yr(:) !! R bootstrap case annual zero rates.
      real(dp), intent(in) :: mt(:) !! Strictly increasing CDS maturities supplied to the authentic R package.
      real(dp), intent(in) :: quotes(:) !! Decimal market spreads supplied to the authentic R package.
      real(dp), intent(in) :: expected_sp(:) !! Survival probabilities returned by R credule 0.1.4.
      real(dp), intent(in) :: expected_h(:) !! Interval hazard rates returned by R credule 0.1.4.
      real(dp), intent(in) :: recovery !! R bootstrap case recovery fraction.
      integer, intent(in) :: pf !! R bootstrap case premium payments per year.
      integer, intent(in) :: df !! R bootstrap case default intervals per year.
      logical, intent(in) :: accrued !! R bootstrap case accrual flag.
      character(*), intent(in) :: label !! Unique R fixture identifier for diagnostic output.
      real(dp) :: sp(size(mt)), h(size(mt))
      integer :: status
      call bootstrap_cds(yt, yr, mt, quotes, recovery, sp, h, status, pf, df, accrued)
      call require(status == credule_success, label // ' status')
      call near(sp, expected_sp, 3.0e-12_dp, label // ' survival')
      call near(h, expected_h, 3.0e-12_dp, label // ' hazards')
      reference_bootstraps = reference_bootstraps + 1
   end subroutine check_bootstrap

   subroutine interpolation_tests()
      real(dp) :: yt(2), yr(2), ct(2), cs(2), nan
      yt = [1.0_dp, 2.0_dp]
      yr = [0.02_dp, 0.04_dp]
      ct = [1.0_dp, 3.0_dp]
      cs = exp(-[0.01_dp, 0.11_dp])
      call near([discount_factor(yt, yr, 0.0_dp)], [1.0_dp], 0.0_dp, 'discount time zero')
      call near([discount_factor(yt, yr, 0.5_dp)], [exp(-0.01_dp)], 1.0e-15_dp, 'discount left extrapolation')
      call near([discount_factor(yt, yr, 1.5_dp)], [exp(-0.045_dp)], 1.0e-15_dp, 'linear zero-rate interpolation')
      call near([discount_factor(yt, yr, 3.0_dp)], [exp(-0.12_dp)], 1.0e-15_dp, 'discount right extrapolation')
      call near([survival_probability(ct, cs, 0.0_dp)], [1.0_dp], 0.0_dp, 'survival origin')
      call near([survival_probability(ct, cs, 0.5_dp)], [exp(-0.005_dp)], 1.0e-15_dp, 'survival left hazard')
      call near([survival_probability(ct, cs, 2.0_dp)], [exp(-0.06_dp)], 1.0e-15_dp, 'piecewise hazard interpolation')
      call near([survival_probability(ct, cs, 4.0_dp)], [exp(-0.16_dp)], 1.0e-15_dp, 'survival terminal hazard')
      call near([survival_probability([2.0_dp], [exp(-0.04_dp)], 3.0_dp)], &
                [exp(-0.06_dp)], 1.0e-15_dp, 'single knot extrapolation')
      call require(ieee_is_nan(discount_factor(yt, yr, -1.0_dp)), 'negative discount time rejected')
      call require(ieee_is_nan(survival_probability(ct, cs, -1.0_dp)), 'negative survival time rejected')
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      call require(ieee_is_nan(discount_factor(yt, yr, nan)), 'NaN discount time rejected')
      call require(ieee_is_nan(survival_probability(ct, cs, nan)), 'NaN survival time rejected')
      call require(ieee_is_nan(survival_probability([1.0_dp], [0.0_dp], 1.0_dp)), 'zero knot survival rejected')
      call require(ieee_is_nan(discount_factor([1.0_dp], [0.1_dp, 0.2_dp], 1.0_dp)), 'curve mismatch rejected')
   end subroutine interpolation_tests

   subroutine analytic_tests()
      real(dp), parameter :: rate = 0.02_dp, hazard = 0.035_dp, maturity = 5.0_dp, recovery = 0.4_dp
      real(dp) :: q, total, a, d, expected_a, expected_d, quoted(1), unaccrued(1)
      real(dp) :: premium_sum, accrual
      integer :: pf, df, status
      do pf = 1, 12
         df = 24
         q = exp(-(rate + hazard) / real(pf, dp))
         premium_sum = q * (1.0_dp - q ** int(maturity * real(pf, dp))) / (1.0_dp - q)
         expected_a = premium_sum / real(pf, dp)
         accrual = 0.5_dp * (exp(hazard / real(pf, dp)) - 1.0_dp) * expected_a
         q = exp(-(rate + hazard) / real(df, dp))
         total = q * (1.0_dp - q ** int(maturity * real(df, dp))) / (1.0_dp - q)
         expected_d = (1.0_dp - recovery) * (exp(hazard / real(df, dp)) - 1.0_dp) * total
         call cds_legs([1.0_dp], [rate], [maturity], [exp(-hazard * maturity)], maturity, recovery, &
                       a, d, status, pf, df, .false.)
         call require(status == credule_success, 'analytic unaccrued leg status')
         call near([a, d], [expected_a, expected_d], 4.0e-14_dp, 'geometric-series independent legs')
         call price_cds([1.0_dp], [rate], [maturity], [exp(-hazard * maturity)], [maturity], recovery, &
                        unaccrued, status, pf, df, .false.)
         call near(unaccrued, [expected_d / expected_a], 2.0e-14_dp, 'geometric-series unaccrued spread')
         call cds_legs([1.0_dp], [rate], [maturity], [exp(-hazard * maturity)], maturity, recovery, &
                       a, d, status, pf, df, .true.)
         call near([a, d], [expected_a + accrual, expected_d], 4.0e-14_dp, 'analytic accrued legs')
         call price_cds([1.0_dp], [rate], [maturity], [exp(-hazard * maturity)], [maturity], recovery, &
                        quoted, status, pf, df, .true.)
         call near(quoted, [expected_d / (expected_a + accrual)], 2.0e-14_dp, 'analytic accrued spread')
         call require(quoted(1) < unaccrued(1), 'accrual lowers par spread')
         call require(abs(quoted(1) * a - d) < 1.0e-14_dp, 'premium and protection balance')
      end do
   end subroutine analytic_tests

   subroutine inversion_tests()
      real(dp) :: mt(5), sp(5), hazard(5), quotes(5), recovered(5), fitted(5), repriced(5), low_recovery(5)
      real(dp) :: yt(3), yr(3)
      integer :: i, pf, status
      logical :: accrued
      mt = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 10.0_dp]
      hazard = [0.012_dp, 0.026_dp, 0.011_dp, 0.078_dp, 0.040_dp]
      sp(1) = exp(-hazard(1) * mt(1))
      do i = 2, 5
         sp(i) = sp(i - 1) * exp(-hazard(i) * (mt(i) - mt(i - 1)))
      end do
      yt = [0.5_dp, 3.0_dp, 12.0_dp]
      yr = [-0.005_dp, 0.02_dp, 0.04_dp]
      do pf = 1, 12
         do i = 0, 1
            accrued = i == 1
            call price_cds(yt, yr, mt, sp, mt, 0.4_dp, quotes, status, pf, 12, accrued)
            call require(status == credule_success, 'inversion pricing status')
            call bootstrap_cds(yt, yr, mt, quotes, 0.4_dp, recovered, fitted, status, pf, 12, accrued)
            call require(status == credule_success, 'inversion bootstrap status')
            call near(recovered, sp, 2.0e-12_dp, 'price/bootstrap survival inversion')
            call near(fitted, hazard, 2.0e-12_dp, 'price/bootstrap hazard inversion')
            call price_cds(yt, yr, mt, recovered, mt, 0.4_dp, repriced, status, pf, 12, accrued)
            call near(repriced, quotes, 2.0e-12_dp, 'calibrated market repricing')
            call require(all(fitted >= 0.0_dp), 'nonnegative interval hazards')
            call require(all(recovered(2:5) <= recovered(1:4)), 'nonincreasing calibrated survival')
            call price_cds(yt, yr, mt, sp, mt, 0.0_dp, low_recovery, status, pf, 12, accrued)
            call near(quotes, 0.6_dp * low_recovery, 2.0e-14_dp, 'loss-given-default linearity')
         end do
      end do
   end subroutine inversion_tests

   subroutine edge_tests()
      real(dp) :: s(2), h(2), p(2), one(1), output_three(3), nan, inf, a, d
      real(dp) :: empty(0)
      integer :: status
      nan = ieee_value(0.0_dp, ieee_quiet_nan)
      inf = ieee_value(0.0_dp, ieee_positive_inf)
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], [0.0_dp, 0.0_dp], 0.4_dp, s, h, status)
      call require(status == credule_success, 'zero-spread endpoint bootstrap accepted')
      call near(s, [1.0_dp, 1.0_dp], 0.0_dp, 'zero spreads imply certain survival')
      call near(h, [0.0_dp, 0.0_dp], 0.0_dp, 'zero spreads imply zero hazard')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], s, [1.0_dp, 2.0_dp], 0.4_dp, p, status)
      call near(p, [0.0_dp, 0.0_dp], 0.0_dp, 'default-free CDS par spreads')
      call price_cds([1.0_dp], [0.02_dp], [2.0_dp], [0.8_dp], [1.0_dp, 2.0_dp], 1.0_dp, p, status)
      call require(status == credule_success, 'full recovery pricing accepted')
      call near(p, [0.0_dp, 0.0_dp], 0.0_dp, 'full recovery zero spread')
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], p, 1.0_dp, s, h, status)
      call require(status == credule_invalid_input, 'unidentifiable full recovery rejected')
      call require(all(ieee_is_nan(s)) .and. all(ieee_is_nan(h)), 'invalid bootstrap outputs are NaN')
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], [0.01_dp, 0.0001_dp], 0.4_dp, s, h, status)
      call require(status == credule_no_bracket, 'negative forward hazard rejected')
      call require(all(ieee_is_nan(s)) .and. all(ieee_is_nan(h)), 'failed bootstrap is transactional')
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], [100.0_dp, 100.0_dp], 0.4_dp, s, h, status)
      call require(status == credule_no_bracket, 'hazard above 30 bracket rejected')
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], [0.01_dp, 0.02_dp], &
                         0.4_dp, s, output_three, status)
      call require(status == credule_invalid_input, 'bootstrap mismatched output extents rejected')
      call require(all(ieee_is_nan(output_three)), 'mismatched hazard output filled with NaN')
      call bootstrap_cds([1.0_dp], [0.02_dp], [2.0_dp, 1.0_dp], [0.01_dp, 0.02_dp], 0.4_dp, s, h, status)
      call require(status == credule_invalid_input, 'unordered bootstrap maturities rejected')
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 1.0_dp], [0.01_dp, 0.02_dp], 0.4_dp, s, h, status)
      call require(status == credule_invalid_input, 'duplicate bootstrap maturities rejected')
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], [-0.01_dp, 0.02_dp], 0.4_dp, s, h, status)
      call require(status == credule_invalid_input, 'negative spread rejected')
      call bootstrap_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], [nan, 0.02_dp], 0.4_dp, s, h, status)
      call require(status == credule_invalid_input, 'NaN spread rejected')
      call bootstrap_cds([1.0_dp], [0.02_dp], empty, empty, 0.4_dp, s(:0), h(:0), status)
      call require(status == credule_success, 'empty bootstrap is valid')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], empty, 0.4_dp, p(:0), status)
      call require(status == credule_success, 'empty pricing vector is valid')
      call price_cds(empty, empty, [1.0_dp], [0.9_dp], [1.0_dp], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'empty yield curve rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], [1.0_dp], 0.4_dp, p, status)
      call require(status == credule_invalid_input, 'price output extent rejected')
      call price_cds([1.0_dp], [inf], [1.0_dp], [0.9_dp], [1.0_dp], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'infinite rate rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp, 2.0_dp], [0.8_dp, 0.9_dp], [1.0_dp], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'increasing survival curve rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [1.1_dp], [1.0_dp], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'survival above one rejected')
      call price_cds([1.0_dp], [0.02_dp], [0.0_dp], [0.9_dp], [1.0_dp], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'zero tenor rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], [1.0_dp], -0.1_dp, one, status)
      call require(status == credule_invalid_input, 'negative recovery rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], [nan], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'NaN maturity rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], [1.0_dp], 0.4_dp, one, status, 0)
      call require(status == credule_invalid_input, 'zero premium frequency rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], [1.0_dp], 0.4_dp, one, status, 4, -1)
      call require(status == credule_invalid_input, 'negative default frequency rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], [0.1_dp], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'maturity without premium payment rejected')
      call price_cds([1.0_dp], [0.02_dp], [1.0_dp], [0.9_dp], [huge(1.0_dp)], 0.4_dp, one, status)
      call require(status == credule_invalid_input, 'integer grid overflow rejected')
      call cds_legs([1.0_dp], [0.02_dp], [1.0_dp], [0.0_dp], 1.0_dp, 0.4_dp, a, d, status)
      call require(status == credule_invalid_input, 'invalid legs rejected')
      call require(ieee_is_nan(a) .and. ieee_is_nan(d), 'invalid legs are NaN')
      ! Preserve upstream's frozen tail when pricing beyond the final credit knot.
      call cds_legs([1.0_dp], [0.0_dp], [1.0_dp], [0.9_dp], 2.0_dp, 0.4_dp, a, d, status, 1, 1, .false.)
      call near([a, d], [1.8_dp, 0.06_dp], 1.0e-15_dp, 'upstream zero-hazard pricing tail')
      call cds_legs([1.0_dp], [0.0_dp], [3.0_dp], [exp(-0.15_dp)], &
                    1.24_dp, 0.4_dp, a, d, status, 4, 4, .false.)
      call price_cds([1.0_dp], [0.0_dp], [3.0_dp], [exp(-0.15_dp)], [1.0_dp], 0.4_dp, one, status, 4, 4, .false.)
      call near([d / a], one, 1.0e-15_dp, 'incomplete terminal intervals omitted')
   end subroutine edge_tests
end program test_credule
