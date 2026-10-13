! SPDX-License-Identifier: GPL-3.0-only
program test_optionstrat
   use optionstrat
   use r_reference_data
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
   implicit none
   type(option_result) :: c, p, rows(3), matrix(2, 2)
   type(dual_option_result) :: both
   type(price_range) :: interval
   type(vertical_result) :: spread
   type(double_vertical_result) :: condor
   real(dp) :: a(6), b(9), cvec(6), pvec(6), s, x, v, t, r, d, h, price, expected, q, inf, trials(5)
   real(dp) :: reference(10), observed(10), values(2, 2), empty(0)
   character(len=7), parameter :: greek(6) = [character(len=7) :: 'premium', 'delta', 'gamma', 'vega', 'theta', 'rho']
   character(len=4) :: family
   integer :: i, j, checks
   checks = 0
   do i = 1, size(analytics_ref, 2)
      a = analytics_ref(1:6, i)
      c = calleval(a(1), a(2), a(3), a(4), a(5), a(6))
      p = puteval(a(1), a(2), a(3), a(4), a(5), a(6))
      both = opteval(a(1), a(2), a(3), a(4), a(5), a(6))
      cvec = unpack(c)
      pvec = unpack(p)
      do j = 1, 6
         call near(cvec(j), analytics_ref(j + 6, i), 3.0e-12_dp, 'genuine R call')
         call near(pvec(j), analytics_ref(j + 12, i), 3.0e-12_dp, 'genuine R put')
         call near(callgreek(trim(greek(j)), a(1), a(2), a(3), a(4), a(5), a(6)), &
            cvec(j), 0.0_dp, 'call selector')
         call near(putgreek(trim(greek(j)), a(1), a(2), a(3), a(4), a(5), a(6)), &
            pvec(j), 0.0_dp, 'put selector')
      end do
      call near(both%call%premium, c%premium, 0.0_dp, 'dual call')
      call near(both%put%rho, p%rho, 0.0_dp, 'dual put')
      call near(c%premium - p%premium, a(1) * exp(-a(6) * a(4)) - a(2) * exp(-a(5) * a(4)), &
         2.0e-12_dp, 'put-call parity')
      call truth(c%delta >= 0.0_dp .and. p%delta <= 0.0_dp, 'signed deltas')
      call truth(c%rho >= 0.0_dp .and. p%rho <= 0.0_dp, 'signed rhos')
      if (i <= 6) then
         call near(lambda('call', a(1), a(2), a(3), a(4), a(5), a(6)), extra_ref(1, i), &
            2.0e-11_dp, 'R call elasticity')
         call near(lambda('put', a(1), a(2), a(3), a(4), a(5), a(6)), extra_ref(2, i), &
            2.0e-11_dp, 'R put elasticity')
         call near(iv_calc('call', c%premium, a(1), a(2), a(4), a(5), a(6)), extra_ref(3, i), &
            1.0e-13_dp, 'R call grid IV')
         call near(iv_calc('put', p%premium, a(1), a(2), a(4), a(5), a(6)), extra_ref(4, i), &
            1.0e-13_dp, 'R put grid IV')
         call near(implied_volatility('call', c%premium, a(1), a(2), a(4), a(5), a(6)), a(3), &
            5.0e-11_dp, 'continuous call IV')
         call near(implied_volatility('put', p%premium, a(1), a(2), a(4), a(5), a(6)), a(3), &
            5.0e-11_dp, 'continuous put IV')
      end if
   end do
   do i = 1, size(probabilities_ref, 2)
      b = probabilities_ref(1:9, i)
      call near(prob_below(b(1), b(3), b(4), b(5), b(6), b(7), b(8)), probabilities_ref(10, i), &
         2.0e-13_dp, 'R probability below')
      call near(prob_above(b(1), b(2), b(4), b(5), b(6), b(7), b(8)), probabilities_ref(11, i), &
         2.0e-13_dp, 'R probability above')
      call near(prob_btwn(b(1), b(2), b(3), b(4), b(5), b(6), b(7), b(8)), probabilities_ref(12, i), &
         2.0e-13_dp, 'R probability between')
      call near(prob_below_quantile(b(1), b(9), b(4), b(5), b(6), b(7), b(8)), probabilities_ref(13, i), &
         2.0e-13_dp, 'R below quantile')
      call near(prob_above_quantile(b(1), b(9), b(4), b(5), b(6), b(7), b(8)), probabilities_ref(14, i), &
         2.0e-13_dp, 'R above quantile convention')
      interval = prob_btwn_quantile(b(1), b(9), b(4), b(5), b(6), b(7), b(8))
      call near(interval%lower_price, probabilities_ref(15, i), 2.0e-13_dp, 'R between lower quantile')
      call near(interval%upper_price, probabilities_ref(16, i), 2.0e-13_dp, 'R between upper quantile')
   end do
   do i = 1, 4
      family = 'call'
      if (i > 2) family = 'put'
      x = 90.0_dp
      s = 110.0_dp
      if (mod(i, 2) == 0) then
         x = 110.0_dp
         s = 90.0_dp
      end if
      spread = vertical(trim(family), 100.0_dp, x, s, 0.25_dp, 0.02_dp, 0.2_dp, 0.25_dp, 0.3_dp, 0.01_dp)
      observed = [spread%spot, spread%short_strike, spread%long_strike, spread%max_profit, spread%max_loss, &
         spread%breakeven, spread%prob_be, spread%prob_max_profit, spread%prob_max_loss, spread%initial_dc]
      reference = spreads_ref(:, i)
      call truth(all(abs(observed - reference) <= 0.005001_dp), 'R rounded vertical analytics')
      call near(spread%max_profit - spread%max_loss, abs(s - x), 2.0e-14_dp, 'vertical payoff width')
   end do
   condor = dv(100.0_dp, 80.0_dp, 90.0_dp, 110.0_dp, 120.0_dp, 0.25_dp, 0.02_dp, &
      0.2_dp, 0.21_dp, 0.22_dp, 0.23_dp, 0.3_dp, 0.01_dp)
   observed = [condor%spot, condor%long_put, condor%short_put, condor%short_call, condor%long_call, &
      condor%lower_be, condor%higher_be, condor%prob_be, condor%prob_max_profit, condor%prob_max_loss]
   call truth(all(abs(observed - spreads_ref(:, 5)) <= 0.005001_dp), 'R rounded double vertical analytics')

   ! Independent Black-Scholes value: N is evaluated directly with ERFC.
   s = 100.0_dp
   x = 100.0_dp
   v = 0.2_dp
   t = 0.25_dp
   r = 0.03_dp
   d = 0.01_dp
   expected = s * exp(-d * t) * 0.5_dp * erfc(-0.1_dp / sqrt(2.0_dp)) &
      - x * exp(-r * t) * 0.5_dp
   call near(callpremium(s, x, v, t, r, d), expected, 2.0e-13_dp, 'independent BSM price')
   h = 1.0e-4_dp
   call near(calldelta(s, x, v, t, r, d), (callpremium(s + h, x, v, t, r, d) &
      - callpremium(s - h, x, v, t, r, d)) / (2.0_dp * h), 2.0e-8_dp, 'call delta finite difference')
   call near(putdelta(s, x, v, t, r, d), (putpremium(s + h, x, v, t, r, d) &
      - putpremium(s - h, x, v, t, r, d)) / (2.0_dp * h), 2.0e-8_dp, 'put delta finite difference')
   call near(optiongamma(s, x, v, t, r, d), (calldelta(s + h, x, v, t, r, d) &
      - calldelta(s - h, x, v, t, r, d)) / (2.0_dp * h), 2.0e-9_dp, 'gamma finite difference')
   h = 1.0e-6_dp
   call near(optionvega(s, x, v, t, r, d), 0.01_dp * (callpremium(s, x, v + h, t, r, d) &
      - callpremium(s, x, v - h, t, r, d)) / (2.0_dp * h), 2.0e-9_dp, 'vega per one percent')
   call near(callrho(s, x, v, t, r, d), 0.01_dp * (callpremium(s, x, v, t, r + h, d) &
      - callpremium(s, x, v, t, r - h, d)) / (2.0_dp * h), 2.0e-9_dp, 'call rho per one percent')
   call near(putrho(s, x, v, t, r, d), 0.01_dp * (putpremium(s, x, v, t, r + h, d) &
      - putpremium(s, x, v, t, r - h, d)) / (2.0_dp * h), 2.0e-9_dp, 'put rho per one percent')
   call near(calltheta_standard(s, x, v, t, r, d), -(callpremium(s, x, v, t + h, r, d) &
      - callpremium(s, x, v, t - h, r, d)) / (730.0_dp * h), 2.0e-9_dp, 'standard call theta derivative')
   call near(puttheta_standard(s, x, v, t, r, d), -(putpremium(s, x, v, t + h, r, d) &
      - putpremium(s, x, v, t - h, r, d)) / (730.0_dp * h), 2.0e-9_dp, 'standard put theta derivative')
   call truth(abs(calltheta(s, x, v, t, r, d) - calltheta_standard(s, x, v, t, r, d)) > 0.001_dp, &
      'upstream theta defect remains explicit')

   call near(calltheta_standard(1.0e200_dp, 1.0e200_dp, v, t, r, d) / 1.0e198_dp, &
      calltheta_standard(100.0_dp, 100.0_dp, v, t, r, d), 2.0e-14_dp, 'large-scale standard call theta')
   call near(puttheta_standard(1.0e200_dp, 1.0e200_dp, v, t, r, d) / 1.0e198_dp, &
      puttheta_standard(100.0_dp, 100.0_dp, v, t, r, d), 2.0e-14_dp, 'large-scale standard put theta')
   trials = [0.0001_dp, 0.0357_dp, 0.2_dp, 1.23456_dp, 6.0_dp]
   do i = 1, size(trials)
      price = callpremium(s, x, trials(i), 1.0_dp, 0.0_dp)
      call near(implied_volatility('call', price, s, x, 1.0_dp, 0.0_dp), trials(i), 2.0e-11_dp, 'ATM IV roundtrip')
   end do
   price = callpremium(s, x, 0.2_dp, 1.0_dp, 0.0_dp)
   call near(iv_calc('call', price, s, x, 1.0_dp, 0.0_dp), 0.0_dp, 0.0_dp, 'R ATM IV grid defect')
   call near(implied_volatility('call', 0.0_dp, s, x, 1.0_dp, 0.0_dp), 0.0_dp, 0.0_dp, 'IV lower bound')
   call truth(ieee_is_nan(implied_volatility('call', s, s, x, 1.0_dp, 0.0_dp)), 'IV upper bound has no finite root')
   call truth(ieee_is_nan(implied_volatility('put', -1.0_dp, s, x, 1.0_dp, 0.0_dp)), 'IV below price bound')
   call truth(ieee_is_nan(implied_volatility('call', 1.0_dp, s, x, 0.0_dp, 0.0_dp)), 'IV expired undefined')
   call near(callpremium(110.0_dp, 100.0_dp, 0.0_dp, 1.0_dp, 0.0_dp), 10.0_dp, 0.0_dp, 'zero volatility call')
   call near(putpremium(90.0_dp, 100.0_dp, 0.0_dp, 1.0_dp, 0.0_dp), 10.0_dp, 0.0_dp, 'zero volatility put')
   call near(calldelta(100.0_dp, 100.0_dp, 0.0_dp, 1.0_dp, 0.0_dp), 0.5_dp, 0.0_dp, 'zero-vol kink delta')
   call truth(ieee_is_nan(optiongamma(100.0_dp, 100.0_dp, 0.0_dp, 1.0_dp, 0.0_dp)), 'zero-vol kink gamma')
   call near(callpremium(100.0_dp, 90.0_dp, 0.2_dp, 0.0_dp, 0.02_dp), 10.0_dp, 0.0_dp, 'zero tenor price')
   call truth(optiongamma(100.0_dp, 100.0_dp, 0.2_dp, 0.0_dp, 0.02_dp) > 1.0e90_dp, 'tiny tenor gamma')
   call near(callpremium(1.0e-100_dp, 1.0e100_dp, 0.2_dp, 1.0_dp, 0.0_dp), 0.0_dp, 0.0_dp, 'extreme OTM')
   call near(putpremium(1.0e-100_dp, 1.0e100_dp, 0.2_dp, 1.0_dp, 0.0_dp), 1.0e100_dp, 1.0e-14_dp, 'extreme ITM')
   rows = calleval([90.0_dp, 100.0_dp, 110.0_dp], 100.0_dp, 0.2_dp, 1.0_dp, 0.02_dp)
   call truth(size(rows) == 3 .and. all(rows%premium >= 0.0_dp), 'elemental rank-one results')
   values = reshape([80.0_dp, 90.0_dp, 100.0_dp, 110.0_dp], [2, 2])
   matrix = puteval(values, 100.0_dp, 0.2_dp, 1.0_dp, 0.02_dp)
   call truth(all(shape(matrix) == [2, 2]), 'elemental rank-two shape')
   call truth(size(callpremium(empty, 100.0_dp, 0.2_dp, 1.0_dp, 0.02_dp)) == 0, 'zero-length shape')
   price = prob_above_quantile_standard(100.0_dp, 0.2_dp, asd = 0.2_dp, dte = 30.0_dp, mean = 0.04_dp)
   call near(prob_above(100.0_dp, price, asd = 0.2_dp, dte = 30.0_dp, mean = 0.04_dp), &
      0.2_dp, 2.0e-12_dp, 'corrected above quantile roundtrip')
   interval = prob_btwn_quantile_standard(100.0_dp, 0.8_dp, asd = 0.2_dp, dte = 30.0_dp, mean = 0.04_dp)
   call near(prob_btwn(100.0_dp, interval%lower_price, interval%upper_price, asd = 0.2_dp, &
      dte = 30.0_dp, mean = 0.04_dp), 0.8_dp, 2.0e-12_dp, 'corrected interval quantile roundtrip')
   call near(prob_below(100.0_dp, 100.0_dp, dte = 10.0_dp), 1.0_dp, 0.0_dp, 'degenerate normal CDF')
   call near(prob_above(100.0_dp, 100.0_dp, dte = 10.0_dp), 0.0_dp, 0.0_dp, 'degenerate upper tail')
   call truth(.not. ieee_is_finite(prob_below_quantile(100.0_dp, 0.0_dp, dsd = 0.01_dp, dte = 1.0_dp)), &
      'quantile endpoint infinity')
   call near(tdiff('2020-02-28', '2020-03-01', 'days'), 2.0_dp, 0.0_dp, 'leap-year date difference')
   call near(tdiff('2019-02-28', '2019-03-01', 'days'), 1.0_dp, 0.0_dp, 'ordinary date difference')
   call near(tdiff('2000-03-01', '2000-02-28', 'years'), -2.0_dp / 365.0_dp, 0.0_dp, 'signed year difference')
   call near(r_cont(0.12_dp, 2.0_dp), 2.0_dp * log(1.06_dp), 1.0e-15_dp, 'nominal rate conversion')
   call near(r_cont(1.0e-15_dp, 12.0_dp), 1.0e-15_dp, 1.0e-28_dp, 'stable tiny nominal rate')

   q = ieee_value(0.0_dp, ieee_quiet_nan)
   inf = ieee_value(0.0_dp, ieee_positive_inf)
   call truth(ieee_is_nan(callpremium(-1.0_dp, x, v, t, r)), 'negative spot invalid')
   call truth(ieee_is_nan(callpremium(s, 0.0_dp, v, t, r)), 'zero strike invalid')
   call truth(ieee_is_nan(putpremium(s, x, -v, t, r)), 'negative volatility invalid')
   call truth(ieee_is_nan(callpremium(s, x, v, -t, r)), 'negative tenor invalid')
   call truth(ieee_is_nan(callpremium(s, x, v, t, q)), 'NaN invalid')
   call truth(ieee_is_nan(callpremium(s, x, v, t, inf)), 'infinity invalid')
   call truth(ieee_is_nan(callpremium(s, x, v, t, -1.0e5_dp)), 'discount overflow invalid')
   call truth(ieee_is_nan(callgreek('unknown', s, x, v, t, r)), 'bad Greek invalid')
   call truth(ieee_is_nan(lambda('stock', s, x, v, t, r)), 'bad type invalid')
   call truth(ieee_is_nan(lambda('call', 1.0_dp, 100.0_dp, 0.0_dp, t, r)), 'zero-price elasticity undefined')
   call truth(ieee_is_nan(iv_calc('bad', 1.0_dp, s, x, t, r)), 'grid bad type invalid')
   call truth(ieee_is_nan(prob_below(s, x)), 'all-zero probability inputs invalid')
   call truth(ieee_is_nan(prob_above(s, x, asd = -1.0_dp)), 'negative probability scale invalid')
   call truth(ieee_is_nan(prob_below_quantile(s, 1.1_dp, dsd = 0.01_dp, dte = 10.0_dp)), 'bad quantile invalid')
   call truth(ieee_is_nan(prob_below(s, x, asd = 0.2_dp, tradedays = 0.0_dp)), 'bad trading days invalid')
   call truth(ieee_is_nan(tdiff('1900-02-29', '1900-03-01', 'days')), 'invalid leap date')
   call truth(ieee_is_nan(tdiff('2000-02-29', '2000-03-01', 'months')), 'invalid period')
   call truth(ieee_is_nan(r_cont(-2.0_dp, 1.0_dp)), 'invalid nominal rate')
   call truth(ieee_is_nan(r_cont(0.1_dp, 0.0_dp)), 'invalid compounding')
   spread = vertical('call', s, x, x, t, r, v)
   call truth(ieee_is_nan(spread%max_profit), 'equal-strike spread invalid')
   print '(a,i0,a)', 'optionstrat: ', checks, ' deterministic checks passed'
contains
   pure function unpack(value) result(vector)
      type(option_result), intent(in) :: value !! Scalar option analytics to place in premium-to-rho column order.
      real(dp) :: vector(6)
      vector = [value%premium, value%delta, value%gamma, value%vega, value%theta, value%rho]
   end function unpack
   subroutine truth(condition, label)
      logical, intent(in) :: condition !! Assertion predicate; false terminates the deterministic test.
      character(len=*), intent(in) :: label !! Descriptive check name printed if the assertion fails.
      checks = checks + 1
      if (.not. condition) then
         print '(a)', 'FAILED: ' // label
         error stop 1
      end if
   end subroutine truth
   subroutine near(actual, wanted, tolerance, label)
      real(dp), intent(in) :: actual !! Computed scalar to compare, required to be finite.
      real(dp), intent(in) :: wanted !! Finite trusted reference scalar.
      real(dp), intent(in) :: tolerance !! Relative tolerance with an absolute floor of the same size.
      character(len=*), intent(in) :: label !! Descriptive comparison label for failure reporting.
      if (.not. ieee_is_finite(actual) .or. abs(actual - wanted) > tolerance * max(1.0_dp, abs(wanted))) then
         print '(a,2es25.16)', trim(label) // ': ', actual, wanted
         error stop 1
      end if
      checks = checks + 1
   end subroutine near
end program test_optionstrat
