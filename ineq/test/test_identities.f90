! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_identities
   use ineq_api
   use test_support
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_nan, ieee_is_finite
   implicit none
   real(dp), parameter :: x(5) = [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 10.0_dp]
   real(dp), parameter :: scales(5) = [1.0e-300_dp, 1.0e-150_dp, 1.0_dp, 1.0e150_dp, 1.0e300_dp]
   real(dp) :: pair_sum, g, scale, area, p, previous(5), current(5), gap(5), nan
   real(dp) :: cv, variance, expected
   real(dp), allocatable :: values(:)
   type(lorenz_curve) :: curve, base, bound
   logical :: dominates
   integer :: i, j, status
   nan = ieee_value(0.0_dp, ieee_quiet_nan)
   pair_sum = 0.0_dp
   do i = 1, size(x)
      do j = 1, size(x)
         pair_sum = pair_sum + abs(x(i) - x(j))
      end do
   end do
   g = pair_sum / (2.0_dp * real(size(x), dp) * sum(x))
   call check(gini(x), g, 'pairwise Gini identity')
   call check(gini(x, corr=.true.), g * 5.0_dp / 4.0_dp, 'finite-sample Gini correction')
   base = lc(x)
   area = sum((base%p(2:) - base%p(:5)) * (base%l(2:) + base%l(:5))) / 2.0_dp
   call check(g, 1.0_dp - 2.0_dp * area, 'Gini from Lorenz trapezoids')
   call check(rs(x), maxval(base%p - base%l), 'RS from maximum Lorenz gap')
   cv = var_coeff(x)
   call check(herfindahl(x), (1.0_dp + cv ** 2) / 5.0_dp, 'Herfindahl variance identity')
   call check(rosenbluth(x), 1.0_dp / (5.0_dp * (1.0_dp - g)), 'Rosenbluth Gini identity')
   do i = 1, size(scales)
      scale = scales(i)
      call check(gini(scale * x), g, 'Gini scale invariance')
      call check(rs(scale * x), rs(x), 'RS scale invariance')
      call check(var_coeff(scale * x), cv, 'CV scale invariance without variance overflow')
      call check(atkinson(scale * x), atkinson(x), 'Atkinson scale invariance')
      call check(atkinson(scale * x, 1.0_dp), atkinson(x, 1.0_dp), 'geometric Atkinson scale invariance')
      call check(theil(scale * x), theil(x), 'Theil T scale invariance')
      call check(theil(scale * x, 1.0_dp), theil(x, 1.0_dp), 'Theil L scale invariance')
      call check(entropy(scale * x), entropy(x), 'entropy scale invariance')
      call check(herfindahl(scale * x), herfindahl(x), 'Herfindahl scale invariance')
      call check(rosenbluth(scale * x), rosenbluth(x), 'Rosenbluth scale invariance')
      curve = lc(scale * x)
      call require(curve%status == ineq_ok, 'scaled Lorenz status')
      call check_array(curve%p, base%p, 'scaled Lorenz population')
      call check_array(curve%l, base%l, 'scaled Lorenz income')
      call check_array(curve%l_general / scale, base%l_general, 'generalized Lorenz scales with income')
      call check(watts(scale * x, 5.0_dp * scale), watts(x, 5.0_dp), 'Watts joint scale invariance')
      call check(sen(scale * x, 5.0_dp * scale), sen(x, 5.0_dp), 'Sen joint scale invariance')
      call check(sst(scale * x, 5.0_dp * scale), sst(x, 5.0_dp), 'SST joint scale invariance')
      call check(foster(scale * x, 5.0_dp * scale, 3.0_dp), foster(x, 5.0_dp, 3.0_dp), &
         'Foster joint scale invariance')
   end do
   curve = lc(1.0e200_dp * x, [1.0e200_dp, 1.0e200_dp, 1.0e200_dp, 1.0e200_dp, 1.0e200_dp])
   call require(curve%status == ineq_ok, 'Lorenz avoids income-times-weight overflow')
   call check_array(curve%l, base%l, 'Lorenz extreme income times frequency')
   call check_array(curve%l_general / 1.0e200_dp, base%l_general, 'generalized extreme Lorenz')
   curve = lc(1.0e-200_dp * x, [1.0e-200_dp, 1.0e-200_dp, 1.0e-200_dp, 1.0e-200_dp, 1.0e-200_dp])
   call require(curve%status == ineq_ok, 'Lorenz avoids income-times-weight underflow')
   call check_array(curve%l, base%l, 'Lorenz tiny income times frequency')
   gap = max(0.0_dp, 1.0_dp - x / 5.0_dp)
   call check(foster(x, 5.0_dp), 4.0_dp / 5.0_dp, 'Foster default is the headcount ratio')
   call check(foster(x, 5.0_dp, 2.0_dp), sum(gap) / 5.0_dp, 'Foster parameter two is mean poverty gap')
   call check(foster(x, 5.0_dp, 3.0_dp), sum(gap ** 2) / 5.0_dp, 'Foster parameter three is squared gap')
   call check(sst(x, 5.0_dp), sum(gap) / 5.0_dp * (1.0_dp + gini(gap)), 'SST gap-Gini identity')
   call check(foster([1.0_dp, 5.0_dp, 5.0_dp], 5.0_dp), 1.0_dp / 3.0_dp, 'poverty line is strict')
   call check(pov(x, 0.0_dp), 0.0_dp, 'no poor gives zero')
   call check(watts([real(dp) :: ], 1.0_dp), 0.0_dp, 'empty poverty sample gives zero')
   call check(foster([nan], 1.0_dp), 0.0_dp, 'all-missing poverty sample gives zero')
   call require(ieee_is_nan(sen([0.0_dp, 0.0_dp], 1.0_dp)), 'upstream all-zero-poor Sen remains NaN')
   call check(theil([0.0_dp, 1.0_dp, 1.0_dp]), 0.0_dp, 'Theil T removes zeros before its mean')
   call check(entropy([0.0_dp, 1.0_dp, 1.0_dp], 1.0_dp), 0.0_dp, 'entropy order one follows Theil zero removal')
   call check(atkinson(x, 0.0_dp), 0.0_dp, 'Atkinson zero aversion')
   expected = 1.0_dp - exp(sum(log(x)) / 5.0_dp) / (sum(x) / 5.0_dp)
   call check(atkinson(x, 1.0_dp - 1.0e-10_dp), expected, 'Atkinson continuous geometric-mean limit', 1.0e-10_dp)
   call check(atkinson(x, 1.0_dp + 1.0e-10_dp), expected, 'Atkinson limit from above', 1.0e-10_dp)
   call check(entropy(x, 1.0e-10_dp), theil(x, 1.0_dp), 'entropy order-zero limit', 1.0e-10_dp)
   call check(entropy(x, 1.0_dp - 1.0e-10_dp), theil(x), 'entropy order-one limit', 1.0e-10_dp)
   call check(herfindahl(x, 1.0e-10_dp), exp(sum(x / sum(x) * log(x / sum(x)))), &
      'Herfindahl nonzero-order limit', 1.0e-10_dp)
   call check(herfindahl(x, 0.0_dp), 1.0_dp, 'Herfindahl exact-zero R convention')
   call check(herfindahl([real(dp) :: ]), 0.0_dp, 'Herfindahl empty convention')
   call check(rosenbluth([real(dp) :: ]), -1.0_dp, 'Rosenbluth empty convention')
   variance = sum((x - sum(x) / 5.0_dp) ** 2) / 5.0_dp
   call check(kolm(x, 1.0e-8_dp) / 1.0e-8_dp, variance / 2.0_dp, 'Kolm small-aversion variance limit', 1.0e-7_dp)
   call check(kolm(x + 1000.0_dp), kolm(x), 'Kolm translation invariance')
   call check(kolm(-x, -2.0_dp), -kolm(x, 2.0_dp), 'Kolm reflection and parameter sign')
   call check(kolm(1.0e300_dp * x) / 1.0e300_dp, sum(x) / 5.0_dp - minval(x), 'Kolm avoids exponential overflow')
   call require(ieee_is_finite(kolm(x, 1.0e300_dp)), 'Kolm huge aversion remains finite')
   call check(theil([1.0e-300_dp, 1.0e300_dp]), log(2.0_dp), 'Theil T heterogeneous 600-order incomes')
   call check(theil([1.0e-300_dp, 1.0e300_dp], 1.0_dp), log(5.0e299_dp), &
      'mean log deviation avoids normalized-income underflow')
   expected = 0.5_dp * 2.0_dp ** (-0.1_dp) * 1.0e60_dp / 0.11_dp
   call check(entropy([1.0e-300_dp, 1.0e300_dp], -0.1_dp) / expected, 1.0_dp, &
      'negative-order entropy avoids normalized-income underflow')
   expected = exp(-600.0_dp * log(10.0_dp) / 11.0_dp)
   call check(herfindahl([1.0e-300_dp, 1.0e300_dp], -1.1_dp) / expected, 1.0_dp, &
      'negative-order concentration avoids normalized-income underflow')
   call check(watts([1.0e-300_dp, 1.0e300_dp], 1.0e300_dp), 300.0_dp * log(10.0_dp), &
      'Watts avoids poverty-ratio overflow')
   call check(atkinson([1.0e-300_dp, 1.0e300_dp], 1.0_dp), 1.0_dp, 'geometric Atkinson extreme spread')
   call require(ieee_is_nan(atkinson(x, nan)), 'NaN index parameter is rejected')
   call require(ieee_is_nan(lc_lognorm(0.5_dp, nan)), 'NaN lognormal parameter is rejected')
   values = lasym([10.0_dp, 51.0_dp, 44.0_dp, 10.0_dp, 87.0_dp, 104.0_dp], interval=.true.)
   call check_array(values, [0.70915032679738554_dp, 1.0424836601307189_dp], 'Lasym exact-mean tie regression')
   values = lasym([0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp], [1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp], interval=.true.)
   call check_array(values, [0.4_dp, 1.0_dp], 'weighted Lasym exact mean tie')
   values = lasym([1.0_dp, nan, 3.0_dp], interval=.true., na_rm=.false.)
   call require(size(values) == 2, 'Lasym keepNA interval length')
   call require(all(ieee_is_nan(values)), 'Lasym keepNA interval values')
   call major([106.0_dp, 228.0_dp, 338.0_dp, 720.0_dp], &
      [119.0_dp, 228.0_dp, 338.0_dp, 707.0_dp], dominates, status)
   call require(dominates .and. status == ineq_ok, 'majorization equal-integer-total regression')
   call major([1.0_dp], [2.0_dp], dominates, status)
   call require(.not. dominates .and. status == ineq_incomparable, 'majorization incomparable total')
   call major([1.0_dp, 2.0_dp], [3.0_dp], dominates, status)
   call require(status == ineq_incomparable, 'majorization incomparable length')
   call major([real(dp) :: ], [real(dp) :: ], dominates, status)
   call require(dominates .and. status == ineq_ok, 'empty vectors majorize each other')
   previous = 0.0_dp
   do i = 0, 100
      p = real(i, dp) / 100.0_dp
      current = [lc_exp(p), lc_lognorm(p), lc_pareto(p), lc_singh(p), lc_dagum(p)]
      call require(all(current >= previous - 1.0e-14_dp), 'theoretical Lorenz monotonicity')
      call require(all(current >= -1.0e-14_dp) .and. all(current <= p + 1.0e-14_dp), 'theoretical Lorenz bounds')
      previous = current
   end do
   call check_array(previous, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp], 'Lorenz endpoint extensions')
   call check(lc_exp(1.0e-100_dp) / 1.0e-200_dp, 0.5_dp, 'exponential Lorenz tiny-population limit')
   call check(lc_exp(1.0e-10_dp) / 1.0e-20_dp, 0.5_dp, 'exponential Lorenz avoids cancellation', 1.0e-10_dp)
   call check(lc_pareto(0.75_dp, 2.0_dp), 0.5_dp, 'Pareto analytic identity')
   call check(lc_dagum(0.6_dp), 1.0_dp - sqrt(1.0_dp - 0.6_dp ** 2), 'Dagum default beta closed form')
   call check(lc_singh(0.5_dp, [1.5_dp, 1.0_dp]), 0.0_dp, 'Singh zero beta shape limit')
   call check(lc_singh(1.0_dp, [1.5_dp, 1.0_dp]), 1.0_dp, 'Singh zero beta shape endpoint')
   call check(lc_dagum(0.0_dp, [2.0_dp, -2.0_dp]), 1.0_dp, 'upstream degenerate beta parameter')
   call require(ieee_is_nan(lc_exp(-0.1_dp)), 'invalid theoretical population probability')
   call require(ieee_is_nan(theor_lc(0.5_dp, 'unknown')), 'invalid theoretical dispatcher')
   call require(ieee_is_nan(ineq(x, type='unknown')), 'invalid inequality dispatcher')
   call require(ieee_is_nan(conc(x, type='unknown')), 'invalid concentration dispatcher')
   call require(ieee_is_nan(pov(x, 1.0_dp, type='unknown')), 'invalid poverty dispatcher')
   curve = lc([0.0_dp, 0.0_dp])
   call require(curve%status == ineq_invalid .and. size(curve%p) == 0, 'undefined empirical curve status')
   curve = lc(x, [1.0_dp])
   call require(curve%status == ineq_invalid, 'mismatched frequencies status')
   bound = lc_mehran([1.0_dp], [1.0_dp])
   call require(bound%status == ineq_invalid, 'Mehran requires at least two classes')
   if (failed > 0) error stop 1
   print *, 'Independent identity and edge checks passed:', checked
end program test_identities
