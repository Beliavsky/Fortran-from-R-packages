! SPDX-License-Identifier: MIT
program test_fit
   use vasicekfit, only : dp, vasicek_fit_result, fit_vasicek
   implicit none
   real(dp), parameter :: u1(20) = [ &
      -1.5_dp, -1.3421052631578947_dp, -1.1842105263157894_dp, -1.0263157894736843_dp, &
      -0.8684210526315790_dp, -0.7105263157894737_dp, -0.5526315789473685_dp, &
      -0.3947368421052633_dp, -0.23684210526315796_dp, -0.07894736842105265_dp, &
       0.07894736842105265_dp, 0.23684210526315774_dp, 0.39473684210526305_dp, &
       0.5526315789473681_dp, 0.7105263157894735_dp, 0.8684210526315788_dp, &
       1.0263157894736840_dp, 1.1842105263157894_dp, 1.3421052631578947_dp, 1.5_dp ]
   real(dp), parameter :: u2(20) = [ &
      0.20_dp, -0.23_dp, 0.26_dp, -0.29_dp, 0.32_dp, -0.35_dp, 0.38_dp, -0.41_dp, &
      0.44_dp, -0.47_dp, 0.50_dp, -0.53_dp, 0.56_dp, -0.59_dp, 0.62_dp, -0.65_dp, &
      0.68_dp, -0.71_dp, 0.74_dp, -0.77_dp ]
   real(dp), parameter :: y(20) = [ &
      0.009600134840871976_dp, 0.029502890965674730_dp, 0.037198232908477340_dp, &
      0.020785442980663090_dp, 0.050685547743244040_dp, 0.017640200610933460_dp, &
      0.032230717502823025_dp, 0.072262483067218770_dp, 0.016096229432628024_dp, &
      0.048090550856388070_dp, 0.031519799988685336_dp, 0.065329288703111180_dp, &
      0.015700650584264047_dp, 0.060401782241646530_dp, 0.039279530963784316_dp, &
      0.033014580456703310_dp, 0.081703730226032780_dp, 0.045890279005243326_dp, &
      0.071152806570663060_dp, 0.033199582822105970_dp ]
   real(dp) :: x(20,2), y_portfolio(12)
   type(vasicek_fit_result) :: fit, fit_bias, fit_portfolio, bad

   x(:,1) = u1
   x(:,2) = u2
   fit = fit_vasicek(y, x)
   call assert_true(fit%ok)
   call assert_close(fit%beta(1), -1.7985201106695505_dp, 2.0e-12_dp)
   call assert_close(fit%beta(2), 0.1314009664866802_dp, 2.0e-12_dp)
   call assert_close(fit%beta(3), -0.02546891665148888_dp, 2.0e-12_dp)
   call assert_close(fit%sigma2, 0.05009917316270995_dp, 2.0e-12_dp)
   call assert_close(fit%p, 0.03962169865773729_dp, 2.0e-12_dp)
   call assert_close(fit%rho, 0.04770899210578393_dp, 2.0e-12_dp)
   call assert_close(fit%kappa(1), 0.12822815730764678_dp, 2.0e-12_dp)
   call assert_close(fit%kappa(2), -0.02485394391047733_dp, 2.0e-12_dp)
   call assert_close(maxval(abs(fit%fitted + fit%residuals - fit%response_probit)), 0.0_dp, 3.0e-15_dp)

   fit_bias = fit_vasicek(y, x, bias_correct=.true.)
   call assert_true(fit_bias%ok)
   call assert_true(fit_bias%sigma2 > fit%sigma2)
   call assert_close(fit_bias%sigma2, fit%sigma2 * 20.0_dp / 17.0_dp, 2.0e-15_dp)

   y_portfolio = [0.0_dp, 0.01_dp, 0.0_dp, 0.02_dp, 0.03_dp, 0.01_dp, &
      0.04_dp, 0.02_dp, 0.0_dp, 0.01_dp, 0.03_dp, 0.02_dp]
   fit_portfolio = fit_vasicek(y_portfolio, portfolio_size=200)
   call assert_true(fit_portfolio%ok)
   call assert_true(all(fit_portfolio%adjusted_response > 0.0_dp))
   call assert_true(all(fit_portfolio%adjusted_response < 1.0_dp))

   bad = fit_vasicek([0.0_dp, 0.2_dp, 0.3_dp])
   call assert_true(.not. bad%ok)

   print '(a)', 'test_fit: PASS'

contains

   subroutine assert_close(actual, expected, tolerance)
      real(dp), intent(in) :: actual, expected, tolerance
      if (abs(actual - expected) > tolerance + tolerance * abs(expected)) then
         print '(a,3es25.16)', 'mismatch: ', actual, expected, abs(actual - expected)
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_true(condition)
      logical, intent(in) :: condition
      if (.not. condition) error stop 1
   end subroutine assert_true

end program test_fit
