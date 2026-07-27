! SPDX-License-Identifier: GPL-2.0-or-later
program test_fitting
   use rnd, only : dp, option_prices, bsm_fit, shimko_fit, rate_result
   use rnd, only : price_bsm_option, extract_bsm_density, compute_implied_volatility
   use rnd, only : extract_shimko_density, extract_rates, get_point_estimate
   implicit none
   real(dp), parameter :: s0=100.0_dp, r=0.03_dp, y=0.01_dp, te=0.75_dp, sigma=0.22_dp
   real(dp) :: strikes(9), initial(2), true_mu, true_zeta
   real(dp), allocatable :: implied(:), point(:)
   logical :: status(9)
   type(option_prices) :: market
   type(bsm_fit) :: fit
   type(shimko_fit) :: shimko
   type(rate_result) :: rates
   integer :: i

   do i = 1, size(strikes)
      strikes(i) = 80.0_dp+5.0_dp*real(i-1,dp)
   end do
   market = price_bsm_option(s0,strikes,r,te,sigma,y)
   true_mu = log(s0)+(r-y-0.5_dp*sigma*sigma)*te
   true_zeta = sigma*sqrt(te)
   initial = [true_mu+0.02_dp,true_zeta*1.1_dp]
   fit = extract_bsm_density(r,y,te,s0,market%call,strikes,market%put,strikes, &
      initial_values=initial,max_iter=3000)
   call assert_close(fit%mu,true_mu,2.0e-5_dp,"fitted mu")
   call assert_close(fit%zeta,true_zeta,2.0e-5_dp,"fitted zeta")

   implied = compute_implied_volatility(r,te,s0,strikes,y,market%call,1.0e-6_dp,2.0_dp,status)
   call assert_true(all(status),"implied volatility convergence")
   call assert_true(maxval(abs(implied-sigma)) < 2.0e-10_dp,"implied volatility recovery")

   shimko = extract_shimko_density(market%call,strikes,r,y,te,s0)
   call assert_close(shimko%a0,sigma,2.0e-9_dp,"Shimko intercept")
   call assert_true(abs(shimko%a1) < 2.0e-10_dp .and. abs(shimko%a2) < 2.0e-12_dp, &
      "Shimko constant curve")

   rates = extract_rates(market%call,market%put,s0,strikes,te)
   call assert_close(rates%risk_free_rate,r,2.0e-13_dp,"rate extraction")
   call assert_close(rates%dividend_yield,y,2.0e-13_dp,"yield extraction")
   point = get_point_estimate(market%call,strikes,r,te)
   call assert_true(size(point) == size(strikes)-2,"point estimate size")
   call assert_true(all(point > 0.0_dp),"point estimate positivity")
   print '(a)', 'test_fitting: PASS'
contains
   subroutine assert_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual, expected, tolerance
      character(*), intent(in) :: message
      if (abs(actual-expected) > tolerance) then
         print '(a,2es20.10)', trim(message)//': ',actual,expected
         error stop 1
      end if
   end subroutine assert_close
   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      if (.not. condition) then
         print '(a)', trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_fitting
