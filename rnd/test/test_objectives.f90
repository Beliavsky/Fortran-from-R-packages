! SPDX-License-Identifier: GPL-2.0-or-later
program test_objectives
   use rnd, only : dp, option_prices
   use rnd, only : price_bsm_option, price_gb_option, price_mln_option, price_ew_option
   use rnd, only : bsm_objective, gb_objective, mln_objective, ew_objective, mln_am_objective
   implicit none
   real(dp), parameter :: s0=100.0_dp, r=0.03_dp, y=0.01_dp, te=0.5_dp
   real(dp) :: strikes(5), weights(5), theta_bsm(2), theta_gb(4), theta_mln(5), theta_ew(3)
   real(dp) :: theta_am(10), sigma, v, skew, kurt, value
   type(option_prices) :: prices
   integer :: i

   do i = 1, size(strikes)
      strikes(i) = 90.0_dp+5.0_dp*real(i-1,dp)
   end do
   weights = 1.0_dp
   sigma = 0.2_dp
   prices = price_bsm_option(s0,strikes,r,te,sigma,y)
   theta_bsm = [log(s0)+(r-y-0.5_dp*sigma*sigma)*te,sigma*sqrt(te)]
   value = bsm_objective(theta_bsm,s0,r,te,y,prices%call,strikes,weights, &
      prices%put,strikes,weights,1.0_dp)
   call assert_true(value < 1.0e-20_dp,"BSM objective at truth")

   theta_gb = [2.0_dp,100.0_dp,3.0_dp,4.0_dp]
   prices = price_gb_option(r,te,s0,strikes,y,theta_gb(1),theta_gb(2),theta_gb(3),theta_gb(4))
   value = gb_objective(theta_gb,r,te,y,s0,prices%call,strikes,weights,prices%put, &
      strikes,weights,0.0_dp)
   call assert_true(value < 1.0e-20_dp,"GB objective at prices")

   theta_mln = [0.4_dp,log(95.0_dp),log(110.0_dp),0.16_dp,0.24_dp]
   prices = price_mln_option(r,te,y,strikes,theta_mln(1),theta_mln(2),theta_mln(3), &
      theta_mln(4),theta_mln(5))
   value = mln_objective(theta_mln,r,y,te,exp((y-r)*te)*(theta_mln(1)* &
      exp(theta_mln(2)+0.5_dp*theta_mln(4)**2)+(1.0_dp-theta_mln(1))* &
      exp(theta_mln(3)+0.5_dp*theta_mln(5)**2)),prices%call,strikes,weights, &
      prices%put,strikes,weights,0.0_dp)
   call assert_true(value < 1.0e-20_dp,"MLN objective at prices")

   v = sqrt(exp(sigma*sigma*te)-1.0_dp)
   skew = 3.0_dp*v+v**3
   kurt = 16.0_dp*v**2+15.0_dp*v**4+6.0_dp*v**6+v**8
   theta_ew = [sigma,skew,kurt]
   prices = price_ew_option(r,te,s0,strikes,sigma,y,skew,kurt)
   value = ew_objective(theta_ew,r,y,te,s0,prices%call,strikes,weights,1.0_dp)
   call assert_true(value < 1.0e-20_dp,"EW objective at truth")

   theta_am = [0.2_dp,0.8_dp,log(85.0_dp),log(100.0_dp),log(120.0_dp), &
      0.12_dp,0.18_dp,0.22_dp,0.2_dp,0.5_dp]
   value = mln_am_objective(theta_am,s0,r,te,prices%call,weights,prices%put,weights, &
      strikes,1.0_dp)
   call assert_true(value >= 0.0_dp .and. value < huge(1.0_dp),"AM objective finite")
   print '(a)', 'test_objectives: PASS'
contains
   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      if (.not. condition) then
         print '(a)', trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_objectives
