! SPDX-License-Identifier: GPL-2.0-or-later
program test_pricing
   use rnd, only : dp, option_prices, am_price_result
   use rnd, only : price_bsm_option, price_mln_option, price_ew_option, price_shimko_option
   use rnd, only : price_am_option
   implicit none
   real(dp), parameter :: tol = 2.0e-11_dp
   real(dp) :: strikes(3), s0, r, y, te, sigma, parity(3)
   type(option_prices) :: bsm, mln, ew, shimko
   type(am_price_result) :: am

   s0 = 100.0_dp
   strikes = [90.0_dp,100.0_dp,110.0_dp]
   r = 0.04_dp
   y = 0.01_dp
   te = 0.5_dp
   sigma = 0.25_dp
   bsm = price_bsm_option(s0,strikes,r,te,sigma,y)
   parity = s0*exp(-y*te)-strikes*exp(-r*te)
   call assert_all_close(bsm%call-bsm%put,parity,tol,"BSM parity")
   call assert_true(all(bsm%call >= 0.0_dp) .and. all(bsm%put >= 0.0_dp),"BSM nonnegative")

   mln = price_mln_option(r,te,y,strikes,1.0_dp, &
      log(s0)+(r-y-0.5_dp*sigma*sigma)*te,log(s0),sigma*sqrt(te),0.4_dp)
   call assert_all_close(mln%call,bsm%call,5.0e-11_dp,"single-component MLN equals BSM")
   call assert_all_close(mln%put,bsm%put,5.0e-11_dp,"single-component MLN put")

   ew = price_ew_option(r,te,s0,strikes,sigma,y, &
      lognormal_skew(sigma,te),lognormal_kurt(sigma,te))
   call assert_all_close(ew%call,bsm%call,2.0e-10_dp,"Edgeworth baseline")

   shimko = price_shimko_option(r,te,s0,strikes,y,sigma,0.0_dp,0.0_dp)
   call assert_all_close(shimko%call,bsm%call,tol,"constant Shimko equals BSM")

   am = price_am_option(100.0_dp,r,te,0.5_dp,log(90.0_dp),log(100.0_dp), &
      log(115.0_dp),0.15_dp,0.20_dp,0.25_dp,0.25_dp,0.50_dp)
   call assert_true(am%call >= 0.0_dp .and. am%put >= 0.0_dp,"AM prices nonnegative")
   call assert_true(am%prob_above >= 0.0_dp .and. am%prob_above <= 1.0_dp,"AM probability")
   print '(a)', 'test_pricing: PASS'
contains
   real(dp) function lognormal_skew(vol, maturity) result(value)
      real(dp), intent(in) :: vol, maturity
      real(dp) :: v
      v = sqrt(exp(vol*vol*maturity)-1.0_dp)
      value = 3.0_dp*v+v**3
   end function lognormal_skew
   real(dp) function lognormal_kurt(vol, maturity) result(value)
      real(dp), intent(in) :: vol, maturity
      real(dp) :: v
      v = sqrt(exp(vol*vol*maturity)-1.0_dp)
      value = 16.0_dp*v**2+15.0_dp*v**4+6.0_dp*v**6+v**8
   end function lognormal_kurt
   subroutine assert_all_close(actual,expected,tolerance,message)
      real(dp), intent(in) :: actual(:), expected(:), tolerance
      character(*), intent(in) :: message
      if (maxval(abs(actual-expected)) > tolerance) then
         print '(a,es14.6)', trim(message)//' error: ',maxval(abs(actual-expected))
         error stop 1
      end if
   end subroutine assert_all_close
   subroutine assert_true(condition,message)
      logical, intent(in) :: condition
      character(*), intent(in) :: message
      if (.not. condition) then
         print '(a)', trim(message)
         error stop 1
      end if
   end subroutine assert_true
end program test_pricing
