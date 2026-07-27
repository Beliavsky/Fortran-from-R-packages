! SPDX-License-Identifier: GPL-2.0-or-later
program basic_pricing
   use rnd, only : dp, option_prices, price_bsm_option
   implicit none
   real(dp) :: strikes(5)
   type(option_prices) :: prices
   integer :: i
   strikes = [80.0_dp,90.0_dp,100.0_dp,110.0_dp,120.0_dp]
   prices = price_bsm_option(100.0_dp,strikes,0.03_dp,0.5_dp,0.2_dp,0.01_dp)
   print '(a)', ' strike       call        put'
   do i = 1, size(strikes)
      print '(3f12.5)',strikes(i),prices%call(i),prices%put(i)
   end do
end program basic_pricing
