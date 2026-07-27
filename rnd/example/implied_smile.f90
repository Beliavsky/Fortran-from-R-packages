! SPDX-License-Identifier: GPL-2.0-or-later
program implied_smile
   use rnd, only : dp, option_prices, shimko_fit
   use rnd, only : price_bsm_option, extract_shimko_density
   implicit none
   real(dp) :: strikes(9), local_volatility
   real(dp), allocatable :: calls(:)
   type(option_prices) :: one
   type(shimko_fit) :: fit
   integer :: i

   allocate(calls(size(strikes)))
   do i = 1, size(strikes)
      strikes(i) = 80.0_dp+5.0_dp*real(i-1,dp)
      local_volatility = 0.20_dp+0.00002_dp*(strikes(i)-100.0_dp)**2
      one = price_bsm_option(100.0_dp,strikes(i:i),0.03_dp,0.5_dp,local_volatility,0.01_dp)
      calls(i) = one%call(1)
   end do
   fit = extract_shimko_density(calls,strikes,0.03_dp,0.01_dp,0.5_dp,100.0_dp)
   print '(a,3es14.6)', 'quadratic volatility coefficients: ',fit%a0,fit%a1,fit%a2
end program implied_smile
