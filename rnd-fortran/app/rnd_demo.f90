! SPDX-License-Identifier: GPL-2.0-or-later
program rnd_demo
   use rnd, only : dp, option_prices, bsm_fit
   use rnd, only : price_bsm_option, extract_bsm_density
   implicit none
   real(dp), parameter :: s0=100.0_dp, r=0.04_dp, y=0.01_dp, te=0.5_dp, sigma=0.25_dp
   real(dp) :: strikes(7), initial(2)
   type(option_prices) :: market
   type(bsm_fit) :: fit
   integer :: i

   do i = 1, size(strikes)
      strikes(i) = 85.0_dp+5.0_dp*real(i-1,dp)
   end do
   market = price_bsm_option(s0,strikes,r,te,sigma,y)
   initial = [log(s0)+(r-y-0.5_dp*0.2_dp**2)*te,0.2_dp*sqrt(te)]
   fit = extract_bsm_density(r,y,te,s0,market%call,strikes,market%put,strikes, &
      initial_values=initial,max_iter=3000)
   print '(a,f10.6)', 'fitted annual volatility: ',fit%zeta/sqrt(te)
   print '(a,es14.6)', 'objective: ',fit%optimizer%value
   print '(a,i0)', 'convergence code: ',fit%optimizer%convergence
end program rnd_demo
