! SPDX-License-Identifier: MIT
program credit_curve
   use credule, only : dp, price_cds, bootstrap_cds, credule_success
   implicit none
   real(dp) :: yield_tenors(4), yield_rates(4), maturities(5), market_spreads(5)
   real(dp) :: survival(5), hazards(5), repriced(5)
   integer :: status, i

   yield_tenors = [0.5_dp, 2.0_dp, 5.0_dp, 10.0_dp]
   yield_rates = [0.015_dp, 0.020_dp, 0.025_dp, 0.030_dp]
   maturities = [1.0_dp, 2.0_dp, 3.0_dp, 5.0_dp, 10.0_dp]
   market_spreads = [0.006_dp, 0.009_dp, 0.012_dp, 0.016_dp, 0.021_dp]

   ! Defaults: quarterly premiums, monthly default integration, accrued premiums.
   call bootstrap_cds(yield_tenors, yield_rates, maturities, market_spreads, 0.4_dp, &
                      survival, hazards, status)
   if (status /= credule_success) error stop 'Calibration failed'
   call price_cds(yield_tenors, yield_rates, maturities, survival, maturities, &
                  0.4_dp, repriced, status)
   if (status /= credule_success) error stop 'Repricing failed'
   if (maxval(abs(repriced - market_spreads)) > 1.0e-12_dp) error stop 'Repricing residual is too large'
   print '(a)', 'Years  Market (bp)  Repriced (bp)  Survival  Interval hazard'
   do i = 1, size(maturities)
      print '(f5.1,2f14.6,2f12.7)', maturities(i), 10000.0_dp * market_spreads(i), &
                                  10000.0_dp * repriced(i), survival(i), hazards(i)
   end do
end program credit_curve
