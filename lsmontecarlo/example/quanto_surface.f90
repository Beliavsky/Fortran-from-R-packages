! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
program quanto_surface
    use lsmontecarlo, only : dp, price_surface, quanto_amer_put_lsm_price_surface
    use lsmontecarlo, only : surface_maximum, surface_mean, surface_minimum
    implicit none

    type(price_surface) :: surface
    real(dp) :: strikes(3)
    real(dp) :: volatilities(3)

    volatilities = [0.15_dp, 0.25_dp, 0.35_dp]
    strikes = [90.0_dp, 100.0_dp, 110.0_dp]
    surface = quanto_amer_put_lsm_price_surface(volatilities, strikes, spot=100.0_dp, n=4000, m=30, &
        rate=0.04_dp, maturity=1.0_dp, spot2=1.10_dp, sigma2=0.15_dp, rho=-0.25_dp, seed=2468)

    print '(a,f10.5)', 'Surface mean: ', surface_mean(surface)
    print '(a,f10.5)', 'Surface min:  ', surface_minimum(surface)
    print '(a,f10.5)', 'Surface max:  ', surface_maximum(surface)
end program quanto_surface
