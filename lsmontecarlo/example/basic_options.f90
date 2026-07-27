! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
program basic_options
    use lsmontecarlo, only : amer_put_lsm, asian_amer_put_lsm, dp, option_result
    implicit none

    type(option_result) :: american
    type(option_result) :: asian

    american = amer_put_lsm(spot=100.0_dp, sigma=0.25_dp, n=12000, m=60, strike=105.0_dp, &
        rate=0.04_dp, dividend=0.01_dp, maturity=1.0_dp, seed=777)
    asian = asian_amer_put_lsm(spot=100.0_dp, sigma=0.25_dp, n=12000, m=60, strike=105.0_dp, &
        rate=0.04_dp, dividend=0.01_dp, maturity=1.0_dp, seed=777)

    print '(a,f10.5)', 'American put price:       ', american%price
    print '(a,f10.5)', 'Asian American put price: ', asian%price
end program basic_options
