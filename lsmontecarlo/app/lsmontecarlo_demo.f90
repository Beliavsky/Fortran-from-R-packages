! SPDX-License-Identifier: GPL-3.0-only
! Derived from LSMonteCarlo 1.0 by Mikhail A. Beketov.
! Copyright (C) 2013 Mikhail A. Beketov.
program lsmontecarlo_demo
    use lsmontecarlo, only : amer_put_lsm_av, amer_put_lsm_cv, dp, eu_put_bs, option_result
    implicit none

    type(option_result) :: antithetic
    type(option_result) :: controlled
    real(dp) :: european

    european = eu_put_bs(100.0_dp, 0.20_dp, 105.0_dp, 0.05_dp, 0.0_dp, 1.0_dp)
    antithetic = amer_put_lsm_av(100.0_dp, 0.20_dp, 10000, 50, 105.0_dp, 0.05_dp, 0.0_dp, 1.0_dp, 12345)
    controlled = amer_put_lsm_cv(100.0_dp, 0.20_dp, 10000, 50, 105.0_dp, 0.05_dp, 0.0_dp, 1.0_dp, 12345)

    print '(a,f10.5)', 'European put:                 ', european
    print '(a,f10.5,a,f9.5)', 'American put, antithetic:    ', antithetic%price, ' +/- ', antithetic%standard_error
    print '(a,f10.5,a,f9.5)', 'American put, control var.:  ', controlled%price, ' +/- ', controlled%standard_error
end program lsmontecarlo_demo
