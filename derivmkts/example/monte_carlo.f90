! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program monte_carlo
    use derivmkts, only: dp,arithasianmc,arithavgpricecv,asian_mc_result
    implicit none
    type(asian_mc_result) :: crude,cv
    crude=arithasianmc(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,12,50000,123)
    cv=arithavgpricecv(40.0_dp,40.0_dp,0.30_dp,0.08_dp,0.25_dp,0.0_dp,12,50000,123)
    print '(a,f12.6)', 'Arithmetic Asian call: ',crude%avg_price_call
    print '(a,f12.6)', 'Control-variate call:  ',cv%avg_price_call
    print '(a,f12.6)', 'Control beta:          ',cv%beta
    print '(a,f12.6)', 'Crude payoff SD:       ',crude%sd_avg_price_call
    print '(a,f12.6)', 'Corrected payoff SD:   ',cv%sd_avg_price_call
end program monte_carlo
