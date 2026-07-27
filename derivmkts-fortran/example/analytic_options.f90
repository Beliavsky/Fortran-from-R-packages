! SPDX-License-Identifier: MIT
! Copyright (c) 2014-2022 Robert L. McDonald
program analytic_options
    use derivmkts, only: dp,bs_call_greeks,greek_result
    use derivmkts, only: calloncall,compound_result,callperpetual,perpetual_result
    implicit none
    type(greek_result) :: g
    type(compound_result) :: compound
    type(perpetual_result) :: perpetual
    g=bs_call_greeks(100.0_dp,105.0_dp,0.25_dp,0.04_dp,0.75_dp,0.01_dp)
    compound=calloncall(100.0_dp,105.0_dp,7.0_dp,0.25_dp,0.04_dp,0.25_dp,0.75_dp,0.01_dp)
    perpetual=callperpetual(100.0_dp,100.0_dp,0.25_dp,0.05_dp,0.02_dp)
    print '(a,f12.6)', 'Call price: ',g%premium
    print '(a,f12.6)', 'Delta:      ',g%delta
    print '(a,f12.6)', 'Gamma:      ',g%gamma
    print '(a,f12.6)', 'Compound:   ',compound%price
    print '(a,f12.6)', 'Perpetual:  ',perpetual%price
end program analytic_options
