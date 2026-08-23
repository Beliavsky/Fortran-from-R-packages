program basic
    use ecpdist, only: dp, decp, pecp, qecp, ecp_integral_result, ecp_kmoment
    implicit none
    type(ecp_integral_result) :: mean_result

    mean_result = ecp_kmoment(1, 1.0_dp, 1.0_dp, 1.0_dp)
    print '(a,f12.8)', 'density at 1 = ', decp(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
    print '(a,f12.8)', 'cdf at 1     = ', pecp(1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp)
    print '(a,f12.8)', 'median       = ', qecp(0.5_dp, 1.0_dp, 1.0_dp, 1.0_dp)
    print '(a,f12.8)', 'mean         = ', mean_result%estimate
end program basic
