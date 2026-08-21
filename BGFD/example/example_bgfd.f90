program example_bgfd
    use bgfd
    implicit none
    real(dp) :: x, p, q

    x = 1.3_dp
    p = p_bell_w(x, 0.9_dp, 1.6_dp, 0.7_dp)
    q = q_bell_w(0.37_dp, 0.9_dp, 1.6_dp, 0.7_dp)

    print '(a,f12.8)', 'Bell-Weibull density = ', d_bell_w(x,0.9_dp,1.6_dp,0.7_dp)
    print '(a,f12.8)', 'Bell-Weibull cdf     = ', p
    print '(a,f12.8)', 'Bell-Weibull q(.37) = ', q
end program example_bgfd
