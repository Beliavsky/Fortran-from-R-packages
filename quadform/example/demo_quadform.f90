program demo_quadform
    use quadform, only : dp, qf, qd, qtr
    implicit none
    real(dp) :: m(2,2), x(2,3)

    m = reshape([4.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2,2])
    x = reshape([1.0_dp,2.0_dp, 2.0_dp,-1.0_dp, -1.0_dp,0.5_dp], [2,3])

    print '(a)', 'x^T M x ='
    print '(3f12.5)', qf(m,x)
    print '(a,3f12.5)', 'diagonal = ', qd(m,x)
    print '(a,f12.5)', 'trace    = ', qtr(m,x)
end program demo_quadform
