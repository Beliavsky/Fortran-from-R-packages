program demo_skewt
    use skewt, only : dp, dskt, pskt, qskt, rskt
    implicit none
    real(dp) :: x(5), p(5), r(8)
    integer :: i

    x = [-2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
    p = [0.025_dp, 0.10_dp, 0.50_dp, 0.90_dp, 0.975_dp]

    print '(a)', 'x, density, cdf for df=5, gamma=2'
    print '(3es18.8)', (x(i), dskt(x(i), 5.0_dp, 2.0_dp), &
        pskt(x(i), 5.0_dp, 2.0_dp), i=1,size(x))
    print '(a)', 'quantiles:'
    print '(5es18.8)', qskt(p, 5.0_dp, 2.0_dp)
    call rskt(r, 5.0_dp, 2.0_dp)
    print '(a)', 'random sample:'
    print '(4es18.8)', r
end program demo_skewt
