program basic
    use zigg, only : dp, zsetseed, zrnorm, zrexp, zrunif
    implicit none
    real(dp), allocatable :: x(:)

    call zsetseed(12345)

    x = zrnorm(5)
    print '(a,5(1x,f10.6))', 'normal:     ', x

    x = zrexp(5)
    print '(a,5(1x,f10.6))', 'exponential:', x

    x = zrunif(5)
    print '(a,5(1x,f10.6))', 'uniform:    ', x
end program basic
