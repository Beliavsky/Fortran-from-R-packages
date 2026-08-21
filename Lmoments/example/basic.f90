program basic
    use lmoments, only: dp, lmoments_sample, lmom2normpoly4, qnormpoly
    implicit none
    real(dp), parameter :: x(8) = [1.2_dp, 0.3_dp, 2.7_dp, 1.8_dp, 4.0_dp, 0.9_dp, 3.1_dp, 2.2_dp]
    real(dp) :: lm(4), param(4)
    integer :: info

    call lmoments_sample(x, lm, info)
    if (info /= 0) error stop 'lmoments_sample failed'
    param = lmom2normpoly4(lm)
    print '(a,4f12.6)', 'L-moments: ', lm
    print '(a,f12.6)', 'fitted median quantile: ', qnormpoly(0.5_dp, param)
end program basic
