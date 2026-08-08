program test_lmm
    use mla_kinds, only : dp
    use mla_lmm, only : loglik_lmm, grad_lmm
    implicit none
    real(dp) :: b(4), y(6), x(6, 2), g(4), gn(4), bp(4), bm(4), h
    integer :: ni(2), i

    b = [1.0_dp, 0.5_dp, 0.8_dp, 1.2_dp]
    y = [1.2_dp, 1.8_dp, 2.1_dp, 0.7_dp, 1.1_dp, 1.6_dp]
    x(:, 1) = 1.0_dp
    x(:, 2) = [0.0_dp, 1.0_dp, 2.0_dp, 0.0_dp, 1.0_dp, 2.0_dp]
    ni = [3, 3]
    call grad_lmm(b, y, x, ni, g)
    do i = 1, 4
        h = 1.0e-6_dp
        bp = b
        bm = b
        bp(i) = bp(i) + h
        bm(i) = bm(i) - h
        gn(i) = (loglik_lmm(bp, y, x, ni) - loglik_lmm(bm, y, x, ni)) / (2.0_dp * h)
    end do
    if (maxval(abs(g - gn)) > 1.0e-6_dp) error stop 1
end program test_lmm
