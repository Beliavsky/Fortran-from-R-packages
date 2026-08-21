program test_vectors
    use goftest, only : dp, p_ad, q_ad, p_cvm, q_cvm
    implicit none
    real(dp) :: x(3), p(3), q(3)

    x = [0.2_dp, 0.5_dp, 1.0_dp]
    p = p_ad(x)
    if (maxval(abs(p - [p_ad(x(1)), p_ad(x(2)), p_ad(x(3))])) > 1.0e-14_dp) error stop 1
    q = q_ad([0.1_dp, 0.5_dp, 0.9_dp])
    if (maxval(abs(p_ad(q) - [0.1_dp, 0.5_dp, 0.9_dp])) > 2.0e-10_dp) error stop 1

    p = p_cvm(x)
    if (maxval(abs(p - [p_cvm(x(1)), p_cvm(x(2)), p_cvm(x(3))])) > 1.0e-14_dp) error stop 1
    q = q_cvm([0.1_dp, 0.5_dp, 0.9_dp])
    if (maxval(abs(p_cvm(q) - [0.1_dp, 0.5_dp, 0.9_dp])) > 2.0e-9_dp) error stop 1

    print '(a)', 'test_vectors: PASS'
end program test_vectors
