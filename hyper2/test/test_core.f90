! SPDX-License-Identifier: GPL-2.0-or-later
program test_core
    use hyper2
    use test_support
    implicit none
    type(hyper2_model) :: h
    real(dp), allocatable :: g(:), gf(:), hs(:,:)
    real(dp) :: p(3), expected
    integer :: fails
    character(len=name_len) :: names(3)

    fails = 0
    names = ['a                                                               ', &
             'b                                                               ', &
             'c                                                               ']
    h = dirichlet(names, alpha=[1.0_dp,2.0_dp,3.0_dp])
    p = [0.2_dp,0.3_dp,0.5_dp]
    expected = log(0.3_dp) + 2.0_dp*log(0.5_dp)
    call check_close(loglik(p,h), expected, 2.0e-14_dp, 'dirichlet loglik', fails)
    gf = gradient_full(h,p)
    call check_vec(gf, [-3.0_dp, 1.0_dp/0.3_dp-3.0_dp, 4.0_dp-3.0_dp], &
        2.0e-13_dp, 'full gradient', fails)
    g = gradient(h,p)
    call check_vec(g, [gf(1)-gf(3),gf(2)-gf(3)], 2.0e-13_dp, 'simplex gradient', fails)
    hs = hessian_independent(h,p)
    call check(all(abs(hs-transpose(hs)) < 1.0e-13_dp), 'symmetric Hessian', fails)
    call check_close(power_sum_h2(h), 0.0_dp, 1.0e-14_dp, 'homogeneous power sum', fails)

    if (fails /= 0) error stop 1
    write(*,'(a)') 'test_core: PASS'
end program test_core
