! SPDX-License-Identifier: GPL-2.0-or-later
program test_integration
    use hyper2
    use test_support
    implicit none
    type(hyper2_model) :: h
    real(dp), allocatable :: mu(:)
    real(dp) :: z
    character(len=name_len) :: names(3)
    integer :: fails

    fails = 0
    names = ['a                                                               ', &
             'b                                                               ', &
             'c                                                               ']
    h = dirichlet(names, alpha=[1.0_dp,2.0_dp,3.0_dp])
    z = hyper2_B(h, rel_tol=2.0e-8_dp, abs_tol=1.0e-11_dp)
    call check_close(z, 1.0_dp/60.0_dp, 2.0e-6_dp, 'Dirichlet normalizer', fails)
    mu = mean_hyper2(h, rel_tol=5.0e-7_dp, abs_tol=1.0e-10_dp)
    call check_vec(mu, [1.0_dp/6.0_dp,1.0_dp/3.0_dp,0.5_dp], 2.0e-5_dp, &
        'Dirichlet mean', fails)

    if (fails /= 0) error stop 1
    write(*,'(a)') 'test_integration: PASS'
end program test_integration
