! SPDX-License-Identifier: GPL-2.0-or-later
program test_sampling
    use hyper2
    use test_support
    implicit none
    real(dp), allocatable :: x(:,:)
    real(dp) :: means(3)
    integer :: fails

    fails = 0
    call seed_rng(271828)
    x = rdirichlet(30000,[1.0_dp,2.0_dp,3.0_dp])
    means = sum(x,dim=1)/real(size(x,1),dp)
    call check_vec(means,[1.0_dp/6.0_dp,1.0_dp/3.0_dp,0.5_dp],1.5e-2_dp, &
        'Dirichlet RNG mean',fails)
    call check(maxval(abs(sum(x,dim=2)-1.0_dp)) < 2.0e-14_dp, 'Dirichlet RNG simplex',fails)

    if (fails /= 0) error stop 1
    write(*,'(a)') 'test_sampling: PASS'
end program test_sampling
