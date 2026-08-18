! SPDX-License-Identifier: GPL-2.0-or-later
program test_sums
    use hyper2
    use test_support
    implicit none
    type(name_group) :: groups(2)
    type(hyper2_list) :: s, s2
    type(lsl_model) :: lx
    type(fit_result) :: fit
    real(dp) :: p(3), expected
    integer :: fails

    fails = 0
    allocate(groups(1)%items(2),groups(2)%items(1))
    groups(1)%items = ['a                                                               ', &
                       'b                                                               ']
    groups(2)%items = ['c                                                               ']
    s = general_grouped_rank_likelihood(groups)
    call check(s%size() == 2, 'two tied first-place permutations', fails)
    p = [0.5_dp,0.3_dp,0.2_dp]
    expected = 0.5_dp*(0.3_dp/0.5_dp) + 0.3_dp*(0.5_dp/0.7_dp)
    call check_close(like_single_list(p,s), expected, 3.0e-13_dp, &
        'grouped-rank likelihood sum', fails)
    s2 = suplist_scale(s,2)
    call check_close(like_single_list(p,s2), expected*expected, 1.0e-12_dp, &
        'suplist repeated-product scaling', fails)
    lx = lsl_create()
    call lsl_add_component(lx,s,power=2.0_dp)
    call check_close(loglik_lsl(p,lx), 2.0_dp*log(expected), 1.0e-12_dp, &
        'lsl log likelihood', fails)
    fit = maxp_lsl(lx, startp=[0.4_dp,0.4_dp,0.2_dp], max_iter=400)
    call check(abs(sum(fit%p)-1.0_dp) < 1.0e-12_dp .and. all(fit%p>0.0_dp), &
        'lsl optimizer simplex', fails)

    if (fails /= 0) error stop 1
    write(*,'(a)') 'test_sums: PASS'
end program test_sums
