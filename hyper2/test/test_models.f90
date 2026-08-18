! SPDX-License-Identifier: GPL-2.0-or-later
program test_models
    use hyper2
    use test_support
    implicit none
    type(hyper2_model) :: h
    type(hyper3_model) :: h3, hm
    type(fit_result) :: fit
    real(dp) :: wins(2,2), p(2), ll3, m(5,5), pp(5), em
    real(dp), allocatable :: z(:)
    character(len=name_len) :: n2(2), v3(3), n5(5)
    integer :: fails

    fails = 0
    n2 = ['a                                                               ', &
          'b                                                               ']
    wins = reshape([0.0_dp,2.0_dp,8.0_dp,0.0_dp],[2,2])
    h = pairwise(wins,n2)
    fit = maxp(h, startp=[0.8_dp,0.2_dp])
    call check(fit%converged, 'pairwise maxp converged', fails)
    call check_vec(fit%p, [0.8_dp,0.2_dp], 1.0e-7_dp, 'pairwise MLE', fails)
    z = zermelo(wins, maxit=1000, tol=1.0e-12_dp)
    call check_vec(z, fit%p, 1.0e-7_dp, 'zermelo agrees with MLE', fails)

    v3 = ['a                                                               ', &
          'b                                                               ', &
          'a                                                               ']
    h3 = ordervec2supp3(v3)
    p = [0.6_dp,0.4_dp]
    ll3 = log(0.6_dp) - log(2.0_dp*0.6_dp+0.4_dp) + &
        log(0.4_dp) - log(0.6_dp+0.4_dp)
    call check_close(loglik(p,h3), ll3, 2.0e-13_dp, 'repeated-player hyper3 rank likelihood', fails)


    n5 = ['a                                                               ', &
          'b                                                               ', &
          'c                                                               ', &
          'x                                                               ', &
          'y                                                               ']
    m = 0.0_dp
    m(1,:) = [1.2_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
    m(2,:) = [1.2_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    m(3,:) = [9.0_dp,3.0_dp,1.0_dp,0.0_dp,0.0_dp]
    m(4,:) = [0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp]
    m(5,:) = [0.0_dp,0.0_dp,0.0_dp,6.0_dp,8.0_dp]
    hm = hyper3_matrix(m,[3.0_dp,4.0_dp,5.0_dp,-7.0_dp,-5.0_dp],n5)
    pp = [0.2_dp,0.1_dp,0.15_dp,0.25_dp,0.3_dp]
    em = 3.0_dp*log(1.2_dp*pp(1)) + 4.0_dp*log(1.2_dp*pp(1)+pp(2)) + &
        5.0_dp*log(9.0_dp*pp(1)+3.0_dp*pp(2)+pp(3)) - 7.0_dp*log(pp(2)) - &
        5.0_dp*log(6.0_dp*pp(4)+8.0_dp*pp(5))
    call check_close(loglik(pp,hm),em,3.0e-13_dp,'hyper3 matrix constructor',fails)

    if (fails /= 0) error stop 1
    write(*,'(a)') 'test_models: PASS'
end program test_models
