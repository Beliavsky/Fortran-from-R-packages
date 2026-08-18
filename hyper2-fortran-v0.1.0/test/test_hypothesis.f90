! SPDX-License-Identifier: GPL-2.0-or-later
program test_hypothesis
    use hyper2
    use test_support
    implicit none
    type(hyper2_model) :: h
    type(hyper2_test_result) :: tr
    real(dp) :: wins(2,2)
    character(len=name_len) :: names(2)
    integer :: fails

    fails = 0
    names = ['a                                                               ', &
             'b                                                               ']
    wins = reshape([0.0_dp,2.0_dp,8.0_dp,0.0_dp],[2,2])
    h = pairwise(wins,names)
    tr = equalp_test(h,startp=[0.8_dp,0.2_dp])
    call check(tr%converged,'equalp test fit',fails)
    call check_close(tr%statistic, 8.0_dp*log(0.8_dp)+2.0_dp*log(0.2_dp)- &
        10.0_dp*log(0.5_dp), 2.0e-9_dp, 'equalp support difference',fails)
    call check(tr%p_value > 0.0_dp .and. tr%p_value < 0.05_dp, 'equalp p-value',fails)

    if (fails /= 0) error stop 1
    write(*,'(a)') 'test_hypothesis: PASS'
end program test_hypothesis
