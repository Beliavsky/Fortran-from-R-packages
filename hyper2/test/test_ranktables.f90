! SPDX-License-Identifier: GPL-2.0-or-later
program test_ranktables
    use hyper2
    use test_support
    implicit none
    integer :: ranks(2,3), order(3,2)
    integer, allocatable :: o(:,:), r(:,:)
    type(hyper2_model) :: h
    character(len=name_len) :: names(3)
    real(dp) :: p(3), expected
    integer :: fails

    fails = 0
    ranks = reshape([1,3,2,1,3,2],[2,3])
    o = ranktable_to_ordertable(ranks)
    r = ordertable_to_ranktable(o)
    call check(all(r == ranks), 'rank/order table round trip', fails)

    order(:,1) = [1,2,3]
    order(:,2) = [2,1,3]
    names = ['a                                                               ', &
             'b                                                               ', &
             'c                                                               ']
    h = ordertable2supp(order,names)
    p = [0.5_dp,0.3_dp,0.2_dp]
    expected = log(0.5_dp*0.3_dp/0.5_dp) + log(0.3_dp*0.5_dp/0.7_dp)
    call check_close(loglik(p,h), expected, 3.0e-13_dp, 'ordertable support', fails)

    if (fails /= 0) error stop 1
    write(*,'(a)') 'test_ranktables: PASS'
end program test_ranktables
