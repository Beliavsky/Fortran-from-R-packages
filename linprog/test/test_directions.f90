program test_directions
    use linprog
    implicit none
    real(dp) :: c(2), b(3), a(3,2)
    character(len=2) :: dir(3)
    type(linprog_result) :: r1, r2
    type(linprog_control) :: ctl

    c = [2.5_dp,2.0_dp]
    b = [10.0_dp,1.5_dp,12.0_dp]
    a = reshape([1.6_dp,0.5_dp,2.0_dp, 2.4_dp,0.2_dp,2.0_dp], [3,2])
    dir = ['>=','>=','<=']
    call solveLP(c, b, a, r1, ctl, dir)
    ctl%use_lpsolve = .true.
    call solveLP(c, b, a, r2, ctl, dir)
    if (r1%status /= 0 .or. r2%status /= 0) error stop 'direction status'
    if (abs(r1%opt-10.454545455_dp) > 1.0e-9_dp) error stop 'direction legacy objective'
    if (abs(r2%opt-10.4545454545_dp) > 1.0e-8_dp) error stop 'direction lpsolve objective'
    if (maxval(abs(r1%solution-r2%solution)) > 1.0e-7_dp) error stop 'direction solutions'
    print *, 'test_directions: PASS'
end program test_directions
