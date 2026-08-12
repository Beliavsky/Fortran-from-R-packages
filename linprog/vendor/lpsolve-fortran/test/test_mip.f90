program test_mip
    use lpsolve
    implicit none
    type(lp_result) :: r
    type(lp_control) :: ctl
    real(dp) :: c(3), a(2,3), b(2)
    integer :: sense(2), ints(3), bins(3)
    real(dp) :: a2(1,3), b2(1)
    integer :: s2(1), k

    c = [1.0_dp, 9.0_dp, 1.0_dp]
    a = reshape([1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 2.0_dp], [2,3])
    b = [9.0_dp, 15.0_dp]
    sense = LP_LE
    ints = [1,2,3]
    call solve_lp(LP_MAX, c, a, sense, b, r, integer_variables=ints)
    if (r%status /= LP_OPTIMAL) error stop 'integer LP status'
    if (abs(r%objective - 37.0_dp) > 1.0e-9_dp) error stop 'integer LP objective'
    if (maxval(abs(r%solution - [1.0_dp,4.0_dp,0.0_dp])) > 1.0e-9_dp) then
        error stop 'integer LP solution'
    end if

    c = 1.0_dp
    a2 = 1.0_dp
    b2 = 2.0_dp
    s2 = LP_LE
    bins = [1,2,3]
    ctl%num_binary_solutions = 3
    call solve_lp(LP_MAX, c, a2, s2, b2, r, ctl, binary_variables=bins)
    if (r%status /= LP_OPTIMAL) error stop 'binary multi status'
    if (r%solution_count /= 3) error stop 'binary multi count'
    if (maxval(abs(r%solution_objectives - 2.0_dp)) > 1.0e-9_dp) then
        error stop 'binary multi objectives'
    end if
    do k = 1, 3
        if (abs(sum(r%solutions(:,k)) - 2.0_dp) > 1.0e-9_dp) then
            error stop 'binary multi solution'
        end if
    end do
end program test_mip
