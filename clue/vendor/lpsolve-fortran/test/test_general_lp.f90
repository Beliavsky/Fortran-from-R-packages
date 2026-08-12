program test_general_lp
    use lpsolve
    implicit none
    type(lp_result) :: r
    real(dp) :: c(3), a(2,3), b(2)
    integer :: sense(2)

    c = [1.0_dp, 9.0_dp, 1.0_dp]
    a = reshape([1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 2.0_dp], [2,3])
    b = [9.0_dp, 15.0_dp]
    sense = LP_LE

    call solve_lp(LP_MAX, c, a, sense, b, r)
    if (r%status /= LP_OPTIMAL) error stop 'continuous LP status'
    if (abs(r%objective - 40.5_dp) > 1.0e-9_dp) error stop 'continuous LP objective'
    if (maxval(abs(r%solution - [0.0_dp,4.5_dp,0.0_dp])) > 1.0e-9_dp) then
        error stop 'continuous LP solution'
    end if
    if (maxval(abs(r%duals - [4.5_dp,0.0_dp])) > 1.0e-9_dp) error stop 'LP duals'
    if (abs(r%reduced_costs(1) + 3.5_dp) > 1.0e-9_dp) error stop 'reduced cost 1'
end program test_general_lp
