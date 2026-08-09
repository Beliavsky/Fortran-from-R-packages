program basic_lp
    use lpsolve
    implicit none
    type(lp_result) :: result
    real(dp) :: objective(3), a(2,3), rhs(2)
    integer :: sense(2), integers(3)

    objective = [1.0_dp, 9.0_dp, 1.0_dp]
    a = reshape([1.0_dp, 3.0_dp, 2.0_dp, 2.0_dp, 3.0_dp, 2.0_dp], [2,3])
    rhs = [9.0_dp, 15.0_dp]
    sense = LP_LE

    call solve_lp(LP_MAX, objective, a, sense, rhs, result)
    print '(a,f12.6)', 'LP objective:  ', result%objective
    print '(a,3f12.6)', 'LP solution:   ', result%solution

    integers = [1,2,3]
    call solve_lp(LP_MAX, objective, a, sense, rhs, result, integer_variables=integers)
    print '(a,f12.6)', 'MILP objective:', result%objective
    print '(a,3f12.6)', 'MILP solution: ', result%solution
end program basic_lp
