program production_lp
    use linprog
    implicit none
    real(dp) :: c(3), b(3), a(3,3)
    type(linprog_result) :: result
    type(linprog_control) :: control

    c = [1800.0_dp,600.0_dp,600.0_dp]
    b = [40.0_dp,90.0_dp,2500.0_dp]
    a = reshape([0.7_dp,1.5_dp,50.0_dp, 0.35_dp,1.0_dp,12.5_dp, &
                 0.0_dp,3.0_dp,20.0_dp], [3,3])
    control%maximum = .true.
    control%solve_dual = .true.
    call solveLP(c, b, a, result, control)
    print '(a,f12.4)', 'objective = ', result%opt
    print '(a,3f12.4)', 'solution  = ', result%solution
    print '(a,3f12.4)', 'duals     = ', result%con_dual
end program production_lp
