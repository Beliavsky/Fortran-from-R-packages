program lpsolve_equality
    use linprog
    implicit none
    real(dp) :: c(2), b(2), a(2,2)
    character(len=2) :: dir(2)
    type(linprog_result) :: result
    type(linprog_control) :: control

    c = [27.0_dp,9.0_dp]
    b = [8.0_dp,74.0_dp]
    a = reshape([1.0_dp,1.0_dp,-1.0_dp,1.0_dp],[2,2])
    dir = ['==','<=']
    control%maximum = .true.
    control%use_lpsolve = .true.
    call solveLP(c, b, a, result, control, dir)
    print '(a,f12.4)', 'objective = ', result%opt
    print '(a,2f12.4)', 'solution  = ', result%solution
end program lpsolve_equality
