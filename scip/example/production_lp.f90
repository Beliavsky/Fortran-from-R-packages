program production_lp
    use scip
    implicit none
    real(dp) :: a(2,2), obj(2), b(2)
    character(len=2) :: sense(2)
    type(scip_control) :: ctrl
    type(scip_result) :: res

    a = reshape([1.0_dp, 2.0_dp, 2.0_dp, 1.0_dp], [2,2])
    obj = [-5.0_dp, -4.0_dp]
    b = [6.0_dp, 8.0_dp]
    sense = ['<=', '<=']
    ctrl%verbose = .false.

    res = scip_solve(obj, a, b, sense, control=ctrl)
    print '(a,a)', 'status: ', trim(res%status)
    print '(a,f10.4)', 'maximized profit: ', -res%objval
    print '(a,2f12.6)', 'x: ', res%x
end program
