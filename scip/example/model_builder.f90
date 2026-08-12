program model_builder
    use scip
    implicit none
    type(scip_model_t) :: model
    type(scip_solution) :: sol
    integer :: er, first, cidx
    real(dp) :: obj(3), lb(3), ub(3)
    character(len=1) :: vt(3)

    model = scip_model('choose_one', er)
    if (er /= 0) error stop 'could not create SCIP model'
    call model%set_objective_sense('maximize')

    obj = [3.0_dp, 2.0_dp, 1.0_dp]
    lb = 0.0_dp
    ub = 1.0_dp
    vt = 'B'
    first = model%add_vars(obj, lb, ub, vt)
    cidx = model%add_linear_cons([1,2,3], [1.0_dp,1.0_dp,1.0_dp], &
                                 lhs=1.0_dp, rhs=1.0_dp)
    if (first /= 1 .or. cidx /= 1) error stop 'unexpected model indices'

    call model%optimize()
    sol = model%get_solution()
    print '(a,a)', 'status: ', trim(model%get_status())
    print '(a,f8.3)', 'objective: ', sol%objval
    print '(a,3f8.3)', 'x: ', sol%x
    call model%free()
end program
