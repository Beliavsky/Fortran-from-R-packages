program test_model_api
    use scip
    implicit none
    type(scip_model_t) :: model
    type(scip_solution) :: sol, sol1
    type(scip_info) :: info
    integer :: er, first, cidx
    real(dp) :: obj(3), lb(3), ub(3)
    character(len=1) :: vt(3)
    character(len=4) :: names(3)

    model = scip_model('model_api', er)
    call require(er == 0, 'model creation')
    call model%set_objective_sense('maximize', er)
    call require(er == 0, 'objective sense')
    call model%set_param('display/verblevel', 0, er)
    call require(er == 0, 'integer parameter')
    call model%set_param('limits/time', 30.0_dp, er)
    call require(er == 0, 'real parameter')
    call model%set_param('branching/scorefunc', 's', er)
    call require(er == 0, 'text parameter')

    obj = [3.0_dp, 2.0_dp, 1.0_dp]
    lb = 0.0_dp
    ub = 1.0_dp
    vt = 'B'
    names = ['x1  ', 'x2  ', 'x3  ']
    first = model%add_vars(obj, lb, ub, vt, names, er)
    call require(er == 0 .and. first == 1, 'add_vars')
    call require(model%get_nvars() == 3, 'variable count')

    cidx = model%add_linear_cons([1,2,3], [1.0_dp,1.0_dp,1.0_dp], &
                                 lhs=1.0_dp, rhs=1.0_dp, name='choose_one', ierr=er)
    call require(er == 0 .and. cidx == 1, 'linear constraint')
    call require(model%get_nconss() == 1, 'constraint count')

    call model%optimize(er)
    call require(er == 0, 'optimize')
    call require(trim(model%get_status()) == 'optimal', 'status')
    sol = model%get_solution()
    call require(sol%available, 'best solution available')
    call require(abs(sol%objval - 3.0_dp) < 1.0e-9_dp, 'objective')
    call require(maxloc(sol%x, dim=1) == 1, 'best variable')
    call require(model%get_nsols() >= 1, 'solution count')
    sol1 = model%get_sol(1)
    call require(abs(sol1%objval - sol%objval) < 1.0e-12_dp, 'solution pool')
    info = model%get_info()
    call require(info%sol_count >= 1, 'info solution count')
    call require(info%solve_time >= 0.0_dp, 'solve time')

    call model%free(er)
    call require(er == 0, 'free')
    print *, 'test_model_api: PASS'
contains
    subroutine require(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            print *, 'FAIL: ', trim(label)
            error stop 1
        end if
    end subroutine
end program
