program test_status_control
    use scip
    implicit none
    type(scip_model_t) :: model
    type(scip_control) :: ctrl
    type(scip_result) :: res
    real(dp) :: a(1,1), obj(1), b(1)
    character(len=2) :: sense(1)
    integer :: er, x, cidx

    model = scip_model('unbounded', er)
    x = model%add_var(-1.0_dp, lb=0.0_dp, name='x', ierr=er)
    call model%optimize(er)
    call require(trim(model%get_status()) == 'unbounded', 'unbounded status')
    call model%free()

    model = scip_model('infeasible', er)
    x = model%add_var(0.0_dp, 0.0_dp, 1.0_dp, 'B', 'x', er)
    cidx = model%add_linear_cons([x], [1.0_dp], lhs=2.0_dp, name='impossible', ierr=er)
    call model%optimize(er)
    call require(trim(model%get_status()) == 'infeasible', 'infeasible status')
    call model%free()

    ctrl%verbose = .false.
    ctrl%presolving = .false.
    ctrl%gap_limit = 0.01_dp
    call ctrl%add_param('display/freq', -1)
    a(1,1) = 1.0_dp
    obj = [1.0_dp]
    b = [2.0_dp]
    sense = ['>=']
    res = scip_solve(obj, a, b, sense, control=ctrl)
    call require(trim(res%status) == 'optimal', 'control solve status')
    call require(abs(res%objval - 2.0_dp) < 1.0e-9_dp, 'control solve objective')

    print *, 'test_status_control: PASS'
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
