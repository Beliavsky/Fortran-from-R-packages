program test_special_constraints
    use scip
    implicit none
    type(scip_model_t) :: model
    type(scip_solution) :: sol
    integer :: er, x, z, a, b, d, cidx
    real(dp), parameter :: tol = 1.0e-5_dp

    model = scip_model('quadratic', er)
    x = model%add_var(1.0_dp, lb=0.0_dp, ub=10.0_dp, name='x', ierr=er)
    cidx = model%add_quadratic_cons([x], [x], [1.0_dp], lhs=4.0_dp, &
                                    name='x_squared', ierr=er)
    call require(er == 0 .and. cidx == 1, 'quadratic constraint creation')
    call model%optimize(er)
    sol = model%get_solution()
    call require(trim(model%get_status()) == 'optimal', 'quadratic status')
    call require(abs(sol%x(1) - 2.0_dp) < 1.0e-4_dp, 'quadratic solution')
    call model%free()

    model = scip_model('sos1', er)
    call model%set_objective_sense('maximize', er)
    a = model%add_var(1.0_dp, 0.0_dp, 1.0_dp, name='a', ierr=er)
    b = model%add_var(1.0_dp, 0.0_dp, 1.0_dp, name='b', ierr=er)
    cidx = model%add_sos1_cons([a,b], name='sos1', ierr=er)
    call model%optimize(er)
    sol = model%get_solution()
    call require(abs(sol%objval - 1.0_dp) < tol, 'SOS1 objective')
    call model%free()

    model = scip_model('sos2', er)
    call model%set_objective_sense('maximize', er)
    a = model%add_var(1.0_dp, 0.0_dp, 1.0_dp, name='a', ierr=er)
    b = model%add_var(1.0_dp, 0.0_dp, 1.0_dp, name='b', ierr=er)
    d = model%add_var(1.0_dp, 0.0_dp, 1.0_dp, name='d', ierr=er)
    cidx = model%add_sos2_cons([a,b,d], name='sos2', ierr=er)
    call model%optimize(er)
    sol = model%get_solution()
    call require(abs(sol%objval - 2.0_dp) < tol, 'SOS2 objective')
    call model%free()

    model = scip_model('indicator', er)
    call model%set_objective_sense('maximize', er)
    z = model%add_var(100.0_dp, 0.0_dp, 1.0_dp, 'B', 'z', er)
    x = model%add_var(1.0_dp, 0.0_dp, 10.0_dp, 'C', 'x', er)
    cidx = model%add_indicator_cons(z, [x], [1.0_dp], 2.0_dp, 'indicator', er)
    call require(er == 0, 'indicator constraint creation')
    call model%optimize(er)
    sol = model%get_solution()
    call require(trim(model%get_status()) == 'optimal', 'indicator status')
    call require(abs(sol%objval - 102.0_dp) < tol, 'indicator objective')
    call require(abs(sol%x(z) - 1.0_dp) < tol, 'indicator binary')
    call require(abs(sol%x(x) - 2.0_dp) < tol, 'indicator continuous')
    call model%free()

    print *, 'test_special_constraints: PASS'
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
