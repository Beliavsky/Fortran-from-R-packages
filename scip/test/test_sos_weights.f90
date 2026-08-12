program test_sos_weights
    use scip
    implicit none
    type(scip_model_t) :: model
    type(scip_solution) :: sol
    integer :: er, i1, i2, i3, cidx

    model = scip_model('weighted_sos2', er)
    call model%set_objective_sense('maximize', er)
    i1 = model%add_var(1.0_dp, 0.0_dp, 1.0_dp, name='x1', ierr=er)
    i2 = model%add_var(10.0_dp, 0.0_dp, 1.0_dp, name='x2', ierr=er)
    i3 = model%add_var(1.0_dp, 0.0_dp, 1.0_dp, name='x3', ierr=er)
    cidx = model%add_sos2_cons([i1,i2,i3], [1.0_dp,2.0_dp,3.0_dp], 'sos2', er)
    call model%optimize(er)
    sol = model%get_solution()
    call require(sol%available, 'solution available')
    call require(abs(sol%objval - 11.0_dp) < 1.0e-8_dp, 'weighted SOS2 objective')
    call require(sol%x(i2) > 0.5_dp, 'middle variable active')
    call model%free()
    print *, 'test_sos_weights: PASS'
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
