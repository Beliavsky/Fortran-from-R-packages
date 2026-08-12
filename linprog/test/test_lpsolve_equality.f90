program test_lpsolve_equality
    use linprog
    implicit none
    real(dp) :: c(2), b(2), a(2,2)
    character(len=2) :: dir(2)
    type(linprog_result) :: legacy, fixed
    type(linprog_control) :: ctl

    c = [27.0_dp, 9.0_dp]
    b = [8.0_dp, 74.0_dp]
    a = reshape([1.0_dp,1.0_dp, -1.0_dp,1.0_dp], [2,2])
    dir = ['==','<=']
    ctl%maximum = .true.
    call solveLP(c, b, a, legacy, ctl, dir)
    call assert_close(legacy%opt, 1998.0_dp, 1.0e-10_dp, 'legacy equality objective')
    call assert_vec(legacy%solution, [74.0_dp,0.0_dp], 1.0e-10_dp, 'legacy equality solution')

    ctl%use_lpsolve = .true.
    call solveLP(c, b, a, fixed, ctl, dir)
    call assert_true(fixed%status == 0 .and. fixed%lp_status == 0, 'lpSolve equality status')
    call assert_close(fixed%opt, 1404.0_dp, 1.0e-9_dp, 'lpSolve equality objective')
    call assert_vec(fixed%solution, [41.0_dp,33.0_dp], 1.0e-9_dp, 'lpSolve equality solution')
    print *, 'test_lpsolve_equality: PASS'
contains
    subroutine assert_true(ok, msg)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: msg
        if (.not. ok) error stop msg
    end subroutine
    subroutine assert_close(x, y, tol, msg)
        real(dp), intent(in) :: x, y, tol
        character(len=*), intent(in) :: msg
        if (abs(x-y) > tol) error stop msg
    end subroutine
    subroutine assert_vec(x, y, tol, msg)
        real(dp), intent(in) :: x(:), y(:), tol
        character(len=*), intent(in) :: msg
        if (size(x) /= size(y) .or. maxval(abs(x-y)) > tol) error stop msg
    end subroutine
end program test_lpsolve_equality
