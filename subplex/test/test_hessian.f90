program test_hessian
    use subplex, only: dp, subplex_control, subplex_result, subplex_minimize
    implicit none

    type(subplex_control) :: control
    type(subplex_result) :: result
    real(dp) :: x(2), expected(2,2)

    x = [11.0_dp, -33.0_dp]
    allocate(control%parscale(1))
    control%parscale = 5.0_dp
    call subplex_minimize(shifted_rosen, x, result, control, .true.)

    expected = reshape([200.0_dp, -400.0_dp, -400.0_dp, 802.0_dp], [2,2])
    call assert_close_vec(result%x, [23.0_dp, 1.0_dp], 1.0e-5_dp, &
        "shifted Rosenbrock solution")
    call assert_true(maxval(abs(result%hessian-expected)) <= 2.0e-4_dp, &
        "shifted Rosenbrock Hessian")

contains

    function shifted_rosen(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f, a, b
        b = x(1) - 22.0_dp
        a = x(2)
        f = 100.0_dp*(b-a*a)**2 + (1.0_dp-a)**2
    end function shifted_rosen

    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            write(*,'(a)') "FAIL: "//label
            error stop 1
        end if
    end subroutine assert_true

    subroutine assert_close_vec(actual, expected, tol, label)
        real(dp), intent(in) :: actual(:), expected(:), tol
        character(len=*), intent(in) :: label
        call assert_true(maxval(abs(actual-expected)) <= tol, label)
    end subroutine assert_close_vec

end program test_hessian
