program test_ripple
    use subplex, only: dp, subplex_result, subplex_minimize
    implicit none

    type(subplex_result) :: result
    real(dp) :: x1(1), x2(2), x3(3)

    x1 = [1.0_dp]
    call subplex_minimize(ripple, x1, result, compute_hessian=.true.)
    call assert_true(result%convergence == 0, "1-D ripple convergence")
    call assert_true(abs(result%x(1)) <= 1.0e-5_dp, "1-D ripple solution")
    call assert_true(abs(result%value) <= 1.0e-10_dp, "1-D ripple value")
    call assert_true(abs(result%hessian(1,1)-202.0_dp) <= 2.0e-4_dp, &
        "1-D ripple Hessian")

    x2 = [0.1_dp, 3.0_dp]
    call subplex_minimize(ripple, x2, result, compute_hessian=.true.)
    call assert_true(abs(result%value) <= 1.0e-5_dp, "2-D ripple value")

    x3 = [0.1_dp, 3.0_dp, 2.0_dp]
    call subplex_minimize(ripple, x3, result)
    call assert_close_vec(result%x, [0.45932_dp, 1.10399_dp, 0.34408_dp], &
        1.0e-4_dp, "3-D ripple local minimum")

contains

    function ripple(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f, r
        r = sqrt(sum(x*x))
        f = 1.0_dp - exp(-r*r)*cos(10.0_dp*r)**2
    end function ripple

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

end program test_ripple
