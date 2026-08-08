program test_numerical_hessian
    use subplex, only: dp, numerical_hessian
    implicit none

    real(dp) :: x(3), h(3), hess(3,3), expected(3,3)

    x = [1.0_dp, -2.0_dp, 0.5_dp]
    h = 1.0e-4_dp
    expected = 0.0_dp
    expected(1,1) = 2.0_dp
    expected(2,2) = 4.0_dp
    expected(3,3) = 6.0_dp
    call numerical_hessian(quadratic, x, h, hess)
    call assert_true(maxval(abs(hess-expected)) <= 1.0e-6_dp, &
        "numerical Hessian")

contains

    function quadratic(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = x(1)**2 + 2.0_dp*x(2)**2 + 3.0_dp*x(3)**2
    end function quadratic

    subroutine assert_true(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            write(*,'(a)') "FAIL: "//label
            error stop 1
        end if
    end subroutine assert_true

end program test_numerical_hessian
