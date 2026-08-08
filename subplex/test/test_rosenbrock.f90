program test_rosenbrock
    use subplex, only: dp, subplex_control, subplex_result, subplex_minimize
    implicit none

    type(subplex_control) :: control
    type(subplex_result) :: result
    real(dp) :: x2(2), x6(6)

    x2 = [11.0_dp, -33.0_dp]
    call subplex_minimize(rosen, x2, result)
    call assert_true(result%convergence == 0, "2-D Rosenbrock convergence")
    call assert_close_vec(result%x, [1.0_dp, 1.0_dp], 1.0e-4_dp, &
        "2-D Rosenbrock solution")

    x6 = [-33.0_dp, 11.0_dp, 14.0_dp, 9.0_dp, 0.0_dp, 12.0_dp]
    control%maxeval = 30000
    call subplex_minimize(rosen_pairs, x6, result, control)
    call assert_true(result%convergence == 0 .or. result%convergence == 1, &
        "6-D Rosenbrock convergence")
    call assert_close_vec(result%x, [1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
        1.0_dp, 1.0_dp], 1.0e-4_dp, "6-D Rosenbrock solution")

contains

    function rosen(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        f = 100.0_dp*(x(2)-x(1)*x(1))**2 + (1.0_dp-x(1))**2
    end function rosen

    function rosen_pairs(x) result(f)
        real(dp), intent(in) :: x(:)
        real(dp) :: f
        integer :: i
        f = 0.0_dp
        do i = 1, size(x), 2
            f = f + 100.0_dp*(x(i+1)-x(i)*x(i))**2 + (1.0_dp-x(i))**2
        end do
    end function rosen_pairs

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
        call assert_true(size(actual) == size(expected), label//" size")
        call assert_true(maxval(abs(actual-expected)) <= tol, label)
    end subroutine assert_close_vec

end program test_rosenbrock
