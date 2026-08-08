program test_controls
    use subplex, only: dp, subplex_control, subplex_result, subplex_minimize
    implicit none

    type(subplex_control) :: control
    type(subplex_result) :: result
    real(dp) :: x(3)

    x = [0.1_dp, 3.0_dp, 2.0_dp]
    control%maxeval = 2
    call subplex_minimize(ripple, x, result, control)
    call assert_true(result%convergence == -1, "maxeval convergence code")
    call assert_true(result%counts >= 2, "maxeval count semantics")

    control%maxeval = 10000
    control%reltol = 1.0e-3_dp
    call subplex_minimize(ripple, x, result, control)
    call assert_true(abs(result%value-0.7906068_dp) <= 1.0e-3_dp, &
        "looser relative tolerance")

    control%reltol = epsilon(1.0_dp)
    allocate(control%parscale(1))
    control%parscale = 1.0e-90_dp
    call subplex_minimize(ripple, x, result, control)
    call assert_true(result%convergence == -2, "tiny scalar parscale")

    deallocate(control%parscale)
    allocate(control%parscale(3))
    control%parscale = [0.01_dp, 0.05_dp, 0.1_dp]
    call subplex_minimize(ripple, x, result, control)
    call assert_true(result%convergence == 0 .or. result%convergence == 1, &
        "vector parscale")

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

end program test_controls
