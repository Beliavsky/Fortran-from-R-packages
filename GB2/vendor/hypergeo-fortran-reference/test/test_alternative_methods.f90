! SPDX-License-Identifier: GPL-2.0-only
program test_alternative_methods
    use hypergeo_fortran, only : dp, crit_points, hypergeo_gosper, hypergeo_buhring, &
        hypergeo_ode_continue, f15_3_1
    implicit none
    complex(dp) :: a, b, c, z, v, cp(2)

    a = cmplx(1.21_dp, 0.0_dp, dp)
    b = cmplx(1.443_dp, 0.0_dp, dp)
    c = cmplx(1.88_dp, 0.0_dp, dp)

    cp = crit_points()
    v = hypergeo_gosper(a, b, c, cp(1), tol=1.0e-13_dp)
    call close(v, cmplx(0.57970091889194191008_dp, 0.84925906374150583119_dp, dp), &
        2.0e-12_dp, 'Gosper near critical point')

    z = cmplx(1.0_dp, 2.0_dp, dp)
    v = hypergeo_buhring(a, b, c, z, tol=1.0e-13_dp)
    call close(v, cmplx(0.03092203489529181523_dp, 0.55287678708672199067_dp, dp), &
        2.0e-11_dp, 'Buhring continuation')

    z = cmplx(2.13_dp, 0.68_dp, dp)
    v = hypergeo_ode_continue(a, b, c, z, rtol=1.0e-11_dp, atol=1.0e-13_dp)
    call close(v, cmplx(-0.71334498434004598167_dp, 0.59647479085505445853_dp, dp), &
        2.0e-9_dp, 'deSolve ODE continuation')

    z = cmplx(0.28_dp, -0.21_dp, dp)
    v = f15_3_1(a, b, c, z)
    call close(v, cmplx(1.26242625279896547055_dp, -0.33515733250598455101_dp, dp), &
        2.0e-12_dp, 'Euler integral')

    print *, 'test_alternative_methods: PASS'
contains
    subroutine close(x, y, tol, label)
        complex(dp), intent(in) :: x, y
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: label
        if (abs(x - y) > tol) then
            write(*, '(a,1x,es24.15)') 'FAIL ' // label // ' error:', abs(x - y)
            error stop 1
        end if
    end subroutine close
end program test_alternative_methods
