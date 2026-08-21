! SPDX-License-Identifier: GPL-2.0-only
program test_special_and_integer
    use hypergeo_fortran, only : dp, pi, complex_gamma, genhypergeo_series, hypergeo, f15_3_11
    implicit none
    complex(dp) :: v, u(2), l(1)

    v = complex_gamma(cmplx(1.0_dp, 1.0_dp, dp))
    call close(v, cmplx(0.49801566811835604271_dp, -0.15494982830181068512_dp, dp), &
        2.0e-13_dp, 'complex gamma')

    u = [cmplx(0.2_dp, 0.0_dp, dp), cmplx(1.3_dp, 0.0_dp, dp)]
    l = [cmplx(2.7_dp, 0.0_dp, dp)]
    v = genhypergeo_series(u, l, cmplx(0.3_dp, -0.2_dp, dp))
    call close(v, cmplx(1.03032894282768039163_dp, -0.02455481878563063736_dp, dp), &
        2.0e-13_dp, 'generalized series')

    v = hypergeo(cmplx(pi, 0.0_dp, dp), cmplx(-4.0_dp, 0.0_dp, dp), &
        cmplx(2.2_dp, 0.0_dp, dp), cmplx(1.0_dp, 5.0_dp, dp))
    call close(v, cmplx(1670.8287595795885335_dp, -204.81995157365381258_dp, dp), &
        2.0e-10_dp, 'terminating polynomial')

    v = f15_3_11(cmplx(1.0_dp, 0.0_dp, dp), cmplx(3.0_dp, 0.0_dp, dp), &
        cmplx(2.0_dp, 0.0_dp, dp), cmplx(0.9_dp, 0.01_dp, dp))
    call close(v, cmplx(2.04925816767572859288_dp, 0.03163091033158608114_dp, dp), &
        1.0e-11_dp, 'A&S 15.3.11')

    print *, 'test_special_and_integer: PASS'
contains
    subroutine close(a, b, tol, label)
        complex(dp), intent(in) :: a, b
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: label
        if (abs(a - b) > tol) then
            write(*, '(a,1x,es24.15)') 'FAIL ' // label // ' error:', abs(a - b)
            error stop 1
        end if
    end subroutine close
end program test_special_and_integer
