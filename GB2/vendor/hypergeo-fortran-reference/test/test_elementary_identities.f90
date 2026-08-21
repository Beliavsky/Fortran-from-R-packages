! SPDX-License-Identifier: GPL-2.0-only
program test_elementary_identities
    use hypergeo_fortran, only : dp, hypergeo
    implicit none
    complex(dp) :: z(6), x, lhs, rhs
    integer :: i

    z = [ &
        cmplx(0.28_dp, -0.21_dp, dp), cmplx(-0.79_dp, -0.40_dp, dp), &
        cmplx(0.56_dp, 0.05_dp, dp), cmplx(2.13_dp, 0.68_dp, dp), &
        cmplx(-0.43_dp, -1.47_dp, dp), cmplx(1.23_dp, 0.48_dp, dp)]

    do i = 1, size(z)
        x = z(i)
        lhs = hypergeo(cmplx(1.0_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            cmplx(2.0_dp, 0.0_dp, dp), x)
        rhs = -log(1.0_dp - x) / x
        call close(lhs, rhs, 2.0e-10_dp, 'A&S 15.1.3')

        lhs = hypergeo(cmplx(0.5_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            cmplx(1.5_dp, 0.0_dp, dp), x * x)
        rhs = 0.5_dp * log((1.0_dp + x) / (1.0_dp - x)) / x
        call close(lhs, rhs, 2.0e-10_dp, 'A&S 15.1.4')

        lhs = hypergeo(cmplx(0.5_dp, 0.0_dp, dp), cmplx(1.0_dp, 0.0_dp, dp), &
            cmplx(1.5_dp, 0.0_dp, dp), -x * x)
        rhs = atan(x) / x
        call close(lhs, rhs, 2.0e-10_dp, 'A&S 15.1.5')

        lhs = hypergeo(cmplx(0.5_dp, 0.0_dp, dp), cmplx(0.5_dp, 0.0_dp, dp), &
            cmplx(1.5_dp, 0.0_dp, dp), x * x)
        rhs = asin(x) / x
        call close(lhs, rhs, 2.0e-10_dp, 'A&S 15.1.6b')
    end do

    print *, 'test_elementary_identities: PASS'
contains
    subroutine close(a, b, tol, label)
        complex(dp), intent(in) :: a, b
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: label
        if (abs(a - b) > tol) then
            write(*, '(a,2(1x,es24.15))') 'FAIL ' // label // ' error:', abs(a - b), tol
            error stop 1
        end if
    end subroutine close
end program test_elementary_identities
