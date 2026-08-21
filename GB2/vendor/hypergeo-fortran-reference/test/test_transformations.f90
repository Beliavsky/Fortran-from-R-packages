! SPDX-License-Identifier: GPL-2.0-only
program test_transformations
    use hypergeo_fortran, only : dp, hypergeo, f15_3_3, f15_3_4, f15_3_5, &
        f15_3_6, f15_3_7, f15_3_8, f15_3_9
    implicit none
    complex(dp) :: a, b, c, z, ref

    a = cmplx(0.73_dp, 0.0_dp, dp)
    b = cmplx(1.17_dp, 0.0_dp, dp)
    c = cmplx(2.31_dp, 0.0_dp, dp)
    z = cmplx(-0.35_dp, 0.22_dp, dp)
    ref = hypergeo(a, b, c, z)
    call close(f15_3_3(a, b, c, z), ref, 2.0e-11_dp, '15.3.3')
    call close(f15_3_4(a, b, c, z), ref, 2.0e-11_dp, '15.3.4')
    call close(f15_3_5(a, b, c, z), ref, 2.0e-11_dp, '15.3.5')

    z = cmplx(0.72_dp, 0.18_dp, dp)
    ref = hypergeo(a, b, c, z)
    call close(f15_3_6(a, b, c, z), ref, 5.0e-10_dp, '15.3.6')

    z = cmplx(2.2_dp, 0.7_dp, dp)
    ref = hypergeo(a, b, c, z)
    call close(f15_3_7(a, b, c, z), ref, 5.0e-10_dp, '15.3.7')
    call close(f15_3_8(a, b, c, z), ref, 5.0e-10_dp, '15.3.8')
    call close(f15_3_9(a, b, c, z), ref, 5.0e-10_dp, '15.3.9')

    print *, 'test_transformations: PASS'
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
end program test_transformations
