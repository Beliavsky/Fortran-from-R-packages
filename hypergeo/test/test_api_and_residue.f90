! SPDX-License-Identifier: GPL-2.0-only
program test_api_and_residue
    use hypergeo_fortran, only : dp, hypergeo, hypergeo_residue_general, to_real, to_complex, &
        hypergeo_shanks, genhypergeo_shanks
    implicit none
    complex(dp) :: a, b, c, z, v, r, zz(2), back(2), s1, s2
    real(dp) :: rv(4)

    a = cmplx(0.7_dp, 0.0_dp, dp)
    b = cmplx(1.1_dp, 0.0_dp, dp)
    c = cmplx(2.3_dp, 0.0_dp, dp)
    z = cmplx(0.2_dp, 0.1_dp, dp)
    v = hypergeo(a, b, c, z)
    r = hypergeo_residue_general(a, b, c, z, r=0.08_dp)
    call close(r, v, 2.0e-8_dp, 'Cauchy/residue evaluation')

    zz = [cmplx(1.0_dp, -2.0_dp, dp), cmplx(3.5_dp, 4.25_dp, dp)]
    rv = to_real(zz)
    back = to_complex(rv)
    call check(maxval(abs(back - zz)) < 1.0e-15_dp, 'complex/real packing')

    s1 = hypergeo_shanks(a, b, c, cmplx(0.2_dp, 0.0_dp, dp), maxiter=10)
    s2 = genhypergeo_shanks([a, b], [c], cmplx(0.2_dp, 0.0_dp, dp), maxiter=10)
    call close(s1, s2, 1.0e-15_dp, 'hypergeo_shanks wrapper')

    print *, 'test_api_and_residue: PASS'
contains
    subroutine check(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            write(*, '(a)') 'FAIL: ' // label
            error stop 1
        end if
    end subroutine check

    subroutine close(x, y, tol, label)
        complex(dp), intent(in) :: x, y
        real(dp), intent(in) :: tol
        character(len=*), intent(in) :: label
        call check(abs(x - y) <= tol, label)
    end subroutine close
end program test_api_and_residue
