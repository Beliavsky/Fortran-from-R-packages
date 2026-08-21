! SPDX-License-Identifier: GPL-2.0-only
program test_upstream_refs
    use hypergeo_fortran, only : dp, pi, hypergeo, genhypergeo_contfrac_single
    implicit none
    complex(dp) :: z(6), got(6), ref(6), v
    complex(dp) :: u(1), l(3)

    z = [ &
        cmplx(0.28_dp, -0.21_dp, dp), cmplx(-0.79_dp, -0.40_dp, dp), &
        cmplx(0.56_dp, 0.05_dp, dp), cmplx(2.13_dp, 0.68_dp, dp), &
        cmplx(-0.43_dp, -1.47_dp, dp), cmplx(1.23_dp, 0.48_dp, dp)]
    ref = [ &
        cmplx(1.26242625279896547050_dp, -0.33515733250598455101_dp, dp), &
        cmplx(0.55218026726586346626_dp, -0.11826012518395685586_dp, dp), &
        cmplx(2.08656513655255552950_dp, 0.21074089910422408972_dp, dp), &
        cmplx(-0.71334498434004598167_dp, 0.59647479085505445853_dp, dp), &
        cmplx(0.36911468617705947291_dp, -0.35906488952504313903_dp, dp), &
        cmplx(-0.44806924103752606401_dp, 1.91140611055324833040_dp, dp)]
    got = hypergeo(cmplx(1.21_dp, 0.0_dp, dp), cmplx(1.443_dp, 0.0_dp, dp), &
        cmplx(1.88_dp, 0.0_dp, dp), z)
    call check(maxval(abs(got - ref)) < 1.0e-9_dp, 'Maple six-point regression')

    u = [cmplx(0.2_dp, 0.0_dp, dp)]
    l = [cmplx(9.9_dp, 0.0_dp, dp), cmplx(2.7_dp, 0.0_dp, dp), cmplx(8.7_dp, 0.0_dp, dp)]
    v = genhypergeo_contfrac_single(u, l, cmplx(1.0_dp, 10.0_dp, dp))
    call close(v, cmplx(1.0007289707983569879_dp, 0.0086250714217251837317_dp, dp), &
        1.0e-11_dp, 'generalized continued fraction')

    v = hypergeo(cmplx(pi, 0.0_dp, dp), cmplx(pi / 2.0_dp, 0.0_dp, dp), &
        cmplx(3.0_dp * pi / 2.0_dp - 4.0_dp, 0.0_dp, dp), cmplx(0.1_dp, 0.2_dp, dp))
    call close(v, cmplx(0.53745229690249593045_dp, 1.8917456473240515664_dp, dp), &
        1.0e-10_dp, 'cover1 negative m')

    v = hypergeo(cmplx(pi, 0.0_dp, dp), cmplx(pi / 2.0_dp, 0.0_dp, dp), &
        cmplx(3.0_dp * pi / 2.0_dp + 4.0_dp, 0.0_dp, dp), cmplx(10.1_dp, 0.2_dp, dp))
    call close(v, cmplx(-0.29639970263878733845_dp, -0.34765230143995441172_dp, dp), &
        1.0e-10_dp, 'cover1 positive m')

    print *, 'test_upstream_refs: PASS'
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
end program test_upstream_refs
