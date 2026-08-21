! SPDX-License-Identifier: GPL-2.0-only
program hypergeo_demo
    use hypergeo_fortran, only : dp, hypergeo, hypergeo_info
    implicit none
    complex(dp) :: value
    type(hypergeo_info) :: info

    value = hypergeo(cmplx(1.21_dp, 0.0_dp, dp), cmplx(1.443_dp, 0.0_dp, dp), &
        cmplx(1.88_dp, 0.0_dp, dp), cmplx(2.13_dp, 0.68_dp, dp), info=info)

    write(*, '(a,2(1x,es22.14))') '2F1 =', real(value), aimag(value)
    write(*, '(a,1x,a)') 'method =', trim(info%method)
end program hypergeo_demo
