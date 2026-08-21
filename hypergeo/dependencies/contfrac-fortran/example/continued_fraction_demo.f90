! SPDX-License-Identifier: GPL-2.0-only
program continued_fraction_demo
    use contfrac, only : dp, cf, gcf, convergents, as_cf
    implicit none

    real(dp), allocatable :: numerators(:), denominators(:), terms(:), a(:), b(:)
    real(dp) :: phi
    complex(dp), allocatable :: ca(:), cb(:)
    complex(dp) :: z
    integer :: i

    phi = 0.5_dp * (1.0_dp + sqrt(5.0_dp))
    a = [(1.0_dp, i = 1, 100)]
    print '(a,f24.16)', "phi from CF       = ", cf(a)
    print '(a,f24.16)', "intrinsic phi     = ", phi

    terms = as_cf(acos(-1.0_dp), 10)
    print '(a,10(f0.0,1x))', "pi coefficients   = ", terms

    call convergents([3.0_dp, 7.0_dp, 15.0_dp, 1.0_dp, 292.0_dp], numerators, denominators)
    print '(a,f24.16)', "pi convergent     = ", numerators(size(numerators)) / denominators(size(denominators))

    deallocate(a)
    allocate(a(30), b(30))
    do i = 1, 30
        a(i) = real(2 * i, dp)
        b(i) = real(2 * i + 1, dp)
    end do
    print '(a,f24.16)', "Euler identity    = ", gcf(a, b, b0=1.0_dp)

    z = cmplx(1.0_dp, 1.0_dp, kind=dp)
    allocate(ca(100), cb(100))
    ca(1) = z
    ca(2:) = -z * z
    do i = 1, 100
        cb(i) = cmplx(real(2 * i - 1, dp), 0.0_dp, kind=dp)
    end do
    print '(a,2f24.16)', "tan(1+i) from GCF = ", gcf(ca, cb)

end program continued_fraction_demo
