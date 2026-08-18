! SPDX-License-Identifier: GPL-2.0-only
program test_cf
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use contfrac, only : dp, cf, gcf, as_cf, contfrac_info
    implicit none

    real(dp), parameter :: tol = 1.0e-10_dp
    real(dp), allocatable :: a(:), terms(:)
    real(dp) :: value, pi, gamma_e
    complex(dp), allocatable :: ca(:), cb(:)
    complex(dp) :: z, cvalue
    type(contfrac_info) :: info
    integer :: i

    pi = acos(-1.0_dp)
    gamma_e = 0.577215664901532860606512090082402431_dp

    allocate(a(61))
    a(1) = 3.0_dp
    do i = 2, 61, 2
        a(i) = 3.0_dp
        a(i + 1) = 6.0_dp
    end do
    value = cf(a)
    call assert_close(value, sqrt(11.0_dp), tol, "sqrt(11)")
    deallocate(a)

    allocate(a(81))
    a(1) = 8.0_dp
    do i = 0, 9
        a(2 + 8 * i:9 + 8 * i) = [2.0_dp, 2.0_dp, 1.0_dp, 7.0_dp, &
                                      1.0_dp, 2.0_dp, 2.0_dp, 16.0_dp]
    end do
    value = cf(a)
    call assert_close(value, sqrt(71.0_dp), tol, "sqrt(71)")
    deallocate(a)

    allocate(a(91))
    a(1) = 2.0_dp
    do i = 1, 30
        a(3 * i - 1) = 1.0_dp
        a(3 * i) = real(2 * i, dp)
        a(3 * i + 1) = 1.0_dp
    end do
    value = cf(a)
    call assert_close(value, exp(1.0_dp), tol, "exp(1)")
    deallocate(a)

    a = [0.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 4.0_dp, &
         3.0_dp, 13.0_dp, 5.0_dp, 1.0_dp, 1.0_dp, 8.0_dp, 1.0_dp, 2.0_dp, &
         4.0_dp, 1.0_dp, 1.0_dp, 40.0_dp, 1.0_dp, 11.0_dp, 3.0_dp, 7.0_dp, &
         1.0_dp, 7.0_dp, 1.0_dp, 1.0_dp, 5.0_dp, 1.0_dp, 49.0_dp, 4.0_dp, &
         1.0_dp, 65.0_dp, 1.0_dp, 4.0_dp, 7.0_dp, 11.0_dp, 1.0_dp, 399.0_dp, &
         2.0_dp, 1.0_dp, 3.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 5.0_dp, &
         3.0_dp, 2.0_dp, 1.0_dp, 10.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, 1.0_dp, &
         2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp, 1.0_dp, 4.0_dp, 1.0_dp, 1.0_dp, &
         2.0_dp, 5.0_dp, 1.0_dp, 3.0_dp, 6.0_dp, 2.0_dp, 1.0_dp, 2.0_dp, &
         1.0_dp, 1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, 3.0_dp, 16.0_dp, 8.0_dp, &
         1.0_dp, 1.0_dp, 2.0_dp, 16.0_dp]
    value = cf(a)
    call assert_close(value, gamma_e, tol, "Euler gamma")

    z = cmplx(1.0_dp, 1.0_dp, kind=dp)
    allocate(ca(100), cb(100))
    ca(1) = z
    ca(2:) = -z * z
    do i = 1, 100
        cb(i) = cmplx(real(2 * i - 1, dp), 0.0_dp, kind=dp)
    end do
    cvalue = gcf(ca, cb)
    call assert_cclose(cvalue, tan(z), tol, "tan(1+i)")

    terms = as_cf(sqrt(2.0_dp), 12)
    do i = 2, size(terms)
        call assert_close(terms(i), 2.0_dp, 0.0_dp, "sqrt(2) expansion")
    end do

    value = gcf([2.0_dp], [3.0_dp], b0=1.0_dp, finite=.true., info=info)
    call assert_close(value, 5.0_dp / 3.0_dp, 5.0e-15_dp, "finite one-level GCF")
    if (info%iterations /= 1) error stop "unexpected iteration count"

    print '(a)', "test_cf: PASS"

contains

    subroutine assert_close(actual, expected, atol, label)
        real(dp), intent(in) :: actual, expected, atol
        character(len=*), intent(in) :: label
        if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > atol) then
            print '(a,2es25.16)', trim(label) // ": ", actual, expected
            error stop 1
        end if
    end subroutine assert_close

    subroutine assert_cclose(actual, expected, atol, label)
        complex(dp), intent(in) :: actual, expected
        real(dp), intent(in) :: atol
        character(len=*), intent(in) :: label
        if (.not. ieee_is_finite(real(actual, dp)) .or. .not. ieee_is_finite(aimag(actual)) .or. &
            abs(actual - expected) > atol) then
            print '(a,4es25.16)' , trim(label) // ": ", real(actual), aimag(actual), &
                real(expected), aimag(expected)
            error stop 1
        end if
    end subroutine assert_cclose

end program test_cf
