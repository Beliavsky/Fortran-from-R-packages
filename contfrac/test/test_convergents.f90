! SPDX-License-Identifier: GPL-2.0-only
program test_convergents
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use contfrac, only : dp, convergents, gconvergents, gcf
    implicit none

    real(dp), allocatable :: a_num(:), a_den(:), a(:), b(:)
    complex(dp), allocatable :: c_num(:), c_den(:), ca(:), cb(:)
    real(dp) :: pi, value
    complex(dp) :: z, cvalue
    integer :: i

    pi = acos(-1.0_dp)
    call convergents([3.0_dp, 7.0_dp, 15.0_dp, 1.0_dp, 292.0_dp], a_num, a_den)
    if (size(a_num) /= 5 .or. size(a_den) /= 5) error stop "wrong convergent size"
    call assert_close(a_num(1) / a_den(1), 3.0_dp, 0.0_dp, "pi conv 1")
    call assert_close(a_num(2) / a_den(2), 22.0_dp / 7.0_dp, 0.0_dp, "pi conv 2")
    call assert_close(a_num(5) / a_den(5), 103993.0_dp / 33102.0_dp, 0.0_dp, "pi conv 5")
    if (abs(a_num(5) / a_den(5) - pi) > abs(a_num(2) / a_den(2) - pi)) then
        error stop "pi convergents did not improve"
    end if

    allocate(a(30), b(30))
    do i = 1, 30
        a(i) = real(2 * i, dp)
        b(i) = real(2 * i + 1, dp)
    end do
    call gconvergents(a, b, 1.0_dp, a_num, a_den)
    value = a_num(size(a_num)) / a_den(size(a_den))
    call assert_close(value, 1.0_dp / (exp(0.5_dp) - 1.0_dp), 2.0e-12_dp, "Euler identity convergent")
    call assert_close(gcf(a, b, b0=1.0_dp), value, 2.0e-12_dp, "GCF/convergent agreement")

    z = cmplx(0.3_dp, -0.4_dp, kind=dp)
    allocate(ca(24), cb(24))
    ca(1) = z
    ca(2:) = -z * z
    do i = 1, size(cb)
        cb(i) = cmplx(real(2 * i - 1, dp), 0.0_dp, kind=dp)
    end do
    call gconvergents(ca, cb, cmplx(0.0_dp, 0.0_dp, kind=dp), c_num, c_den)
    cvalue = c_num(size(c_num)) / c_den(size(c_den))
    call assert_cclose(cvalue, tan(z), 1.0e-13_dp, "complex convergents")

    print '(a)', "test_convergents: PASS"

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

end program test_convergents
