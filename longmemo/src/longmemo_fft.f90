! Part of the modern Fortran translation of longmemo 1.1-4.
! This file was created or modified for the Fortran project on 2026-07-23.
! Original longmemo authors retain copyright; see ORIGINAL_PACKAGE.txt.
! SPDX-License-Identifier: GPL-2.0-or-later

module longmemo_fft
    use longmemo_kinds, only : dp, pi, two_pi
    implicit none
    private

    public :: fft_forward, fft_inverse, fft_inverse_raw

contains

    subroutine fft_forward(x, y)
        complex(dp), intent(in) :: x(:)
        complex(dp), allocatable, intent(out) :: y(:)

        allocate(y(size(x)))
        y = x
        if (size(x) <= 1) return

        if (is_power_of_two(size(x))) then
            call radix2_transform(y, inverse=.false., normalize=.false.)
        else
            call bluestein_forward(x, y)
        end if
    end subroutine fft_forward


    subroutine fft_inverse_raw(x, y)
        complex(dp), intent(in) :: x(:)
        complex(dp), allocatable, intent(out) :: y(:)
        complex(dp), allocatable :: work(:)

        allocate(work(size(x)))
        work = conjg(x)
        call fft_forward(work, y)
        y = conjg(y)
    end subroutine fft_inverse_raw


    subroutine fft_inverse(x, y)
        complex(dp), intent(in) :: x(:)
        complex(dp), allocatable, intent(out) :: y(:)

        call fft_inverse_raw(x, y)
        if (size(x) > 0) y = y/real(size(x), dp)
    end subroutine fft_inverse


    subroutine bluestein_forward(x, y)
        complex(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: y(:)
        complex(dp), allocatable :: a(:), b(:)
        complex(dp) :: chirp
        real(dp) :: angle
        integer :: n, m, j

        n = size(x)
        m = next_power_of_two(2*n - 1)
        allocate(a(m), b(m))
        a = (0.0_dp, 0.0_dp)
        b = (0.0_dp, 0.0_dp)

        do j = 0, n - 1
            angle = pi*real(j, dp)*real(j, dp)/real(n, dp)
            chirp = cmplx(cos(angle), -sin(angle), dp)
            a(j + 1) = x(j + 1)*chirp

            chirp = cmplx(cos(angle), sin(angle), dp)
            b(j + 1) = chirp
            if (j > 0) b(m - j + 1) = chirp
        end do

        call radix2_transform(a, inverse=.false., normalize=.false.)
        call radix2_transform(b, inverse=.false., normalize=.false.)
        a = a*b
        call radix2_transform(a, inverse=.true., normalize=.true.)

        do j = 0, n - 1
            angle = pi*real(j, dp)*real(j, dp)/real(n, dp)
            chirp = cmplx(cos(angle), -sin(angle), dp)
            y(j + 1) = a(j + 1)*chirp
        end do
    end subroutine bluestein_forward


    subroutine radix2_transform(x, inverse, normalize)
        complex(dp), intent(inout) :: x(:)
        logical, intent(in) :: inverse, normalize
        complex(dp) :: temp, w, wm
        real(dp) :: angle
        integer :: n, i, j, bit, len, half, k

        n = size(x)
        if (.not. is_power_of_two(n)) error stop "radix2_transform: size is not a power of two"

        j = 0
        do i = 1, n - 1
            bit = n/2
            do while (iand(j, bit) /= 0)
                j = ieor(j, bit)
                bit = bit/2
            end do
            j = ieor(j, bit)
            if (i < j) then
                temp = x(i + 1)
                x(i + 1) = x(j + 1)
                x(j + 1) = temp
            end if
        end do

        len = 2
        do while (len <= n)
            angle = merge(two_pi/real(len, dp), -two_pi/real(len, dp), inverse)
            wm = cmplx(cos(angle), sin(angle), dp)
            half = len/2
            do i = 1, n, len
                w = (1.0_dp, 0.0_dp)
                do k = 0, half - 1
                    temp = w*x(i + k + half)
                    x(i + k + half) = x(i + k) - temp
                    x(i + k) = x(i + k) + temp
                    w = w*wm
                end do
            end do
            len = 2*len
        end do

        if (inverse .and. normalize) x = x/real(n, dp)
    end subroutine radix2_transform


    pure logical function is_power_of_two(n)
        integer, intent(in) :: n

        is_power_of_two = n > 0 .and. iand(n, n - 1) == 0
    end function is_power_of_two


    pure integer function next_power_of_two(n) result(p)
        integer, intent(in) :: n

        p = 1
        do while (p < n)
            p = 2*p
        end do
    end function next_power_of_two

end module longmemo_fft
