! SPDX-License-Identifier: GPL-2.0-or-later
module ltsa_fft
    use ltsa_kinds, only : dp, two_pi
    implicit none
    private

    public :: fft_inplace, next_power_of_two

contains

    integer function next_power_of_two(n) result(p)
        integer, intent(in) :: n
        p = 1
        do while (p < max(1,n))
            p = 2*p
        end do
    end function next_power_of_two

    subroutine fft_inplace(x, inverse)
        complex(dp), intent(inout) :: x(:)
        logical, intent(in), optional :: inverse
        logical :: inv
        integer :: n, i, j, m, len, half, k
        complex(dp) :: temp, w, wlen
        real(dp) :: angle
        inv = .false.
        if (present(inverse)) inv = inverse
        n = size(x)
        j = 1
        do i = 1, n
            if (i < j) then
                temp = x(i)
                x(i) = x(j)
                x(j) = temp
            end if
            m = n/2
            do while (m >= 1 .and. j > m)
                j = j-m
                m = m/2
            end do
            j = j+m
        end do
        len = 2
        do while (len <= n)
            angle = two_pi/real(len,dp)
            if (.not. inv) angle = -angle
            wlen = cmplx(cos(angle), sin(angle), dp)
            half = len/2
            do i = 1, n, len
                w = cmplx(1.0_dp,0.0_dp,dp)
                do k = 0, half-1
                    temp = w*x(i+k+half)
                    x(i+k+half) = x(i+k)-temp
                    x(i+k) = x(i+k)+temp
                    w = w*wlen
                end do
            end do
            len = 2*len
        end do
        if (inv) x = x/real(n,dp)
    end subroutine fft_inplace

end module ltsa_fft
