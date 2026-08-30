! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_fft
    use classint_kinds, only: dp
    implicit none
    private

    public :: fft_forward_radix2
    public :: fft_inverse_radix2
    public :: next_power_of_two

contains

    pure subroutine fft_forward_radix2(values)
        !! Replaces a power-of-two complex sequence with its unnormalized forward DFT.
        complex(dp), intent(inout) :: values(:) !! Sequence transformed in place; its size must be a power of two.
        integer :: bit, block_length, i, j, n, offset
        complex(dp) :: temp, twiddle, twiddle_step
        real(dp) :: pi

        n = size(values)
        if (.not. is_power_of_two(n)) error stop "fft_forward_radix2: sequence length must be a power of two"
        if (n == 1) return
        pi = acos(-1.0_dp)

        j = 1
        do i = 2, n
            bit = n / 2
            do while (j > bit)
                j = j - bit
                bit = bit / 2
                if (bit == 0) exit
            end do
            j = j + bit
            if (i < j) then
                temp = values(i)
                values(i) = values(j)
                values(j) = temp
            end if
        end do

        block_length = 2
        do while (block_length <= n)
            twiddle_step = exp(cmplx(0.0_dp, -2.0_dp * pi / real(block_length, dp), kind=dp))
            do i = 1, n, block_length
                twiddle = cmplx(1.0_dp, 0.0_dp, kind=dp)
                do offset = 0, block_length / 2 - 1
                    temp = twiddle * values(i + offset + block_length / 2)
                    values(i + offset + block_length / 2) = values(i + offset) - temp
                    values(i + offset) = values(i + offset) + temp
                    twiddle = twiddle * twiddle_step
                end do
            end do
            block_length = block_length * 2
        end do
    end subroutine fft_forward_radix2

    pure subroutine fft_inverse_radix2(values)
        !! Replaces a power-of-two spectrum with its normalized inverse DFT.
        complex(dp), intent(inout) :: values(:) !! Spectrum transformed in place; its size must be a power of two.
        integer :: n

        n = size(values)
        if (.not. is_power_of_two(n)) error stop "fft_inverse_radix2: sequence length must be a power of two"
        values = conjg(values)
        call fft_forward_radix2(values)
        values = conjg(values) / real(n, dp)
    end subroutine fft_inverse_radix2

    pure elemental function next_power_of_two(n) result(power)
        !! Returns the smallest power of two greater than or equal to `n`.
        integer, intent(in) :: n !! Positive lower bound for the returned power of two.
        integer :: power

        if (n < 1) then
            power = 1
            return
        end if
        power = 1
        do while (power < n)
            if (power > ishft(huge(power), -1)) error stop "next_power_of_two: requested transform is too large"
            power = power * 2
        end do
    end function next_power_of_two

    pure elemental function is_power_of_two(n) result(ok)
        !! Reports whether `n` is a positive integral power of two.
        integer, intent(in) :: n !! Integer tested for positive power-of-two structure.
        logical :: ok
        integer :: reduced

        if (n < 1) then
            ok = .false.
            return
        end if
        reduced = n
        do while (modulo(reduced, 2) == 0)
            reduced = reduced / 2
        end do
        ok = reduced == 1
    end function is_power_of_two

end module classint_fft
