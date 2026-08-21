! FFT support for pmultinom-fortran.
! This is a clean Fortran implementation used in place of the R package's
! fftw dependency. The surrounding translation remains AGPL-3.0-only.
module pmultinom_fft
    use pmultinom_kinds, only : dp
    implicit none
    private

    integer, parameter :: hp_try = selected_real_kind(18, 300)
    integer, parameter :: hp = merge(hp_try, dp, hp_try > 0)
    real(hp), parameter :: pi = acos(-1.0_hp)

    public :: convolve_prefix

contains

    pure integer function next_power_of_two(n) result(p)
        integer, intent(in) :: n

        p = 1
        do while (p < n)
            p = 2 * p
        end do
    end function next_power_of_two

    subroutine fft_inplace(a, inverse)
        complex(hp), intent(inout) :: a(0:)
        logical, intent(in) :: inverse

        integer :: n, i, j, bit, len, half, k
        real(hp) :: angle
        complex(hp) :: w, wlen, u, v, tmp

        n = size(a)
        if (n <= 1) return

        j = 0
        do i = 1, n - 1
            bit = n / 2
            do while (iand(j, bit) /= 0)
                j = ieor(j, bit)
                bit = bit / 2
            end do
            j = ieor(j, bit)
            if (i < j) then
                tmp = a(i)
                a(i) = a(j)
                a(j) = tmp
            end if
        end do

        len = 2
        do while (len <= n)
            angle = 2.0_hp * pi / real(len, hp)
            if (.not. inverse) angle = -angle
            wlen = cmplx(cos(angle), sin(angle), kind=hp)
            half = len / 2
            do i = 0, n - 1, len
                w = cmplx(1.0_hp, 0.0_hp, kind=hp)
                do k = 0, half - 1
                    u = a(i + k)
                    v = a(i + k + half) * w
                    a(i + k) = u + v
                    a(i + k + half) = u - v
                    w = w * wlen
                end do
            end do
            len = 2 * len
        end do

        if (inverse) a = a / real(n, hp)
    end subroutine fft_inplace

    subroutine convolve_prefix(x, y, nmax, z)
        real(dp), intent(in) :: x(0:), y(0:)
        integer, intent(in) :: nmax
        real(dp), intent(out) :: z(0:nmax)

        integer :: nfft, nx, ny, i, ncopy
        complex(hp), allocatable :: fx(:), fy(:)

        nx = min(size(x), nmax + 1)
        ny = min(size(y), nmax + 1)
        nfft = next_power_of_two(nx + ny - 1)
        allocate(fx(0:nfft-1), fy(0:nfft-1))
        fx = cmplx(0.0_hp, 0.0_hp, kind=hp)
        fy = cmplx(0.0_hp, 0.0_hp, kind=hp)
        do i = 0, nx - 1
            fx(i) = cmplx(real(x(i), hp), 0.0_hp, kind=hp)
        end do
        do i = 0, ny - 1
            fy(i) = cmplx(real(y(i), hp), 0.0_hp, kind=hp)
        end do

        call fft_inplace(fx, .false.)
        call fft_inplace(fy, .false.)
        fx = fx * fy
        call fft_inplace(fx, .true.)

        z = 0.0_dp
        ncopy = min(nmax, nx + ny - 2)
        do i = 0, ncopy
            z(i) = real(fx(i), dp)
        end do

        ! FFT roundoff can produce tiny negative values for an exact
        ! convolution of nonnegative sequences.
        where (z < 0.0_dp .and. abs(z) <= 256.0_dp * epsilon(1.0_dp)) z = 0.0_dp
    end subroutine convolve_prefix

end module pmultinom_fft
