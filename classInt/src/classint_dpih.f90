! SPDX-License-Identifier: GPL-2.0-or-later
! SPDX-FileComment: Modern Fortran translation of R package classInt computational code.
module classint_dpih
    use classint_kinds, only: dp
    use classint_utils, only: mean_dp, sample_sd, quantile_r
    use classint_fft, only: fft_forward_radix2, fft_inverse_radix2, next_power_of_two
    implicit none
    private

    public :: dpih_bandwidth
    public :: dpih_bandwidth_direct_reference

contains

    function dpih_bandwidth(x, scale_method, level, gridsize, range_x, truncate) result(h)
        real(dp), intent(in) :: x(:) !! Finite sample used by the KernSmooth direct plug-in histogram bandwidth rule.
        character(len=*), intent(in) :: scale_method !! Scale estimator: "minim", "stdev", or "iqr".
        integer, intent(in) :: level !! Plug-in depth from zero through five; classInt/KernSmooth defaults to two.
        integer, intent(in) :: gridsize !! Number of equally spaced linear-binning grid points; must be at least two.
        real(dp), intent(in) :: range_x(2) !! Original-scale range over which binned functional estimates are formed.
        logical, intent(in) :: truncate !! If true, observations outside range_x are omitted during linear binning.
        real(dp) :: h

        h = dpih_bandwidth_impl(x, scale_method, level, gridsize, range_x, truncate, .true.)
    end function dpih_bandwidth

    function dpih_bandwidth_direct_reference(x, scale_method, level, gridsize, range_x, truncate) result(h)
        real(dp), intent(in) :: x(:) !! Finite sample used to validate FFT results against direct convolution.
        character(len=*), intent(in) :: scale_method !! Scale estimator: "minim", "stdev", or "iqr".
        integer, intent(in) :: level !! Plug-in depth from zero through five; classInt/KernSmooth defaults to two.
        integer, intent(in) :: gridsize !! Number of equally spaced grid points used by the direct reference calculation.
        real(dp), intent(in) :: range_x(2) !! Original-scale range over which binned functional estimates are formed.
        logical, intent(in) :: truncate !! If true, observations outside range_x are omitted during linear binning.
        real(dp) :: h

        h = dpih_bandwidth_impl(x, scale_method, level, gridsize, range_x, truncate, .false.)
    end function dpih_bandwidth_direct_reference

    function dpih_bandwidth_impl(x, scale_method, level, gridsize, range_x, truncate, use_fft) result(h)
        real(dp), intent(in) :: x(:) !! Finite sample used by the direct plug-in histogram bandwidth rule.
        character(len=*), intent(in) :: scale_method !! Scale estimator: "minim", "stdev", or "iqr".
        integer, intent(in) :: level !! Plug-in depth from zero through five; classInt/KernSmooth defaults to two.
        integer, intent(in) :: gridsize !! Number of equally spaced linear-binning grid points; must be at least two.
        real(dp), intent(in) :: range_x(2) !! Original-scale range over which binned functional estimates are formed.
        logical, intent(in) :: truncate !! If true, observations outside range_x are omitted during linear binning.
        logical, intent(in) :: use_fft !! Select KernSmooth FFT convolution; false retains the direct test oracle.
        real(dp) :: h
        real(dp), allocatable :: sx(:)
        real(dp), allocatable :: counts(:)
        real(dp) :: scale
        real(dp) :: mu
        real(dp) :: sa
        real(dp) :: sb
        real(dp) :: alpha
        real(dp) :: psi
        real(dp) :: root2pi
        integer :: n

        if (size(x) < 2) error stop "dpih_bandwidth: at least two observations are required"
        if (level < 0 .or. level > 5) error stop "dpih_bandwidth: level must be 0..5"
        if (gridsize < 2) error stop "dpih_bandwidth: gridsize must be at least two"
        if (range_x(2) <= range_x(1)) error stop "dpih_bandwidth: range must increase"
        n = size(x)
        select case (trim(adjustl(scale_method)))
        case ("stdev")
            scale = sample_sd(x)
        case ("iqr")
            scale = (quantile_r(x, 0.75_dp, 7) - quantile_r(x, 0.25_dp, 7)) / 1.349_dp
        case ("minim")
            scale = min(sample_sd(x), &
                        (quantile_r(x, 0.75_dp, 7) - quantile_r(x, 0.25_dp, 7)) / 1.349_dp)
        case default
            error stop "dpih_bandwidth: unknown scale estimator"
        end select
        if (scale <= 0.0_dp) error stop "dpih_bandwidth: scale estimate is zero"
        mu = mean_dp(x)
        sx = (x - mu) / scale
        sa = (range_x(1) - mu) / scale
        sb = (range_x(2) - mu) / scale
        call linear_bin(sx, sa, sb, gridsize, truncate, counts)
        root2pi = sqrt(2.0_dp / acos(-1.0_dp))

        select case (level)
        case (0)
            h = (24.0_dp * sqrt(acos(-1.0_dp)) / real(n, dp))**(1.0_dp / 3.0_dp)
        case (1)
            alpha = (2.0_dp / (3.0_dp * real(n, dp)))**(1.0_dp / 5.0_dp) * sqrt(2.0_dp)
            psi = bkfe_binned_dispatch(counts, 2, alpha, sa, sb, use_fft)
            h = (6.0_dp / (-psi * real(n, dp)))**(1.0_dp / 3.0_dp)
        case (2)
            alpha = (2.0_dp / (5.0_dp * real(n, dp)))**(1.0_dp / 7.0_dp) * sqrt(2.0_dp)
            psi = bkfe_binned_dispatch(counts, 4, alpha, sa, sb, use_fft)
            alpha = (root2pi / (psi * real(n, dp)))**(1.0_dp / 5.0_dp)
            psi = bkfe_binned_dispatch(counts, 2, alpha, sa, sb, use_fft)
            h = (6.0_dp / (-psi * real(n, dp)))**(1.0_dp / 3.0_dp)
        case (3)
            alpha = (2.0_dp / (7.0_dp * real(n, dp)))**(1.0_dp / 9.0_dp) * sqrt(2.0_dp)
            psi = bkfe_binned_dispatch(counts, 6, alpha, sa, sb, use_fft)
            alpha = (-3.0_dp * root2pi / (psi * real(n, dp)))**(1.0_dp / 7.0_dp)
            psi = bkfe_binned_dispatch(counts, 4, alpha, sa, sb, use_fft)
            alpha = (root2pi / (psi * real(n, dp)))**(1.0_dp / 5.0_dp)
            psi = bkfe_binned_dispatch(counts, 2, alpha, sa, sb, use_fft)
            h = (6.0_dp / (-psi * real(n, dp)))**(1.0_dp / 3.0_dp)
        case (4)
            alpha = (2.0_dp / (9.0_dp * real(n, dp)))**(1.0_dp / 11.0_dp) * sqrt(2.0_dp)
            psi = bkfe_binned_dispatch(counts, 8, alpha, sa, sb, use_fft)
            alpha = (15.0_dp * root2pi / (psi * real(n, dp)))**(1.0_dp / 9.0_dp)
            psi = bkfe_binned_dispatch(counts, 6, alpha, sa, sb, use_fft)
            alpha = (-3.0_dp * root2pi / (psi * real(n, dp)))**(1.0_dp / 7.0_dp)
            psi = bkfe_binned_dispatch(counts, 4, alpha, sa, sb, use_fft)
            alpha = (root2pi / (psi * real(n, dp)))**(1.0_dp / 5.0_dp)
            psi = bkfe_binned_dispatch(counts, 2, alpha, sa, sb, use_fft)
            h = (6.0_dp / (-psi * real(n, dp)))**(1.0_dp / 3.0_dp)
        case (5)
            alpha = (2.0_dp / (11.0_dp * real(n, dp)))**(1.0_dp / 13.0_dp) * sqrt(2.0_dp)
            psi = bkfe_binned_dispatch(counts, 10, alpha, sa, sb, use_fft)
            alpha = (-105.0_dp * root2pi / (psi * real(n, dp)))**(1.0_dp / 11.0_dp)
            psi = bkfe_binned_dispatch(counts, 8, alpha, sa, sb, use_fft)
            alpha = (15.0_dp * root2pi / (psi * real(n, dp)))**(1.0_dp / 9.0_dp)
            psi = bkfe_binned_dispatch(counts, 6, alpha, sa, sb, use_fft)
            alpha = (-3.0_dp * root2pi / (psi * real(n, dp)))**(1.0_dp / 7.0_dp)
            psi = bkfe_binned_dispatch(counts, 4, alpha, sa, sb, use_fft)
            alpha = (root2pi / (psi * real(n, dp)))**(1.0_dp / 5.0_dp)
            psi = bkfe_binned_dispatch(counts, 2, alpha, sa, sb, use_fft)
            h = (6.0_dp / (-psi * real(n, dp)))**(1.0_dp / 3.0_dp)
        end select
        h = scale * h
    end function dpih_bandwidth_impl

    subroutine linear_bin(x, a, b, m, truncate, counts)
        real(dp), intent(in) :: x(:) !! Standardized observations to distribute linearly between neighboring grid points.
        real(dp), intent(in) :: a !! Lower endpoint of the equally spaced binning grid.
        real(dp), intent(in) :: b !! Upper endpoint of the equally spaced binning grid; must exceed a.
        integer, intent(in) :: m !! Number of grid points; must be at least two.
        logical, intent(in) :: truncate !! If true, discard values outside [a,b]; otherwise clamp them to the nearest endpoint.
        real(dp), allocatable, intent(out) :: counts(:) !! Fractional grid counts summing to the retained observation count.
        real(dp) :: delta
        real(dp) :: t
        real(dp) :: frac
        real(dp) :: value
        integer :: i
        integer :: left

        allocate(counts(m), source=0.0_dp)
        delta = (b - a) / real(m - 1, dp)
        do i = 1, size(x)
            value = x(i)
            if (value < a) then
                if (truncate) cycle
                counts(1) = counts(1) + 1.0_dp
                cycle
            end if
            if (value > b) then
                if (truncate) cycle
                counts(m) = counts(m) + 1.0_dp
                cycle
            end if
            if (value >= b) then
                counts(m) = counts(m) + 1.0_dp
                cycle
            end if
            t = (value - a) / delta
            left = floor(t) + 1
            left = max(1, min(m - 1, left))
            frac = t - real(left - 1, dp)
            counts(left) = counts(left) + 1.0_dp - frac
            counts(left + 1) = counts(left + 1) + frac
        end do
    end subroutine linear_bin

    function bkfe_binned_dispatch(counts, drv, bandwidth, a, b, use_fft) result(psi)
        real(dp), intent(in) :: counts(:) !! Linear-binned sample counts on an equally spaced grid.
        integer, intent(in) :: drv !! Even Gaussian-kernel derivative order requested by the plug-in recursion.
        real(dp), intent(in) :: bandwidth !! Positive standardized kernel bandwidth.
        real(dp), intent(in) :: a !! Lower grid endpoint on the standardized scale.
        real(dp), intent(in) :: b !! Upper grid endpoint on the standardized scale.
        logical, intent(in) :: use_fft !! If true use KernSmooth FFT convolution; otherwise use the direct test oracle.
        real(dp) :: psi

        if (use_fft) then
            psi = bkfe_binned_fft(counts, drv, bandwidth, a, b)
        else
            psi = bkfe_binned_direct(counts, drv, bandwidth, a, b)
        end if
    end function bkfe_binned_dispatch

    function bkfe_binned_fft(counts, drv, bandwidth, a, b) result(psi)
        real(dp), intent(in) :: counts(:) !! Linear-binned sample counts on an equally spaced grid.
        integer, intent(in) :: drv !! Even Gaussian-kernel derivative order requested by the direct plug-in recursion.
        real(dp), intent(in) :: bandwidth !! Positive standardized kernel bandwidth.
        real(dp), intent(in) :: a !! Lower grid endpoint on the standardized scale.
        real(dp), intent(in) :: b !! Upper grid endpoint on the standardized scale.
        real(dp) :: psi
        real(dp), allocatable :: kernel(:)
        complex(dp), allocatable :: kernel_fft(:)
        complex(dp), allocatable :: counts_fft(:)
        real(dp) :: delta
        real(dp) :: arg
        real(dp) :: herm0
        real(dp) :: herm1
        real(dp) :: herm
        real(dp) :: n
        integer :: m
        integer :: lmax
        integer :: lag
        integer :: order
        integer :: transform_size

        if (bandwidth <= 0.0_dp) error stop "bkfe_binned_fft: bandwidth must be positive"
        m = size(counts)
        delta = (b - a) / real(m - 1, dp)
        lmax = min(floor(real(4 + drv, dp) * bandwidth / delta), m)
        allocate(kernel(0:lmax))
        do lag = 0, lmax
            arg = real(lag, dp) * delta / bandwidth
            herm0 = 1.0_dp
            herm1 = arg
            herm = 1.0_dp
            if (drv == 1) herm = herm1
            if (drv >= 2) then
                do order = 2, drv
                    herm = arg * herm1 - real(order - 1, dp) * herm0
                    herm0 = herm1
                    herm1 = herm
                end do
            end if
            kernel(lag) = herm * exp(-0.5_dp * arg * arg) / &
                          (sqrt(2.0_dp * acos(-1.0_dp)) * bandwidth**(drv + 1))
        end do

        n = sum(counts)
        if (n <= 0.0_dp) error stop "bkfe_binned_fft: no observations inside range"
        transform_size = next_power_of_two(m + lmax + 1)
        allocate(kernel_fft(transform_size), source=cmplx(0.0_dp, 0.0_dp, kind=dp))
        allocate(counts_fft(transform_size), source=cmplx(0.0_dp, 0.0_dp, kind=dp))
        kernel_fft(1:lmax + 1) = cmplx(kernel(0:lmax), 0.0_dp, kind=dp)
        do lag = 1, lmax
            kernel_fft(transform_size - lmax + lag) = cmplx(kernel(lmax - lag + 1), 0.0_dp, kind=dp)
        end do
        counts_fft(1:m) = cmplx(counts, 0.0_dp, kind=dp)

        call fft_forward_radix2(kernel_fft)
        call fft_forward_radix2(counts_fft)
        kernel_fft = kernel_fft * counts_fft
        call fft_inverse_radix2(kernel_fft)
        psi = dot_product(counts, real(kernel_fft(1:m), dp)) / (n * n)
    end function bkfe_binned_fft

    function bkfe_binned_direct(counts, drv, bandwidth, a, b) result(psi)
        real(dp), intent(in) :: counts(:) !! Linear-binned sample counts on an equally spaced grid.
        integer, intent(in) :: drv !! Even Gaussian-kernel derivative order requested by the direct plug-in recursion.
        real(dp), intent(in) :: bandwidth !! Positive standardized kernel bandwidth.
        real(dp), intent(in) :: a !! Lower grid endpoint on the standardized scale.
        real(dp), intent(in) :: b !! Upper grid endpoint on the standardized scale.
        real(dp) :: psi
        real(dp), allocatable :: kernel(:)
        real(dp) :: delta
        real(dp) :: arg
        real(dp) :: herm0
        real(dp) :: herm1
        real(dp) :: herm
        real(dp) :: conv
        real(dp) :: n
        integer :: m
        integer :: lmax
        integer :: lag
        integer :: i
        integer :: j
        integer :: order

        if (bandwidth <= 0.0_dp) error stop "bkfe_binned_direct: bandwidth must be positive"
        m = size(counts)
        delta = (b - a) / real(m - 1, dp)
        lmax = min(floor(real(4 + drv, dp) * bandwidth / delta), m)
        allocate(kernel(0:lmax))
        do lag = 0, lmax
            arg = real(lag, dp) * delta / bandwidth
            herm0 = 1.0_dp
            herm1 = arg
            herm = 1.0_dp
            if (drv == 1) herm = herm1
            if (drv >= 2) then
                do order = 2, drv
                    herm = arg * herm1 - real(order - 1, dp) * herm0
                    herm0 = herm1
                    herm1 = herm
                end do
            end if
            kernel(lag) = herm * exp(-0.5_dp * arg * arg) / &
                          (sqrt(2.0_dp * acos(-1.0_dp)) * bandwidth**(drv + 1))
        end do
        n = sum(counts)
        if (n <= 0.0_dp) error stop "bkfe_binned_direct: no observations inside range"
        psi = 0.0_dp
        do i = 1, m
            conv = 0.0_dp
            do j = max(1, i - lmax), min(m, i + lmax)
                lag = abs(i - j)
                conv = conv + kernel(lag) * counts(j)
            end do
            psi = psi + counts(i) * conv
        end do
        psi = psi / (n * n)
    end function bkfe_binned_direct

end module classint_dpih
