! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_ad
    use goftest_kinds, only : dp
    use goftest_utils, only : sort_real, clamp01
    implicit none
    private

    public :: ad_cdf
    public :: ad_quantile
    public :: ad_statistic
    public :: ad_test_uniform
    public :: ad_inf_fast
    public :: ad_inf_exact

contains

    pure real(dp) function ad_inf_fast(z) result(p)
        real(dp), intent(in) :: z

        if (z <= 0.0_dp) then
            p = 0.0_dp
        else if (z < 2.0_dp) then
            p = exp(-1.2337141_dp / z) / sqrt(z) * &
                (2.00012_dp + (0.247105_dp - (0.0649821_dp - &
                (0.0347962_dp - (0.011672_dp - 0.00168691_dp * z) * z) * z) * z) * z)
        else
            p = exp(-exp(1.0776_dp - (2.30695_dp - (0.43424_dp - &
                (0.082433_dp - (0.008056_dp - 0.0003146_dp * z) * z) * z) * z) * z))
        end if
        p = clamp01(p)
    end function ad_inf_fast

    pure real(dp) function ad_error_fix(n, x) result(v)
        integer, intent(in) :: n
        real(dp), intent(in) :: x
        real(dp) :: c, t

        if (x > 0.8_dp) then
            v = (-130.2137_dp + (745.2337_dp - (1705.091_dp - &
                (1950.646_dp - (1116.360_dp - 255.7844_dp * x) * x) * x) * x) * x) / real(n, dp)
            return
        end if
        c = 0.01265_dp + 0.1757_dp / real(n, dp)
        if (x < c) then
            t = x / c
            t = sqrt(t) * (1.0_dp - t) * (49.0_dp * t - 102.0_dp)
            v = t * (0.0037_dp / real(n * n, dp) + 0.00078_dp / real(n, dp) + 0.00006_dp) / real(n, dp)
            return
        end if
        t = (x - c) / (0.8_dp - c)
        t = -0.00022633_dp + (6.54034_dp - (14.6538_dp - &
            (14.458_dp - (8.259_dp - 1.91864_dp * t) * t) * t) * t) * t
        v = t * (0.04213_dp + 0.01365_dp / real(n, dp)) / real(n, dp)
    end function ad_error_fix

    pure real(dp) function ad_cdf_finite(z, n) result(p)
        real(dp), intent(in) :: z
        integer, intent(in) :: n
        real(dp) :: x

        if (z <= 0.0_dp) then
            p = 0.0_dp
            return
        end if
        x = ad_inf_fast(z)
        p = clamp01(x + ad_error_fix(n, x))
    end function ad_cdf_finite

    real(dp) function ad_inf_exact(z) result(p)
        real(dp), intent(in) :: z
        real(dp) :: r, pnew
        integer :: j

        if (z < 0.01_dp) then
            p = 0.0_dp
            return
        end if
        r = 1.0_dp / z
        p = r * ad_f(z, 0)
        do j = 1, 99
            r = r * (0.5_dp - real(j, dp)) / real(j, dp)
            pnew = p + real(4 * j + 1, dp) * r * ad_f(z, j)
            if (abs(pnew - p) <= epsilon(p) * max(1.0_dp, abs(p))) exit
            p = pnew
        end do
        p = clamp01(p)
    end function ad_inf_exact

    real(dp) function ad_f(z, j) result(f)
        real(dp), intent(in) :: z
        integer, intent(in) :: j
        real(dp) :: t, fnew, a, b, c, r
        integer :: i

        t = real((4 * j + 1) * (4 * j + 1), dp) * 1.23370055013617_dp / z
        if (t > 150.0_dp) then
            f = 0.0_dp
            return
        end if
        a = 2.22144146907918_dp * exp(-t) / sqrt(t)
        b = 3.93740248643060_dp * erfc(sqrt(t))
        r = z * 0.125_dp
        f = a + b * r
        do i = 1, 199
            c = ((real(i, dp) - 0.5_dp - t) * b + t * a) / real(i, dp)
            a = b
            b = c
            r = r * z / real(8 * i + 8, dp)
            if (abs(r) < 1.0e-40_dp .or. abs(c) < 1.0e-40_dp) return
            fnew = f + c * r
            if (abs(fnew - f) <= epsilon(f) * max(1.0_dp, abs(f))) return
            f = fnew
        end do
    end function ad_f

    real(dp) function ad_cdf(q, n, fast, lower_tail) result(p)
        real(dp), intent(in) :: q
        integer, intent(in), optional :: n
        logical, intent(in), optional :: fast, lower_tail
        logical :: use_fast, lower

        use_fast = .true.
        if (present(fast)) use_fast = fast
        lower = .true.
        if (present(lower_tail)) lower = lower_tail

        if (q <= 0.0_dp) then
            p = 0.0_dp
        else if (present(n)) then
            if (n <= 0) then
                p = 0.0_dp
            else
                p = ad_cdf_finite(q, n)
            end if
        else if (use_fast) then
            p = ad_inf_fast(q)
        else
            p = ad_inf_exact(q)
        end if
        if (.not. lower) p = 1.0_dp - p
    end function ad_cdf

    real(dp) function ad_quantile(prob, n, fast, lower_tail) result(q)
        real(dp), intent(in) :: prob
        integer, intent(in), optional :: n
        logical, intent(in), optional :: fast, lower_tail
        real(dp) :: target, lo, hi, mid, fm
        logical :: lower
        integer :: iter

        lower = .true.
        if (present(lower_tail)) lower = lower_tail
        target = prob
        if (.not. lower) target = 1.0_dp - target
        if (target <= 0.0_dp) then
            q = 0.0_dp
            return
        end if
        if (target >= 1.0_dp) then
            q = huge(1.0_dp)
            return
        end if

        lo = 0.0_dp
        hi = 1.0_dp
        do
            if (present(n)) then
                fm = ad_cdf(hi, n=n, fast=fast)
            else
                fm = ad_cdf(hi, fast=fast)
            end if
            if (fm >= target) exit
            hi = 2.0_dp * hi
            if (hi > 1.0e8_dp) exit
        end do
        do iter = 1, 100
            mid = 0.5_dp * (lo + hi)
            if (present(n)) then
                fm = ad_cdf(mid, n=n, fast=fast)
            else
                fm = ad_cdf(mid, fast=fast)
            end if
            if (fm < target) then
                lo = mid
            else
                hi = mid
            end if
        end do
        q = 0.5_dp * (lo + hi)
    end function ad_quantile

    real(dp) function ad_statistic(u) result(a)
        real(dp), intent(in) :: u(:)
        real(dp), allocatable :: x(:)
        integer :: i, n

        n = size(u)
        if (n == 0) then
            a = 0.0_dp
            return
        end if
        allocate(x(n))
        x = u
        call sort_real(x)
        a = 0.0_dp
        do i = 1, n
            if (x(i) <= 0.0_dp .or. x(n + 1 - i) >= 1.0_dp) then
                a = huge(1.0_dp)
                return
            end if
            a = a - real(2 * i - 1, dp) * log(x(i) * (1.0_dp - x(n + 1 - i)))
        end do
        a = -real(n, dp) + a / real(n, dp)
    end function ad_statistic

    subroutine ad_test_uniform(u, statistic, p_value)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: statistic, p_value
        integer :: n

        n = size(u)
        statistic = ad_statistic(u)
        p_value = 1.0_dp - ad_cdf(statistic, n=n)
    end subroutine ad_test_uniform

end module goftest_ad
