! Computational translation of the R package bivpois.
! License: GPL-2.0-or-later.
module bivpois_math
    use iso_fortran_env, only : int64
    use bivpois_kinds, only : dp
    implicit none
    private

    public :: seed_rng, rpois_one
    public :: mean_counts, variance_counts, covariance_counts, correlation_counts
    public :: normal_upper_tail, chisq1_upper_tail


contains

    subroutine seed_rng(seed)
        integer, intent(in) :: seed
        integer :: n, i
        integer, allocatable :: put(:)
        integer(int64) :: s

        call random_seed(size=n)
        allocate(put(n))
        s = int(seed, int64)
        if (s == 0_int64) s = 104729_int64
        do i = 1, n
            s = modulo(1103515245_int64 * s + 12345_int64 + 97_int64 * i, 2147483647_int64)
            if (s <= 0_int64) s = s + 2147483646_int64
            put(i) = int(s)
        end do
        call random_seed(put=put)
    end subroutine seed_rng

    integer function rpois_one(lambda) result(k)
        real(dp), intent(in) :: lambda
        real(dp) :: l, p, u, v, us
        real(dp) :: slam, loglam, b, a, inv_alpha, vr, lhs, rhs

        if (lambda < 0.0_dp) error stop "rpois_one: lambda must be nonnegative"
        if (lambda <= 0.0_dp) then
            k = 0
            return
        end if

        if (lambda < 30.0_dp) then
            l = exp(-lambda)
            p = 1.0_dp
            k = 0
            do
                call random_number(u)
                p = p * u
                if (p <= l) exit
                k = k + 1
            end do
            return
        end if

        ! PTRS transformed rejection method of Hoermann for large Poisson means.
        slam = sqrt(lambda)
        loglam = log(lambda)
        b = 0.931_dp + 2.53_dp * slam
        a = -0.059_dp + 0.02483_dp * b
        inv_alpha = 1.1239_dp + 1.1328_dp / (b - 3.4_dp)
        vr = 0.9277_dp - 3.6224_dp / (b - 2.0_dp)

        do
            call random_number(u)
            call random_number(v)
            u = u - 0.5_dp
            us = 0.5_dp - abs(u)
            if (us <= 0.0_dp) cycle
            k = floor((2.0_dp * a / us + b) * u + lambda + 0.43_dp)
            if (us >= 0.07_dp .and. v <= vr .and. k >= 0) return
            if (k < 0) cycle
            if (us < 0.013_dp) then
                if (v > us) cycle
            end if
            lhs = log(v * inv_alpha / (a / (us * us) + b))
            rhs = -lambda + real(k, dp) * loglam - log_gamma(real(k + 1, dp))
            if (lhs <= rhs) return
        end do
    end function rpois_one

    pure real(dp) function mean_counts(x) result(m)
        integer, intent(in) :: x(:)
        if (size(x) == 0) error stop "mean_counts: empty input"
        m = sum(real(x, dp)) / real(size(x), dp)
    end function mean_counts

    pure real(dp) function variance_counts(x) result(v)
        integer, intent(in) :: x(:)
        real(dp) :: m
        integer :: n
        n = size(x)
        if (n < 2) error stop "variance_counts: need at least two observations"
        m = mean_counts(x)
        v = sum((real(x, dp) - m)**2) / real(n - 1, dp)
    end function variance_counts

    pure real(dp) function covariance_counts(x, y) result(c)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: mx, my
        integer :: n
        n = size(x)
        if (size(y) /= n) error stop "covariance_counts: size mismatch"
        if (n < 2) error stop "covariance_counts: need at least two observations"
        mx = mean_counts(x)
        my = mean_counts(y)
        c = sum((real(x, dp) - mx) * (real(y, dp) - my)) / real(n - 1, dp)
    end function covariance_counts

    pure real(dp) function correlation_counts(x, y) result(r)
        integer, intent(in) :: x(:), y(:)
        real(dp) :: vx, vy
        vx = variance_counts(x)
        vy = variance_counts(y)
        if (vx <= 0.0_dp .or. vy <= 0.0_dp) then
            r = 0.0_dp
        else
            r = covariance_counts(x, y) / sqrt(vx * vy)
        end if
    end function correlation_counts

    pure real(dp) function normal_upper_tail(z) result(p)
        real(dp), intent(in) :: z
        p = 0.5_dp * erfc(z / sqrt(2.0_dp))
    end function normal_upper_tail

    pure real(dp) function chisq1_upper_tail(x) result(p)
        real(dp), intent(in) :: x
        if (x <= 0.0_dp) then
            p = 1.0_dp
        else
            p = erfc(sqrt(0.5_dp * x))
        end if
    end function chisq1_upper_tail


end module bivpois_math
