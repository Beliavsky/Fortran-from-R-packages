! SPDX-License-Identifier: GPL-2.0-or-later
module evir_math
    use, intrinsic :: iso_fortran_env, only : int64
    use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan, ieee_is_finite
    use evir_kinds, only : dp
    use evir_types, only : evir_rng
    implicit none
    private

    public :: mean_value, sample_variance, sort_ascending, sort_descending
    public :: normal_quantile, chi_square_quantile_df1, ppoints
    public :: seed_rng, random_uniform, random_normal
    public :: invert_matrix, harmonic_numbers, cumulative_max
    public :: is_finite_vector, safe_nan

contains

    pure real(dp) function safe_nan() result(x)
        x = ieee_value(0.0_dp, ieee_quiet_nan)
    end function safe_nan

    pure logical function is_finite_vector(x) result(ok)
        real(dp), intent(in) :: x(:)
        integer :: i
        ok = .true.
        do i = 1, size(x)
            if (.not. ieee_is_finite(x(i))) then
                ok = .false.
                return
            end if
        end do
    end function is_finite_vector

    pure real(dp) function mean_value(x) result(m)
        real(dp), intent(in) :: x(:)
        if (size(x) == 0) then
            m = safe_nan()
        else
            m = sum(x) / real(size(x), dp)
        end if
    end function mean_value

    pure real(dp) function sample_variance(x) result(v)
        real(dp), intent(in) :: x(:)
        real(dp) :: m
        if (size(x) < 2) then
            v = safe_nan()
            return
        end if
        m = mean_value(x)
        v = sum((x - m)**2) / real(size(x) - 1, dp)
    end function sample_variance

    subroutine sort_ascending(x)
        real(dp), intent(inout) :: x(:)
        if (size(x) > 1) call quicksort(x, 1, size(x))
    contains
        recursive subroutine quicksort(a, left, right)
            real(dp), intent(inout) :: a(:)
            integer, intent(in) :: left, right
            integer :: i, j
            real(dp) :: pivot, tmp
            if (left >= right) return
            i = left
            j = right
            pivot = a((left + right) / 2)
            do
                do while (a(i) < pivot)
                    i = i + 1
                end do
                do while (a(j) > pivot)
                    j = j - 1
                end do
                if (i <= j) then
                    tmp = a(i)
                    a(i) = a(j)
                    a(j) = tmp
                    i = i + 1
                    j = j - 1
                end if
                if (i > j) exit
            end do
            if (left < j) call quicksort(a, left, j)
            if (i < right) call quicksort(a, i, right)
        end subroutine quicksort
    end subroutine sort_ascending

    subroutine sort_descending(x)
        real(dp), intent(inout) :: x(:)
        real(dp) :: tmp
        integer :: i, n
        call sort_ascending(x)
        n = size(x)
        do i = 1, n / 2
            tmp = x(i)
            x(i) = x(n + 1 - i)
            x(n + 1 - i) = tmp
        end do
    end subroutine sort_descending

    pure real(dp) function normal_quantile(p) result(x)
        real(dp), intent(in) :: p
        real(dp), parameter :: a(6) = [ &
            -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
            -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
            -3.066479806614716e+01_dp, 2.506628277459239e+00_dp]
        real(dp), parameter :: b(5) = [ &
            -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
            -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
            -1.328068155288572e+01_dp]
        real(dp), parameter :: c(6) = [ &
            -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
            -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
             4.374664141464968e+00_dp, 2.938163982698783e+00_dp]
        real(dp), parameter :: d(4) = [ &
             7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
             2.445134137142996e+00_dp, 3.754408661907416e+00_dp]
        real(dp), parameter :: plow = 0.02425_dp
        real(dp), parameter :: phigh = 1.0_dp - plow
        real(dp) :: q, r

        if (p <= 0.0_dp) then
            x = -huge(1.0_dp)
        else if (p >= 1.0_dp) then
            x = huge(1.0_dp)
        else if (p < plow) then
            q = sqrt(-2.0_dp * log(p))
            x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
                ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
        else if (p <= phigh) then
            q = p - 0.5_dp
            r = q*q
            x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
                (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
        else
            q = sqrt(-2.0_dp * log(1.0_dp-p))
            x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
                ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
        end if
    end function normal_quantile

    pure real(dp) function chi_square_quantile_df1(p) result(x)
        real(dp), intent(in) :: p
        real(dp) :: z
        z = normal_quantile(0.5_dp * (1.0_dp + p))
        x = z*z
    end function chi_square_quantile_df1

    pure real(dp) function ppoints(i, n) result(p)
        integer, intent(in) :: i, n
        real(dp) :: a
        if (n <= 10) then
            a = 0.375_dp
        else
            a = 0.5_dp
        end if
        p = (real(i, dp) - a) / (real(n, dp) + 1.0_dp - 2.0_dp*a)
    end function ppoints

    subroutine seed_rng(rng, seed)
        type(evir_rng), intent(inout) :: rng
        integer(int64), intent(in) :: seed
        if (seed == 0_int64) then
            rng%state = 88172645463325252_int64
        else
            rng%state = seed
        end if
    end subroutine seed_rng

    real(dp) function random_uniform(rng) result(u)
        type(evir_rng), intent(inout) :: rng
        integer(int64) :: x
        x = rng%state
        x = ieor(x, shiftl(x, 13))
        x = ieor(x, shiftr(x, 7))
        x = ieor(x, shiftl(x, 17))
        rng%state = x
        u = real(iand(x, int(z'7FFFFFFFFFFFFFFF', int64)), dp) / &
            real(huge(1_int64), dp)
        if (u <= 0.0_dp) u = epsilon(1.0_dp)
        if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
    end function random_uniform

    real(dp) function random_normal(rng) result(z)
        type(evir_rng), intent(inout) :: rng
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp) :: u1, u2
        u1 = random_uniform(rng)
        u2 = random_uniform(rng)
        z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp*pi*u2)
    end function random_normal

    subroutine invert_matrix(a, ainv, status)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: ainv(:, :)
        integer, intent(out) :: status
        real(dp), allocatable :: aug(:, :), rowtmp(:)
        real(dp) :: pivot, factor
        integer :: n, i, k, imax

        n = size(a, 1)
        status = 0
        ainv = 0.0_dp
        if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
            status = 1
            return
        end if
        allocate(aug(n, 2*n), rowtmp(2*n))
        aug(:, 1:n) = a
        aug(:, n+1:2*n) = 0.0_dp
        do i = 1, n
            aug(i, n+i) = 1.0_dp
        end do
        do k = 1, n
            imax = k
            do i = k+1, n
                if (abs(aug(i, k)) > abs(aug(imax, k))) imax = i
            end do
            if (abs(aug(imax, k)) <= 100.0_dp*epsilon(1.0_dp)) then
                status = 2
                ainv = safe_nan()
                return
            end if
            if (imax /= k) then
                rowtmp = aug(k, :)
                aug(k, :) = aug(imax, :)
                aug(imax, :) = rowtmp
            end if
            pivot = aug(k, k)
            aug(k, :) = aug(k, :) / pivot
            do i = 1, n
                if (i == k) cycle
                factor = aug(i, k)
                aug(i, :) = aug(i, :) - factor*aug(k, :)
            end do
        end do
        ainv = aug(:, n+1:2*n)
    end subroutine invert_matrix

    pure subroutine harmonic_numbers(n, h, h2)
        integer, intent(in) :: n
        real(dp), intent(out) :: h(n), h2(n)
        integer :: i
        h(1) = 1.0_dp
        h2(1) = 1.0_dp
        do i = 2, n
            h(i) = h(i-1) + 1.0_dp/real(i, dp)
            h2(i) = h2(i-1) + 1.0_dp/real(i*i, dp)
        end do
    end subroutine harmonic_numbers

    pure subroutine cumulative_max(x, out)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: out(size(x))
        integer :: i
        if (size(x) == 0) return
        out(1) = x(1)
        do i = 2, size(x)
            out(i) = max(out(i-1), x(i))
        end do
    end subroutine cumulative_max

end module evir_math
