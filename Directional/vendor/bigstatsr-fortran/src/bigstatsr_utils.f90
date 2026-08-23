! SPDX-License-Identifier: GPL-3.0-only
module bigstatsr_utils
    use bigstatsr_kinds, only: dp
    implicit none
    private
    public :: sigmoid, soft_threshold, normal_cdf, normal_quantile, sort_pairs
    public :: mean_dp, sample_var, correlation, median_dp

contains

    pure elemental real(dp) function sigmoid(x) result(p)
        real(dp), intent(in) :: x
        if (x >= 0.0_dp) then
            p = 1.0_dp / (1.0_dp + exp(-x))
        else
            p = exp(x) / (1.0_dp + exp(x))
        end if
    end function sigmoid

    pure elemental real(dp) function soft_threshold(z, l1, l2, v) result(ans)
        real(dp), intent(in) :: z, l1, l2
        real(dp), intent(in), optional :: v
        real(dp) :: den
        den = 1.0_dp + l2
        if (present(v)) den = v * (1.0_dp + l2)
        if (z > l1) then
            ans = (z - l1) / den
        else if (z < -l1) then
            ans = (z + l1) / den
        else
            ans = 0.0_dp
        end if
    end function soft_threshold

    pure elemental real(dp) function normal_cdf(x) result(p)
        real(dp), intent(in) :: x
        p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
    end function normal_cdf

    pure elemental real(dp) function normal_quantile(p) result(x)
        real(dp), intent(in) :: p
        real(dp), parameter :: a(6) = [ &
            -3.969683028665376d+01, 2.209460984245205d+02, &
            -2.759285104469687d+02, 1.383577518672690d+02, &
            -3.066479806614716d+01, 2.506628277459239d+00 ]
        real(dp), parameter :: b(5) = [ &
            -5.447609879822406d+01, 1.615858368580409d+02, &
            -1.556989798598866d+02, 6.680131188771972d+01, &
            -1.328068155288572d+01 ]
        real(dp), parameter :: c(6) = [ &
            -7.784894002430293d-03, -3.223964580411365d-01, &
            -2.400758277161838d+00, -2.549732539343734d+00, &
             4.374664141464968d+00,  2.938163982698783d+00 ]
        real(dp), parameter :: d(4) = [ &
             7.784695709041462d-03, 3.224671290700398d-01, &
             2.445134137142996d+00, 3.754408661907416d+00 ]
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

    pure real(dp) function mean_dp(x) result(ans)
        real(dp), intent(in) :: x(:)
        if (size(x) == 0) then
            ans = 0.0_dp
        else
            ans = sum(x) / real(size(x), dp)
        end if
    end function mean_dp

    pure real(dp) function sample_var(x) result(ans)
        real(dp), intent(in) :: x(:)
        real(dp) :: mu
        if (size(x) < 2) then
            ans = 0.0_dp
            return
        end if
        mu = mean_dp(x)
        ans = sum((x - mu)**2) / real(size(x)-1, dp)
    end function sample_var

    pure real(dp) function correlation(x, y) result(r)
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: mx, my, sx, sy
        mx = mean_dp(x)
        my = mean_dp(y)
        sx = sqrt(sum((x-mx)**2))
        sy = sqrt(sum((y-my)**2))
        if (sx <= 0.0_dp .or. sy <= 0.0_dp) then
            r = 0.0_dp
        else
            r = dot_product(x-mx, y-my) / (sx*sy)
        end if
    end function correlation

    function median_dp(x) result(ans)
        real(dp), intent(in) :: x(:)
        real(dp) :: ans
        real(dp), allocatable :: work(:)
        integer :: n
        n = size(x)
        if (n == 0) then
            ans = 0.0_dp
            return
        end if
        work = x
        call sort_real(work)
        if (mod(n,2) == 1) then
            ans = work((n+1)/2)
        else
            ans = 0.5_dp * (work(n/2) + work(n/2+1))
        end if
    end function median_dp

    subroutine sort_pairs(x, y)
        real(dp), intent(inout) :: x(:)
        integer, intent(inout) :: y(:)
        integer :: i, j, yi
        real(dp) :: xi
        do i = 2, size(x)
            xi = x(i)
            yi = y(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) < xi) exit
                if (abs(x(j) - xi) <= 0.0_dp .and. y(j) <= yi) exit
                x(j+1) = x(j)
                y(j+1) = y(j)
                j = j - 1
            end do
            x(j+1) = xi
            y(j+1) = yi
        end do
    end subroutine sort_pairs

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: v
        do i = 2, size(x)
            v = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= v) exit
                x(j+1) = x(j)
                j = j - 1
            end do
            x(j+1) = v
        end do
    end subroutine sort_real

end module bigstatsr_utils
