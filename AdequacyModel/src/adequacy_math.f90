! SPDX-License-Identifier: GPL-2.0-or-later
module adequacy_math
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use adequacy_kinds, only: dp
    use adequacy_interfaces, only: objective_fn
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)

    public :: normal_cdf, normal_quantile, sort_real, median_real
    public :: mean_real, sample_variance, numeric_gradient, numeric_hessian
    public :: invert_matrix, kolmogorov_pvalue

contains

    pure function mean_real(x) result(ans)
        real(dp), intent(in) :: x(:)
        real(dp) :: ans

        if (size(x) == 0) then
            ans = 0.0_dp
        else
            ans = sum(x) / real(size(x), dp)
        end if
    end function mean_real

    pure function sample_variance(x) result(ans)
        real(dp), intent(in) :: x(:)
        real(dp) :: ans, mu
        integer :: n

        n = size(x)
        if (n < 2) then
            ans = 0.0_dp
            return
        end if
        mu = mean_real(x)
        ans = sum((x - mu)**2) / real(n - 1, dp)
    end function sample_variance

    subroutine sort_real(x)
        real(dp), intent(inout) :: x(:)
        integer :: i, j
        real(dp) :: key

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j + 1) = x(j)
                j = j - 1
            end do
            x(j + 1) = key
        end do
    end subroutine sort_real

    function median_real(x) result(ans)
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
        if (mod(n, 2) == 1) then
            ans = work((n + 1) / 2)
        else
            ans = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
        end if
    end function median_real

    pure elemental function normal_cdf(x) result(p)
        real(dp), intent(in) :: x
        real(dp) :: p
        p = 0.5_dp * erfc(-x / sqrt(2.0_dp))
    end function normal_cdf

    pure elemental function normal_quantile(p) result(x)
        real(dp), intent(in) :: p
        real(dp) :: x, q, r
        real(dp), parameter :: a(6) = [ &
            -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
            -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
            -3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
        real(dp), parameter :: b(5) = [ &
            -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
            -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
            -1.328068155288572e+01_dp ]
        real(dp), parameter :: c(6) = [ &
            -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
            -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
             4.374664141464968e+00_dp, 2.938163982698783e+00_dp ]
        real(dp), parameter :: d(4) = [ &
             7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
             2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
        real(dp), parameter :: plow = 0.02425_dp
        real(dp), parameter :: phigh = 1.0_dp - plow

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
            r = q * q
            x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
                (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
        else
            q = sqrt(-2.0_dp * log(1.0_dp - p))
            x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
                ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
        end if
    end function normal_quantile

    subroutine numeric_gradient(func, par, data, grad)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: par(:), data(:)
        real(dp), intent(out) :: grad(:)
        real(dp), allocatable :: xp(:), xm(:)
        real(dp) :: h
        integer :: i

        allocate(xp(size(par)), xm(size(par)))
        do i = 1, size(par)
            xp = par
            xm = par
            h = sqrt(epsilon(1.0_dp)) * max(1.0_dp, abs(par(i)))
            xp(i) = xp(i) + h
            xm(i) = xm(i) - h
            grad(i) = (func(xp, data) - func(xm, data)) / (2.0_dp * h)
        end do
    end subroutine numeric_gradient

    subroutine numeric_hessian(func, par, data, hess)
        procedure(objective_fn) :: func
        real(dp), intent(in) :: par(:), data(:)
        real(dp), intent(out) :: hess(:, :)
        real(dp), allocatable :: xp(:), xm(:), xpp(:), xpm(:), xmp(:), xmm(:)
        real(dp) :: hi, hj, f0
        integer :: i, j, n

        n = size(par)
        allocate(xp(n), xm(n), xpp(n), xpm(n), xmp(n), xmm(n))
        f0 = func(par, data)
        do i = 1, n
            hi = epsilon(1.0_dp)**0.25_dp * max(1.0_dp, abs(par(i)))
            xp = par
            xm = par
            xp(i) = xp(i) + hi
            xm(i) = xm(i) - hi
            hess(i, i) = (func(xp, data) - 2.0_dp*f0 + func(xm, data)) / (hi*hi)
            do j = i + 1, n
                hj = epsilon(1.0_dp)**0.25_dp * max(1.0_dp, abs(par(j)))
                xpp = par
                xpm = par
                xmp = par
                xmm = par
                xpp(i) = xpp(i) + hi
                xpp(j) = xpp(j) + hj
                xpm(i) = xpm(i) + hi
                xpm(j) = xpm(j) - hj
                xmp(i) = xmp(i) - hi
                xmp(j) = xmp(j) + hj
                xmm(i) = xmm(i) - hi
                xmm(j) = xmm(j) - hj
                hess(i, j) = (func(xpp, data) - func(xpm, data) - &
                              func(xmp, data) + func(xmm, data)) / (4.0_dp*hi*hj)
                hess(j, i) = hess(i, j)
            end do
        end do
    end subroutine numeric_hessian

    subroutine invert_matrix(a, ainv, ok)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: ainv(:, :)
        logical, intent(out) :: ok
        real(dp), allocatable :: aug(:, :), tmp(:)
        real(dp) :: pivot
        integer :: n, i, j, k

        n = size(a, 1)
        if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
            ok = .false.
            return
        end if
        allocate(aug(n, 2*n), tmp(2*n))
        aug(:, 1:n) = a
        aug(:, n+1:2*n) = 0.0_dp
        do i = 1, n
            aug(i, n+i) = 1.0_dp
        end do

        do i = 1, n
            k = i
            do j = i + 1, n
                if (abs(aug(j, i)) > abs(aug(k, i))) k = j
            end do
            if (abs(aug(k, i)) <= 100.0_dp*epsilon(1.0_dp)) then
                ok = .false.
                ainv = 0.0_dp
                return
            end if
            if (k /= i) then
                tmp = aug(i, :)
                aug(i, :) = aug(k, :)
                aug(k, :) = tmp
            end if
            pivot = aug(i, i)
            aug(i, :) = aug(i, :) / pivot
            do j = 1, n
                if (j == i) cycle
                pivot = aug(j, i)
                aug(j, :) = aug(j, :) - pivot * aug(i, :)
            end do
        end do
        ainv = aug(:, n+1:2*n)
        ok = .true.
    end subroutine invert_matrix

    function kolmogorov_pvalue(d, n) result(p)
        real(dp), intent(in) :: d
        integer, intent(in) :: n
        real(dp) :: p, cdf, lambda, term, s
        integer :: k, m

        if (n <= 0 .or. d <= 0.0_dp) then
            p = 1.0_dp
            return
        end if
        if (d >= 1.0_dp) then
            p = 0.0_dp
            return
        end if

        k = int(real(n,dp)*d) + 1
        m = 2*k - 1
        if (n <= 140 .and. m <= 60) then
            cdf = ks_exact_cdf(d, n)
            p = min(1.0_dp, max(0.0_dp, 1.0_dp-cdf))
            return
        end if

        lambda = (sqrt(real(n, dp)) + 0.12_dp + 0.11_dp/sqrt(real(n, dp))) * d
        s = 0.0_dp
        do k = 1, 1000
            term = 2.0_dp * (-1.0_dp)**(k - 1) * exp(-2.0_dp*real(k*k, dp)*lambda*lambda)
            s = s + term
            if (abs(term) < 1.0e-14_dp) exit
        end do
        if (.not. ieee_is_finite(s)) s = 0.0_dp
        p = min(1.0_dp, max(0.0_dp, s))
    end function kolmogorov_pvalue

    function ks_exact_cdf(d, n) result(cdf)
        real(dp), intent(in) :: d
        integer, intent(in) :: n
        real(dp) :: cdf, nd, hval, scale
        real(dp), allocatable :: h(:, :), q(:, :)
        integer :: i, j, g, k, m

        nd = real(n,dp)*d
        k = int(nd) + 1
        m = 2*k - 1
        hval = real(k,dp) - nd
        allocate(h(m,m), q(m,m))
        do i = 1, m
            do j = 1, m
                if (i-j+1 < 0) then
                    h(i,j) = 0.0_dp
                else
                    h(i,j) = 1.0_dp
                end if
            end do
        end do
        do i = 1, m
            h(i,1) = h(i,1) - hval**i
            h(m,i) = h(m,i) - hval**(m-i+1)
        end do
        if (2.0_dp*hval - 1.0_dp > 0.0_dp) then
            h(m,1) = h(m,1) + (2.0_dp*hval - 1.0_dp)**m
        end if
        do i = 1, m
            do j = 1, m
                if (i-j+1 > 0) then
                    do g = 1, i-j+1
                        h(i,j) = h(i,j) / real(g,dp)
                    end do
                end if
            end do
        end do
        call matrix_power(h, n, q)
        scale = 1.0_dp
        do i = 1, n
            scale = scale * real(i,dp) / real(n,dp)
        end do
        cdf = min(1.0_dp, max(0.0_dp, q(k,k)*scale))
    end function ks_exact_cdf

    subroutine matrix_power(a, power, out)
        real(dp), intent(in) :: a(:, :)
        integer, intent(in) :: power
        real(dp), intent(out) :: out(:, :)
        real(dp), allocatable :: base(:, :), tmp(:, :)
        integer :: e, i, n

        n = size(a,1)
        allocate(base(n,n), tmp(n,n))
        out = 0.0_dp
        do i = 1, n
            out(i,i) = 1.0_dp
        end do
        base = a
        e = power
        do while (e > 0)
            if (mod(e,2) == 1) then
                tmp = matmul(out, base)
                out = tmp
            end if
            e = e / 2
            if (e > 0) then
                tmp = matmul(base, base)
                base = tmp
            end if
        end do
    end subroutine matrix_power

end module adequacy_math
