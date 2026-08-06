! SPDX-License-Identifier: GPL-3.0-only
module rsdc_linalg
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use rsdc_kinds, only: dp
    implicit none
    private
    public :: cholesky_lower, inverse_spd, inverse_matrix, logdet_spd
    public :: is_positive_definite, solve_linear, solve_spd
    public :: sample_covariance, sample_correlation, mean_columns
    public :: standard_deviation, quantile_type7, project_simplex

contains

    subroutine cholesky_lower(a, l, ok)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: l(:, :)
        logical, intent(out) :: ok
        integer :: i, j, k, n
        real(dp) :: s
        n = size(a, 1)
        l = 0.0_dp
        ok = size(a, 2) == n .and. size(l, 1) == n .and. size(l, 2) == n
        if (.not. ok) return
        do i = 1, n
            do j = 1, i
                s = a(i, j)
                do k = 1, j - 1
                    s = s - l(i, k) * l(j, k)
                end do
                if (i == j) then
                    if (.not. ieee_is_finite(s) .or. s <= 0.0_dp) then
                        ok = .false.
                        return
                    end if
                    l(i, j) = sqrt(s)
                else
                    if (l(j, j) <= 0.0_dp) then
                        ok = .false.
                        return
                    end if
                    l(i, j) = s / l(j, j)
                end if
            end do
        end do
    end subroutine cholesky_lower

    logical function is_positive_definite(a, tolerance) result(ok)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(in), optional :: tolerance
        real(dp), allocatable :: l(:, :), b(:, :)
        real(dp) :: tol
        integer :: n
        logical :: chol_ok
        n = size(a, 1)
        if (size(a, 2) /= n) then
            ok = .false.
            return
        end if
        tol = 0.0_dp
        if (present(tolerance)) tol = tolerance
        allocate(l(n, n), b(n, n))
        b = 0.5_dp * (a + transpose(a))
        if (tol > 0.0_dp) b = b - tol * identity_local(n)
        call cholesky_lower(b, l, chol_ok)
        ok = chol_ok
    end function is_positive_definite

    function identity_local(n) result(a)
        integer, intent(in) :: n
        real(dp) :: a(n, n)
        integer :: i
        a = 0.0_dp
        do i = 1, n
            a(i, i) = 1.0_dp
        end do
    end function identity_local

    subroutine solve_lower(l, b, x, ok)
        real(dp), intent(in) :: l(:, :), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        integer :: i, j, n
        real(dp) :: s
        n = size(b)
        ok = size(l, 1) == n .and. size(l, 2) == n .and. size(x) == n
        if (.not. ok) return
        do i = 1, n
            s = b(i)
            do j = 1, i - 1
                s = s - l(i, j) * x(j)
            end do
            if (abs(l(i, i)) <= tiny(1.0_dp)) then
                ok = .false.
                return
            end if
            x(i) = s / l(i, i)
        end do
    end subroutine solve_lower

    subroutine solve_upper(u, b, x, ok)
        real(dp), intent(in) :: u(:, :), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        integer :: i, j, n
        real(dp) :: s
        n = size(b)
        ok = size(u, 1) == n .and. size(u, 2) == n .and. size(x) == n
        if (.not. ok) return
        do i = n, 1, -1
            s = b(i)
            do j = i + 1, n
                s = s - u(i, j) * x(j)
            end do
            if (abs(u(i, i)) <= tiny(1.0_dp)) then
                ok = .false.
                return
            end if
            x(i) = s / u(i, i)
        end do
    end subroutine solve_upper

    subroutine solve_spd(a, b, x, ok)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: l(:, :), y(:)
        integer :: n
        logical :: ok1
        n = size(b)
        allocate(l(n, n), y(n))
        call cholesky_lower(a, l, ok)
        if (.not. ok) return
        call solve_lower(l, b, y, ok1)
        if (.not. ok1) then
            ok = .false.
            return
        end if
        call solve_upper(transpose(l), y, x, ok)
    end subroutine solve_spd

    subroutine inverse_spd(a, ainv, ok)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: ainv(:, :)
        logical, intent(out) :: ok
        real(dp), allocatable :: e(:), x(:)
        integer :: i, n
        logical :: ok1
        n = size(a, 1)
        ok = size(a, 2) == n .and. size(ainv, 1) == n .and. size(ainv, 2) == n
        if (.not. ok) return
        allocate(e(n), x(n))
        do i = 1, n
            e = 0.0_dp
            e(i) = 1.0_dp
            call solve_spd(a, e, x, ok1)
            if (.not. ok1) then
                ok = .false.
                return
            end if
            ainv(:, i) = x
        end do
        ainv = 0.5_dp * (ainv + transpose(ainv))
    end subroutine inverse_spd

    real(dp) function logdet_spd(a, ok) result(value)
        real(dp), intent(in) :: a(:, :)
        logical, intent(out) :: ok
        real(dp), allocatable :: l(:, :)
        integer :: i, n
        n = size(a, 1)
        allocate(l(n, n))
        call cholesky_lower(a, l, ok)
        if (.not. ok) then
            value = huge(1.0_dp)
            return
        end if
        value = 0.0_dp
        do i = 1, n
            value = value + 2.0_dp * log(l(i, i))
        end do
    end function logdet_spd

    subroutine inverse_matrix(a, ainv, ok)
        real(dp), intent(in) :: a(:, :)
        real(dp), intent(out) :: ainv(:, :)
        logical, intent(out) :: ok
        real(dp), allocatable :: aug(:, :), rowtmp(:)
        real(dp) :: pivot, factor
        integer :: i, j, k, p, n
        n = size(a, 1)
        ok = size(a, 2) == n .and. size(ainv, 1) == n .and. size(ainv, 2) == n
        if (.not. ok) return
        allocate(aug(n, 2 * n), rowtmp(2 * n))
        aug = 0.0_dp
        aug(:, 1:n) = a
        do i = 1, n
            aug(i, n + i) = 1.0_dp
        end do
        do i = 1, n
            p = i
            do k = i + 1, n
                if (abs(aug(k, i)) > abs(aug(p, i))) p = k
            end do
            if (abs(aug(p, i)) <= 100.0_dp * epsilon(1.0_dp)) then
                ok = .false.
                return
            end if
            if (p /= i) then
                rowtmp = aug(i, :)
                aug(i, :) = aug(p, :)
                aug(p, :) = rowtmp
            end if
            pivot = aug(i, i)
            aug(i, :) = aug(i, :) / pivot
            do j = 1, n
                if (j == i) cycle
                factor = aug(j, i)
                aug(j, :) = aug(j, :) - factor * aug(i, :)
            end do
        end do
        ainv = aug(:, n + 1:2 * n)
    end subroutine inverse_matrix

    subroutine solve_linear(a, b, x, ok)
        real(dp), intent(in) :: a(:, :), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: ai(:, :)
        integer :: n
        n = size(b)
        allocate(ai(n, n))
        call inverse_matrix(a, ai, ok)
        if (ok) x = matmul(ai, b)
    end subroutine solve_linear

    function mean_columns(x) result(mu)
        real(dp), intent(in) :: x(:, :)
        real(dp) :: mu(size(x, 2))
        if (size(x, 1) == 0) then
            mu = 0.0_dp
        else
            mu = sum(x, dim=1) / real(size(x, 1), dp)
        end if
    end function mean_columns

    subroutine sample_covariance(x, cov)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: cov(:, :)
        real(dp) :: mu(size(x, 2))
        real(dp), allocatable :: xc(:, :)
        integer :: n
        n = size(x, 1)
        mu = mean_columns(x)
        allocate(xc(n, size(x, 2)))
        xc = x - spread(mu, 1, n)
        if (n > 1) then
            cov = matmul(transpose(xc), xc) / real(n - 1, dp)
        else
            cov = 0.0_dp
        end if
    end subroutine sample_covariance

    subroutine sample_correlation(x, cor)
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: cor(:, :)
        real(dp), allocatable :: cov(:, :)
        real(dp), allocatable :: sd(:)
        integer :: i, j, k
        k = size(x, 2)
        allocate(cov(k, k), sd(k))
        call sample_covariance(x, cov)
        do i = 1, k
            sd(i) = sqrt(max(cov(i, i), 0.0_dp))
        end do
        cor = 0.0_dp
        do i = 1, k
            cor(i, i) = 1.0_dp
            do j = i + 1, k
                if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) then
                    cor(i, j) = cov(i, j) / (sd(i) * sd(j))
                    cor(j, i) = cor(i, j)
                end if
            end do
        end do
    end subroutine sample_correlation

    real(dp) function standard_deviation(x) result(sd)
        real(dp), intent(in) :: x(:)
        real(dp) :: mu
        integer :: n
        n = size(x)
        if (n <= 1) then
            sd = 0.0_dp
            return
        end if
        mu = sum(x) / real(n, dp)
        sd = sqrt(sum((x - mu) ** 2) / real(n - 1, dp))
    end function standard_deviation

    real(dp) function quantile_type7(x, prob) result(q)
        real(dp), intent(in) :: x(:), prob
        real(dp), allocatable :: y(:)
        real(dp) :: h, g
        integer :: i, j, n
        n = size(x)
        if (n == 0) then
            q = 0.0_dp
            return
        end if
        allocate(y(n)); y = x
        do i = 2, n
            g = y(i); j = i - 1
            do while (j >= 1)
                if (y(j) <= g) exit
                y(j + 1) = y(j)
                j = j - 1
            end do
            y(j + 1) = g
        end do
        if (prob <= 0.0_dp) then
            q = y(1); return
        else if (prob >= 1.0_dp) then
            q = y(n); return
        end if
        h = 1.0_dp + real(n - 1, dp) * prob
        i = int(floor(h)); g = h - real(i, dp)
        if (i >= n) then
            q = y(n)
        else
            q = (1.0_dp - g) * y(i) + g * y(i + 1)
        end if
    end function quantile_type7

    subroutine project_simplex(v, w)
        real(dp), intent(in) :: v(:)
        real(dp), intent(out) :: w(:)
        real(dp), allocatable :: u(:)
        real(dp) :: cssv, theta, tmp
        integer :: i, j, rho, n
        n = size(v)
        allocate(u(n)); u = v
        do i = 2, n
            tmp = u(i); j = i - 1
            do while (j >= 1)
                if (u(j) >= tmp) exit
                u(j + 1) = u(j)
                j = j - 1
            end do
            u(j + 1) = tmp
        end do
        cssv = 0.0_dp; rho = 1
        do i = 1, n
            cssv = cssv + u(i)
            if (u(i) - (cssv - 1.0_dp) / real(i, dp) > 0.0_dp) rho = i
        end do
        theta = (sum(u(1:rho)) - 1.0_dp) / real(rho, dp)
        w = max(v - theta, 0.0_dp)
    end subroutine project_simplex
end module rsdc_linalg
