! SPDX-License-Identifier: GPL-3.0-only
module rsdc_parameters
    use rsdc_kinds, only: dp
    use rsdc_types, only: rsdc_const, rsdc_nox, rsdc_tvtp
    implicit none
    private
    public :: lower_tri_size, correlations_to_matrix, matrix_to_correlations
    public :: partial_to_correlation, correlation_to_partial
    public :: transition_at, fixed_transition_from_parameters
    public :: expected_parameter_count, unpack_natural_parameters

contains

    pure integer function lower_tri_size(k) result(n)
        integer, intent(in) :: k
        n = k * (k - 1) / 2
    end function lower_tri_size

    subroutine correlations_to_matrix(rho, k, r)
        real(dp), intent(in) :: rho(:)
        integer, intent(in) :: k
        real(dp), intent(out) :: r(k, k)
        integer :: i, j, q
        r = 0.0_dp
        do i = 1, k
            r(i, i) = 1.0_dp
        end do
        q = 0
        do j = 1, k - 1
            do i = j + 1, k
                q = q + 1
                r(i, j) = rho(q)
                r(j, i) = rho(q)
            end do
        end do
    end subroutine correlations_to_matrix

    subroutine matrix_to_correlations(r, rho)
        real(dp), intent(in) :: r(:, :)
        real(dp), intent(out) :: rho(:)
        integer :: i, j, q, k
        k = size(r, 1)
        q = 0
        do j = 1, k - 1
            do i = j + 1, k
                q = q + 1
                rho(q) = r(i, j)
            end do
        end do
    end subroutine matrix_to_correlations

    subroutine partial_to_correlation(z, k, r)
        real(dp), intent(in) :: z(:)
        integer, intent(in) :: k
        real(dp), intent(out) :: r(k, k)
        real(dp) :: zz(k, k), rho
        integer :: i, j, l, q
        zz = 0.0_dp
        q = 0
        do i = 1, k - 1
            do j = i + 1, k
                q = q + 1
                zz(j, i) = max(min(z(q), 0.999999_dp), -0.999999_dp)
            end do
        end do
        r = 0.0_dp
        do i = 1, k
            r(i, i) = 1.0_dp
        end do
        do j = 2, k
            do i = 1, j - 1
                rho = zz(j, i)
                if (i > 1) then
                    do l = i - 1, 1, -1
                        rho = zz(i, l) * zz(j, l) + rho * &
                            sqrt(max(1.0_dp - zz(i, l) ** 2, 0.0_dp) * &
                                 max(1.0_dp - zz(j, l) ** 2, 0.0_dp))
                    end do
                end if
                r(j, i) = rho
                r(i, j) = rho
            end do
        end do
    end subroutine partial_to_correlation

    subroutine correlation_to_partial(r, z)
        real(dp), intent(in) :: r(:, :)
        real(dp), intent(out) :: z(:)
        real(dp), allocatable :: s(:, :), s2(:, :), zz(:, :)
        real(dp) :: denom
        integer :: i, j, a, b, q, k
        k = size(r, 1)
        allocate(s(k, k), s2(k, k), zz(k, k))
        s = r
        zz = 0.0_dp
        do i = 1, k - 1
            do j = i + 1, k
                zz(j, i) = max(min(s(j, i), 0.999999_dp), -0.999999_dp)
            end do
            if (i < k - 1) then
                s2 = s
                do a = i + 1, k
                    do b = i + 1, k
                        denom = sqrt(max((1.0_dp - s(a, i) ** 2) * &
                                         (1.0_dp - s(b, i) ** 2), tiny(1.0_dp)))
                        s2(a, b) = (s(a, b) - s(a, i) * s(b, i)) / denom
                    end do
                end do
                s = s2
            end if
        end do
        q = 0
        do i = 1, k - 1
            do j = i + 1, k
                q = q + 1
                z(q) = zz(j, i)
            end do
        end do
    end subroutine correlation_to_partial

    pure real(dp) function logistic(x) result(p)
        real(dp), intent(in) :: x
        if (x >= 0.0_dp) then
            p = 1.0_dp / (1.0_dp + exp(-x))
        else
            p = exp(x) / (1.0_dp + exp(x))
        end if
    end function logistic

    subroutine transition_at(beta, x, pmat)
        real(dp), intent(in) :: beta(:, :), x(:)
        real(dp), intent(out) :: pmat(:, :)
        integer :: i, j, n, p, lo, hi
        real(dp) :: pii, c, denom
        real(dp), allocatable :: logits(:), e(:)
        n = size(beta, 1)
        p = size(x)
        pmat = 0.0_dp
        if (n == 1) then
            pmat(1, 1) = 1.0_dp
            return
        end if
        if (n == 2) then
            do i = 1, n
                pii = logistic(dot_product(x, beta(i, 1:p)))
                pmat(i, :) = 1.0_dp - pii
                pmat(i, i) = pii
            end do
        else
            allocate(logits(n), e(n))
            do i = 1, n
                logits = 0.0_dp
                do j = 1, n - 1
                    lo = (j - 1) * p + 1
                    hi = j * p
                    logits(j) = dot_product(x, beta(i, lo:hi))
                end do
                c = maxval(logits)
                e = exp(logits - c)
                denom = sum(e)
                pmat(i, :) = e / denom
            end do
        end if
    end subroutine transition_at

    subroutine fixed_transition_from_parameters(trans, n, pmat, ok)
        real(dp), intent(in) :: trans(:)
        integer, intent(in) :: n
        real(dp), intent(out) :: pmat(n, n)
        logical, intent(out) :: ok
        real(dp) :: last
        integer :: i, lo, hi
        ok = .true.
        pmat = 0.0_dp
        if (n == 1) then
            pmat(1, 1) = 1.0_dp
        else if (n == 2) then
            if (size(trans) < 2 .or. any(trans(1:2) < 0.0_dp) .or. any(trans(1:2) > 1.0_dp)) then
                ok = .false.; return
            end if
            pmat(1, :) = [trans(1), 1.0_dp - trans(1)]
            pmat(2, :) = [1.0_dp - trans(2), trans(2)]
        else
            if (size(trans) /= n * (n - 1)) then
                ok = .false.; return
            end if
            do i = 1, n
                lo = (i - 1) * (n - 1) + 1
                hi = i * (n - 1)
                last = 1.0_dp - sum(trans(lo:hi))
                if (any(trans(lo:hi) < 0.0_dp) .or. last < 0.0_dp) then
                    ok = .false.; return
                end if
                pmat(i, 1:n - 1) = trans(lo:hi)
                pmat(i, n) = last
            end do
        end if
    end subroutine fixed_transition_from_parameters

    pure integer function expected_parameter_count(method, n, k, p) result(np)
        integer, intent(in) :: method, n, k, p
        integer :: c
        c = lower_tri_size(k)
        select case (method)
        case (rsdc_const)
            np = c
        case (rsdc_nox)
            np = n * (n - 1) + n * c
        case (rsdc_tvtp)
            if (n == 2) then
                np = n * p + n * c
            else
                np = n * (n - 1) * p + n * c
            end if
        case default
            np = 0
        end select
    end function expected_parameter_count

    subroutine unpack_natural_parameters(parameters, method, n, k, p, beta, rho, pmat, ok)
        real(dp), intent(in) :: parameters(:)
        integer, intent(in) :: method, n, k, p
        real(dp), allocatable, intent(out) :: beta(:, :), rho(:, :), pmat(:, :)
        logical, intent(out) :: ok
        integer :: c, nh, nb
        c = lower_tri_size(k)
        ok = size(parameters) == expected_parameter_count(method, n, k, p)
        if (.not. ok) return
        select case (method)
        case (rsdc_const)
            allocate(beta(0, 0), rho(1, c), pmat(1, 1))
            rho(1, :) = parameters
            pmat(1, 1) = 1.0_dp
        case (rsdc_nox)
            nh = n * (n - 1)
            allocate(beta(0, 0), rho(n, c), pmat(n, n))
            call fixed_transition_from_parameters(parameters(1:nh), n, pmat, ok)
            if (.not. ok) return
            rho = reshape(parameters(nh + 1:), [n, c], order=[2, 1])
        case (rsdc_tvtp)
            nb = merge(n * p, n * (n - 1) * p, n == 2)
            allocate(beta(n, nb / n), rho(n, c), pmat(0, 0))
            beta = reshape(parameters(1:nb), shape(beta), order=[2, 1])
            rho = reshape(parameters(nb + 1:), [n, c], order=[2, 1])
        case default
            ok = .false.
        end select
    end subroutine unpack_natural_parameters
end module rsdc_parameters
