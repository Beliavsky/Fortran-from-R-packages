! SPDX-License-Identifier: GPL-2.0-or-later
module mla_lmm
    use mla_kinds, only : dp
    implicit none
    private

    real(dp), parameter :: pi = acos(-1.0_dp)
    public :: loglik_lmm, grad_lmm

contains

    function loglik_lmm(b, y, x, ni) result(loglik)
        real(dp), intent(in) :: b(:), y(:), x(:, :)
        integer, intent(in) :: ni(:)
        real(dp) :: loglik
        real(dp), allocatable :: r(:), vinvr(:)
        real(dp) :: alpha, sigma, s2, a2, den, q, logdet
        integer :: i, j0, n_i, p

        p = size(x, 2)
        alpha = b(p + 1)
        sigma = b(p + 2)
        s2 = sigma * sigma
        a2 = alpha * alpha
        loglik = 0.0_dp
        j0 = 0
        if (s2 <= tiny(1.0_dp)) then
            loglik = -huge(1.0_dp)
            return
        end if

        do i = 1, size(ni)
            n_i = ni(i)
            allocate(r(n_i), vinvr(n_i))
            r = y(j0 + 1:j0 + n_i) - matmul(x(j0 + 1:j0 + n_i, :), b(1:p))
            den = s2 + real(n_i, dp) * a2
            vinvr = r / s2 - (a2 / (s2 * den)) * sum(r)
            q = dot_product(r, vinvr)
            logdet = real(n_i - 1, dp) * log(s2) + log(den)
            loglik = loglik - 0.5_dp * (real(n_i, dp) * log(2.0_dp * pi) + logdet + q)
            deallocate(r, vinvr)
            j0 = j0 + n_i
        end do
    end function loglik_lmm

    subroutine grad_lmm(b, y, x, ni, grad)
        real(dp), intent(in) :: b(:), y(:), x(:, :)
        integer, intent(in) :: ni(:)
        real(dp), intent(out) :: grad(:)
        real(dp), allocatable :: r(:), z(:), v2r(:)
        real(dp) :: alpha, sigma, s2, a2, den
        real(dp) :: tr_vinv, tr_vinv_j, sumz
        integer :: i, j0, n_i, p

        p = size(x, 2)
        alpha = b(p + 1)
        sigma = b(p + 2)
        s2 = sigma * sigma
        a2 = alpha * alpha
        grad = 0.0_dp
        j0 = 0
        if (s2 <= tiny(1.0_dp)) return

        do i = 1, size(ni)
            n_i = ni(i)
            allocate(r(n_i), z(n_i), v2r(n_i))
            r = y(j0 + 1:j0 + n_i) - matmul(x(j0 + 1:j0 + n_i, :), b(1:p))
            den = s2 + real(n_i, dp) * a2
            z = r / s2 - (a2 / (s2 * den)) * sum(r)
            grad(1:p) = grad(1:p) + matmul(transpose(x(j0 + 1:j0 + n_i, :)), z)

            sumz = sum(z)
            tr_vinv_j = real(n_i, dp) / den
            grad(p + 1) = grad(p + 1) - alpha * tr_vinv_j + alpha * sumz * sumz

            tr_vinv = real(n_i - 1, dp) / s2 + 1.0_dp / den
            v2r = z / s2 - (a2 / (s2 * den)) * sumz
            grad(p + 2) = grad(p + 2) - sigma * tr_vinv + sigma * dot_product(r, v2r)

            deallocate(r, z, v2r)
            j0 = j0 + n_i
        end do
    end subroutine grad_lmm

end module mla_lmm
