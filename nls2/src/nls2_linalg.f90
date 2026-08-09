! SPDX-License-Identifier: GPL-2.0-only
module nls2_linalg
    use nls2_kinds, only : dp
    implicit none
    private
    public :: least_squares, invert_spd, norm2_sq

contains

    pure real(dp) function norm2_sq(x) result(v)
        real(dp), intent(in) :: x(:)
        v = dot_product(x, x)
    end function norm2_sq

    subroutine least_squares(a, b, x, rank_ok, rmat)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: rank_ok
        real(dp), intent(out), optional :: rmat(:,:)
        integer :: n, p, i, j
        real(dp), allocatable :: q(:,:), r(:,:), v(:), qt_b(:)
        real(dp) :: nv, tol

        n = size(a, 1)
        p = size(a, 2)
        x = 0.0_dp
        rank_ok = .false.
        if (size(b) /= n .or. size(x) /= p .or. n < p) return

        allocate(q(n,p), r(p,p), v(n), qt_b(p))
        q = 0.0_dp
        r = 0.0_dp
        tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(a))) * real(max(n,p), dp)

        do j = 1, p
            v = a(:,j)
            do i = 1, j - 1
                r(i,j) = dot_product(q(:,i), v)
                v = v - r(i,j) * q(:,i)
            end do
            ! One reorthogonalization pass is cheap and improves robustness.
            do i = 1, j - 1
                nv = dot_product(q(:,i), v)
                r(i,j) = r(i,j) + nv
                v = v - nv * q(:,i)
            end do
            r(j,j) = sqrt(max(0.0_dp, dot_product(v, v)))
            if (r(j,j) <= tol) return
            q(:,j) = v / r(j,j)
        end do

        qt_b = matmul(transpose(q), b)
        do i = p, 1, -1
            x(i) = qt_b(i)
            if (i < p) x(i) = x(i) - dot_product(r(i,i+1:p), x(i+1:p))
            x(i) = x(i) / r(i,i)
        end do
        rank_ok = .true.
        if (present(rmat)) rmat = r
    end subroutine least_squares

    subroutine invert_spd(a, ainv, ok)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        logical, intent(out) :: ok
        integer :: n, i, j, k
        real(dp), allocatable :: l(:,:), y(:), x(:)
        real(dp) :: s, tol

        n = size(a,1)
        ainv = 0.0_dp
        ok = .false.
        if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) return
        allocate(l(n,n), y(n), x(n))
        l = 0.0_dp
        tol = epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a))) * real(n, dp)

        do i = 1, n
            do j = 1, i
                s = a(i,j)
                do k = 1, j - 1
                    s = s - l(i,k) * l(j,k)
                end do
                if (i == j) then
                    if (s <= tol) return
                    l(i,j) = sqrt(s)
                else
                    l(i,j) = s / l(j,j)
                end if
            end do
        end do

        do j = 1, n
            y = 0.0_dp
            do i = 1, n
                s = merge(1.0_dp, 0.0_dp, i == j)
                if (i > 1) s = s - dot_product(l(i,1:i-1), y(1:i-1))
                y(i) = s / l(i,i)
            end do
            x = 0.0_dp
            do i = n, 1, -1
                s = y(i)
                if (i < n) s = s - dot_product(l(i+1:n,i), x(i+1:n))
                x(i) = s / l(i,i)
            end do
            ainv(:,j) = x
        end do
        ok = .true.
    end subroutine invert_spd

end module nls2_linalg
