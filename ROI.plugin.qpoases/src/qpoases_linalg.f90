! SPDX-License-Identifier: LGPL-2.1-or-later
module qpoases_linalg
    use qpoases_kinds, only : dp
    implicit none
    private
    public :: solve_linear, norm_inf, norm2_vec, symmetric_part, matrix_rank

contains

    pure real(dp) function norm_inf(x) result(v)
        real(dp), intent(in) :: x(:)
        if (size(x) == 0) then
            v = 0.0_dp
        else
            v = maxval(abs(x))
        end if
    end function norm_inf

    pure real(dp) function norm2_vec(x) result(v)
        real(dp), intent(in) :: x(:)
        v = sqrt(max(0.0_dp, dot_product(x, x)))
    end function norm2_vec

    pure subroutine symmetric_part(a, s)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: s(size(a,1), size(a,2))
        s = 0.5_dp * (a + transpose(a))
    end subroutine symmetric_part

    subroutine solve_linear(a, b, x, info)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: m(:,:), rhs(:), tmp_row(:)
        real(dp) :: pivot, factor, scale
        integer :: n, i, k, p

        n = size(b)
        info = 0
        if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
            info = -1
            return
        end if
        if (n == 0) return

        allocate(m(n,n), rhs(n), tmp_row(n))
        m = a
        rhs = b
        scale = max(1.0_dp, maxval(abs(m)))

        do k = 1, n - 1
            p = k - 1 + maxloc(abs(m(k:n,k)), dim=1)
            pivot = abs(m(p,k))
            if (pivot <= 100.0_dp * epsilon(1.0_dp) * scale) then
                info = k
                return
            end if
            if (p /= k) then
                tmp_row = m(k,:)
                m(k,:) = m(p,:)
                m(p,:) = tmp_row
                pivot = rhs(k)
                rhs(k) = rhs(p)
                rhs(p) = pivot
            end if
            do i = k + 1, n
                factor = m(i,k) / m(k,k)
                m(i,k) = 0.0_dp
                m(i,k+1:n) = m(i,k+1:n) - factor * m(k,k+1:n)
                rhs(i) = rhs(i) - factor * rhs(k)
            end do
        end do

        if (abs(m(n,n)) <= 100.0_dp * epsilon(1.0_dp) * scale) then
            info = n
            return
        end if

        x(n) = rhs(n) / m(n,n)
        do i = n - 1, 1, -1
            x(i) = (rhs(i) - dot_product(m(i,i+1:n), x(i+1:n))) / m(i,i)
        end do
    end subroutine solve_linear

    integer function matrix_rank(a, tol) result(r)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in), optional :: tol
        real(dp), allocatable :: m(:,:), tmp(:)
        real(dp) :: t, pivot, factor, scale
        integer :: nr, nc, row, col, p, i

        nr = size(a,1)
        nc = size(a,2)
        allocate(m(nr,nc), tmp(nc))
        m = a
        scale = max(1.0_dp, maxval(abs(m)))
        t = sqrt(epsilon(1.0_dp)) * scale
        if (present(tol)) t = tol
        r = 0
        row = 1
        do col = 1, nc
            if (row > nr) exit
            p = row - 1 + maxloc(abs(m(row:nr,col)), dim=1)
            pivot = abs(m(p,col))
            if (pivot <= t) cycle
            if (p /= row) then
                tmp = m(row,:)
                m(row,:) = m(p,:)
                m(p,:) = tmp
            end if
            do i = row + 1, nr
                factor = m(i,col) / m(row,col)
                m(i,col:nc) = m(i,col:nc) - factor * m(row,col:nc)
            end do
            r = r + 1
            row = row + 1
        end do
    end function matrix_rank
end module qpoases_linalg
