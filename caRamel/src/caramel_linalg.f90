module caramel_linalg
    use caramel_kinds, only: dp
    implicit none
    private
    public :: determinant, solve_linear, cholesky_lower

contains

    real(dp) function determinant(a) result(det)
        real(dp), intent(in) :: a(:,:)
        real(dp), allocatable :: b(:,:)
        real(dp) :: factor, pivot_abs
        integer :: n, i, j, k, pivot

        n = size(a, 1)
        if (size(a, 2) /= n) error stop "determinant: matrix must be square"
        if (n == 0) then
            det = 1.0_dp
            return
        end if

        b = a
        det = 1.0_dp
        do k = 1, n
            pivot = k
            pivot_abs = abs(b(k,k))
            do i = k + 1, n
                if (abs(b(i,k)) > pivot_abs) then
                    pivot = i
                    pivot_abs = abs(b(i,k))
                end if
            end do
            if (pivot_abs <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(b)))) then
                det = 0.0_dp
                return
            end if
            if (pivot /= k) then
                do j = k, n
                    factor = b(k,j)
                    b(k,j) = b(pivot,j)
                    b(pivot,j) = factor
                end do
                det = -det
            end if
            det = det * b(k,k)
            do i = k + 1, n
                factor = b(i,k) / b(k,k)
                b(i,k:n) = b(i,k:n) - factor * b(k,k:n)
            end do
        end do
    end function determinant

    subroutine solve_linear(a, b, x, ok)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: aa(:,:), bb(:)
        real(dp) :: factor, tmp, pivot_abs, scale
        integer :: n, i, j, k, pivot

        n = size(a,1)
        if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) then
            error stop "solve_linear: inconsistent dimensions"
        end if
        aa = a
        bb = b
        scale = max(1.0_dp, maxval(abs(aa)))
        ok = .true.

        do k = 1, n
            pivot = k
            pivot_abs = abs(aa(k,k))
            do i = k + 1, n
                if (abs(aa(i,k)) > pivot_abs) then
                    pivot = i
                    pivot_abs = abs(aa(i,k))
                end if
            end do
            if (pivot_abs <= 100.0_dp * epsilon(1.0_dp) * scale) then
                ok = .false.
                x = 0.0_dp
                return
            end if
            if (pivot /= k) then
                do j = k, n
                    tmp = aa(k,j)
                    aa(k,j) = aa(pivot,j)
                    aa(pivot,j) = tmp
                end do
                tmp = bb(k)
                bb(k) = bb(pivot)
                bb(pivot) = tmp
            end if
            do i = k + 1, n
                factor = aa(i,k) / aa(k,k)
                aa(i,k:n) = aa(i,k:n) - factor * aa(k,k:n)
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do

        do i = n, 1, -1
            if (i < n) then
                x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n))) / aa(i,i)
            else
                x(i) = bb(i) / aa(i,i)
            end if
        end do
    end subroutine solve_linear

    subroutine cholesky_lower(a, l, ok)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: l(:,:)
        logical, intent(out) :: ok
        real(dp) :: s
        integer :: n, i, j, k

        n = size(a,1)
        if (size(a,2) /= n .or. any(shape(l) /= [n,n])) then
            error stop "cholesky_lower: inconsistent dimensions"
        end if
        l = 0.0_dp
        ok = .true.
        do i = 1, n
            do j = 1, i
                s = a(i,j)
                do k = 1, j - 1
                    s = s - l(i,k) * l(j,k)
                end do
                if (i == j) then
                    if (s <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(a(i,i)))) then
                        ok = .false.
                        return
                    end if
                    l(i,j) = sqrt(s)
                else
                    l(i,j) = s / l(j,j)
                end if
            end do
        end do
    end subroutine cholesky_lower

end module caramel_linalg
