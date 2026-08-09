! SPDX-License-Identifier: GPL-2.0-or-later
module onls_linalg
    use onls_kinds, only : dp
    implicit none
    private
    public :: solve_spd, invert_spd
contains
    subroutine solve_spd(a, b, x, ok)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        logical, intent(out) :: ok
        real(dp), allocatable :: l(:,:), y(:)
        real(dp) :: s, scale, ridge
        integer :: n, i, j, k, attempt

        n = size(b)
        if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
            ok = .false.
            return
        end if
        allocate(l(n,n), y(n))
        scale = max(1.0_dp, maxval(abs(a)))
        ridge = 0.0_dp
        do attempt = 0, 10
            l = 0.0_dp
            ok = .true.
            do i = 1, n
                do j = 1, i
                    s = a(i,j)
                    if (i == j) s = s + ridge
                    do k = 1, j - 1
                        s = s - l(i,k) * l(j,k)
                    end do
                    if (i == j) then
                        if (s <= 100.0_dp * epsilon(1.0_dp) * scale) then
                            ok = .false.
                            exit
                        end if
                        l(i,j) = sqrt(s)
                    else
                        l(i,j) = s / l(j,j)
                    end if
                end do
                if (.not. ok) exit
            end do
            if (ok) exit
            if (abs(ridge) <= tiny(1.0_dp)) then
                ridge = 1.0e-12_dp * scale
            else
                ridge = 10.0_dp * ridge
            end if
        end do
        if (.not. ok) return
        do i = 1, n
            s = b(i)
            do k = 1, i - 1
                s = s - l(i,k) * y(k)
            end do
            y(i) = s / l(i,i)
        end do
        do i = n, 1, -1
            s = y(i)
            do k = i + 1, n
                s = s - l(k,i) * x(k)
            end do
            x(i) = s / l(i,i)
        end do
    end subroutine solve_spd

    subroutine invert_spd(a, ainv, ok)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        logical, intent(out) :: ok
        real(dp), allocatable :: e(:), col(:)
        integer :: n, j

        n = size(a,1)
        if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) then
            ok = .false.
            return
        end if
        allocate(e(n), col(n))
        ainv = 0.0_dp
        do j = 1, n
            e = 0.0_dp
            e(j) = 1.0_dp
            call solve_spd(a, e, col, ok)
            if (.not. ok) return
            ainv(:,j) = col
        end do
        ainv = 0.5_dp * (ainv + transpose(ainv))
    end subroutine invert_spd
end module onls_linalg
