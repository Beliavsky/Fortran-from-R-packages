! SPDX-License-Identifier: GPL-2.0-or-later
! Translation of dirmult 0.1.3-5 by Torben Tvedebrink.
! See LICENSE and provenance/upstream/DESCRIPTION.

module dirmult_linalg
    use dirmult_types, only : dp
    implicit none
    private
    public :: solve_linear, invert_matrix

contains

    subroutine solve_linear(a, b, x, info)
        real(dp), intent(in) :: a(:,:), b(:)
        real(dp), intent(out) :: x(:)
        integer, intent(out) :: info
        real(dp), allocatable :: aa(:,:), bb(:), rowtmp(:)
        real(dp) :: factor, scale, tmp
        integer :: n, i, k, p, idx

        n = size(b)
        info = 0
        x = 0.0_dp
        if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
            info = -1
            return
        end if
        if (n == 0) return

        allocate(aa(n,n), bb(n), rowtmp(n))
        aa = a
        bb = b
        scale = max(1.0_dp, maxval(abs(aa)))

        do k = 1, n - 1
            idx = maxloc(abs(aa(k:n,k)), dim=1)
            p = k + idx - 1
            if (abs(aa(p,k)) <= 100.0_dp * epsilon(1.0_dp) * scale) then
                info = k
                return
            end if
            if (p /= k) then
                rowtmp = aa(k,:)
                aa(k,:) = aa(p,:)
                aa(p,:) = rowtmp
                tmp = bb(k)
                bb(k) = bb(p)
                bb(p) = tmp
            end if
            do i = k + 1, n
                factor = aa(i,k) / aa(k,k)
                aa(i,k) = 0.0_dp
                aa(i,k+1:n) = aa(i,k+1:n) - factor * aa(k,k+1:n)
                bb(i) = bb(i) - factor * bb(k)
            end do
        end do

        if (abs(aa(n,n)) <= 100.0_dp * epsilon(1.0_dp) * scale) then
            info = n
            return
        end if

        do i = n, 1, -1
            if (i < n) then
                x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n))) / aa(i,i)
            else
                x(i) = bb(i) / aa(i,i)
            end if
        end do
    end subroutine solve_linear

    subroutine invert_matrix(a, ainv, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(out) :: ainv(:,:)
        integer, intent(out) :: info
        real(dp), allocatable :: e(:), x(:)
        integer :: n, j, stat

        n = size(a,1)
        info = 0
        ainv = 0.0_dp
        if (size(a,2) /= n .or. size(ainv,1) /= n .or. size(ainv,2) /= n) then
            info = -1
            return
        end if
        allocate(e(n), x(n))
        do j = 1, n
            e = 0.0_dp
            e(j) = 1.0_dp
            call solve_linear(a, e, x, stat)
            if (stat /= 0) then
                info = stat
                return
            end if
            ainv(:,j) = x
        end do
    end subroutine invert_matrix

end module dirmult_linalg
