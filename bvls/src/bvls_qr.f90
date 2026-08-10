! SPDX-License-Identifier: GPL-2.0-or-later
!
! Householder least-squares kernel derived from the qr subroutine in
! Stark & Parker's BVLS source. The original source is preserved under
! original/bvls-master/src/bvls.f.
module bvls_qr
    use bvls_kinds, only : dp
    implicit none
    private
    public :: householder_lsq

contains

    subroutine householder_lsq(a, b, x, resq, info)
        real(dp), intent(in) :: a(:,:)
        real(dp), intent(in) :: b(:)
        real(dp), intent(out) :: x(:)
        real(dp), intent(out) :: resq
        integer, intent(out) :: info

        real(dp), allocatable :: aw(:,:), bw(:)
        real(dp) :: sq, qv1, u1, dotv, const, sumv
        integer :: m, n, i, j, jj

        m = size(a, 1)
        n = size(a, 2)
        x = 0.0_dp
        resq = -2.0_dp
        info = -2
        if (size(b) /= m .or. size(x) /= n .or. m < n) return

        allocate(aw(m,n), bw(m))
        aw = a
        bw = b

        resq = -1.0_dp
        info = -1
        do j = 1, n
            sq = sum(aw(j:m,j)**2)
            if (.not. (sq > 0.0_dp .or. sq < 0.0_dp)) return
            qv1 = -sign(sqrt(sq), aw(j,j))
            u1 = aw(j,j) - qv1
            aw(j,j) = qv1

            do jj = j + 1, n
                dotv = u1 * aw(j,jj)
                if (j < m) dotv = dotv + dot_product(aw(j+1:m,jj), aw(j+1:m,j))
                const = dotv / abs(qv1 * u1)
                if (j < m) aw(j+1:m,jj) = aw(j+1:m,jj) - const * aw(j+1:m,j)
                aw(j,jj) = aw(j,jj) - const * u1
            end do

            dotv = u1 * bw(j)
            if (j < m) dotv = dotv + dot_product(bw(j+1:m), aw(j+1:m,j))
            const = dotv / abs(qv1 * u1)
            bw(j) = bw(j) - const * u1
            if (j < m) bw(j+1:m) = bw(j+1:m) - const * aw(j+1:m,j)
        end do

        do i = n, 1, -1
            sumv = bw(i)
            if (i < n) sumv = sumv - dot_product(aw(i,i+1:n), x(i+1:n))
            if (.not. (aw(i,i) > 0.0_dp .or. aw(i,i) < 0.0_dp)) return
            x(i) = sumv / aw(i,i)
        end do

        resq = 0.0_dp
        if (n < m) resq = sum(bw(n+1:m)**2)
        info = 0
    end subroutine householder_lsq

end module bvls_qr
