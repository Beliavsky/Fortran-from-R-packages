! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/artransform.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

!transformation of unconstrained parameters to stationary region
subroutine artransform(phi,p)
    use kfas_kinds, only: dp
    implicit none

    integer, intent(in) :: p
    integer :: i, j
    real(dp), intent(inout), dimension(p) :: phi
    real(dp), dimension(p,p) :: u

    u = 0.0_dp
    do i = 1,p
        u(i,i) = phi(i)
    end do
    do i = 2, p
        do j = 1, i - 1
            u(i,j) = u(i - 1,j) - phi(i) * u(i - 1,i - j)
        end do
    end do
    phi = u(p,:)
end subroutine artransform
