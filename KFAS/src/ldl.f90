! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/ldl.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! LDL decomposition
subroutine ldl(a,n,tol,info)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: n
    integer, intent(inout) :: info
    integer :: i,j,k
    real(dp) :: di,tmp
    real(dp), intent(inout), dimension(n,n) :: a
    real(dp), intent(in) :: tol

    do i = 1, n
        di = a(i,i)
        if(abs(di) <= tol) then
            a(:,i) = 0.0_dp
        else
            do j = i + 1, n
                tmp = a(j,i) / di
                a(j,i) = tmp
                a(j,j) = a(j,j) - tmp**2 * di
                do k = j + 1, n
                    a(k,j) = a(k,j) - tmp * a(k,i)
                end do
            end do
        end if
    end do
    do i = 1,n
        a(i,(i + 1):n) = 0.0_dp
        if(a(i,i) < 0.0_dp) then
            info = -1
        end if
    end do

end subroutine ldl
