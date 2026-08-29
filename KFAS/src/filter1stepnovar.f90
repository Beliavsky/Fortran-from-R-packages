! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/filter1stepnovar.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

!diffuse filtering for single time point
subroutine dfilter1stepnv(ymiss, yt, zt, tt, at, vt, ft,kt,&
finf,kinf,p,m,j,lik)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, j
    integer :: i
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(p) :: yt
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), intent(inout), dimension(m) :: at
    real(dp), intent(inout), dimension(p) :: vt
    real(dp), intent(in), dimension(p) :: ft,finf
    real(dp), intent(in), dimension(m,p) :: kt,kinf
    real(dp), intent(inout) :: lik
    real(dp), dimension(m) :: ahelp

    real(dp), external :: ddot
    external dgemv

    do i = 1, j
        if(ymiss(i) == 0) then
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (finf(i) > 0.0_dp) then
                at = at + vt(i) / finf(i) * kinf(:,i)
                lik = lik - 0.5_dp * log(finf(i))
            else
                if (ft(i) > 0.0_dp) then
                    at = at + vt(i) / ft(i) * kt(:,i)
                    lik = lik - 0.5_dp * (log(ft(i)) + vt(i)**2 / ft(i))
                end if
            end if
        end if
    end do
    if(j == p) then
        call dgemv('n',m,m,1.0_dp,tt,m,at,1,0.0_dp,ahelp,1)
        at = ahelp
    end if
end subroutine dfilter1stepnv


!non-diffuse filtering for single time point
subroutine filter1stepnv(ymiss, yt, zt, tt, at,vt, ft,kt,p,m,j,lik)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m,j
    integer :: i
    integer, intent(in), dimension(p) :: ymiss
    real(dp), intent(in), dimension(p) :: yt
    real(dp), intent(in), dimension(m,p) :: zt
    real(dp), intent(in), dimension(m,m) :: tt
    real(dp), intent(inout), dimension(m) :: at
    real(dp), intent(in), dimension(p) :: ft
    real(dp), intent(inout), dimension(p) :: vt
    real(dp), intent(in), dimension(m,p) :: kt
    real(dp), intent(inout) :: lik
    real(dp), dimension(m) :: ahelp

    real(dp), external :: ddot

    external dgemv

    do i = j + 1, p
        if(ymiss(i) == 0) then
            vt(i) = yt(i) - ddot(m,zt(:,i),1,at,1)
            if (ft(i) > 0.0_dp) then
                at = at + vt(i) / ft(i) * kt(:,i)
                lik = lik - 0.5_dp * (log(ft(i)) + vt(i)**2 / ft(i))
            end if
        end if
    end do

    call dgemv('n',m,m,1.0_dp,tt,m,at,1,0.0_dp,ahelp,1)
    at = ahelp

end subroutine filter1stepnv
