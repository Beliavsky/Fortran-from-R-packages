! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/filtersimfast.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

!fast filtering algorithm used in simulation filter
subroutine filtersimfast(yt, ymiss, timevar, zt,tt, &
a1, ft,kt,finf, kinf, dt, jt, p, m, n,at)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m,n,dt,jt
    integer :: t
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(p,n) :: ft,finf
    real(dp), intent(in), dimension(m,p,n) :: kt,kinf
    real(dp), intent(inout), dimension(m,n + 1) :: at
    real(dp), dimension(p,n) :: vt
    real(dp) :: lik
    real(dp), external :: ddot

    external dgemv

    lik = 0.0_dp

    at(:,1) = a1
    if(dt > 0) then
        !diffuse filtering begins
        do t = 1, (dt - 1)
            at(:,t + 1) = at(:,t)
            call dfilter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
            tt(:,:,(t - 1) * timevar(3) + 1),at(:,t + 1),vt(:,t),ft(:,t),kt(:,:,t),&
            finf(:,t),kinf(:,:,t),p,m,p,lik)
        end do

        t = dt
        at(:,t + 1) = at(:,t)
        call dfilter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
        tt(:,:,(t - 1) * timevar(3) + 1),at(:,t + 1),vt(:,t),ft(:,t),kt(:,:,t),&
        finf(:,t),kinf(:,:,t),p,m,jt,lik)
        !non-diffuse filtering begins
        if(jt < p) then
            call filter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
            tt(:,:,(t - 1) * timevar(3) + 1),at(:,t + 1),vt(:,t),ft(:,t),kt(:,:,t),p,m,jt,lik)
        end if
    end if
    !Non-diffuse filtering continues from t=d+1, i=1
    do t = dt + 1, n
        at(:,t + 1) = at(:,t)
        call filter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
        tt(:,:,(t - 1) * timevar(3) + 1),at(:,t + 1),vt(:,t),ft(:,t),kt(:,:,t),p,m,0,lik)
    end do


end subroutine filtersimfast
