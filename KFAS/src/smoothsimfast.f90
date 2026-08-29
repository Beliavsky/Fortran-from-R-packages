! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/smoothsimfast.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

!fast disturbance smoother for simulation
subroutine smoothsimfast(yt, ymiss, timevar, zt, ht,tt, rtv,qt,a1, ft,kt,&
finf, kinf, dt, jt, p, m, n,r,epshat,etahat,rt0,rt1,needeps)
    use kfas_kinds, only: dp

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) :: p, m, r,n,dt,jt
    integer :: t, i
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p,p,(n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rtv
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(p,n) :: ft,finf
    real(dp), intent(in), dimension(m,p,n) :: kt,kinf
    real(dp), intent(inout), dimension(p,n) :: epshat
    real(dp), intent(inout), dimension(r,n) :: etahat
    real(dp), dimension(p,n) :: vt
    real(dp), dimension(m) :: at
    real(dp), dimension(m,m) :: im
    real(dp), intent(inout), dimension(m) :: rt0,rt1
    real(dp) :: lik
    lik = 0.0_dp
    at = a1
    if(dt > 0) then
        !diffuse filtering begins
        do t = 1, dt - 1
            call dfilter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
            tt(:,:,(t - 1) * timevar(3) + 1),at,vt(:,t),ft(:,t),kt(:,:,t),&
            finf(:,t),kinf(:,:,t),p,m,p,lik)
        end do

        t = dt
        call dfilter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
        tt(:,:,(t - 1) * timevar(3) + 1),at,vt(:,t),ft(:,t),kt(:,:,t),&
        finf(:,t),kinf(:,:,t),p,m,jt,lik)
        !non-diffuse filtering begins
        if(jt < p) then
            call filter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
            tt(:,:,(t - 1) * timevar(3) + 1),at,vt(:,t),ft(:,t),kt(:,:,t),p,m,jt,lik)
        end if
    end if
    !Non-diffuse filtering continues from t=d+1, i=1
    do t = dt + 1, n
        call filter1stepnv(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
        tt(:,:,(t - 1) * timevar(3) + 1),at,vt(:,t),ft(:,t),kt(:,:,t),p,m,0,lik)
    end do


    !smoothing begins

    im = 0.0_dp
    do i = 1, m
        im(i,i) = 1.0_dp
    end do

    rt0 = 0.0_dp

    do t = n, dt + 1, -1
        call smooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
        tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,1,rt0,etahat(:,t),epshat(:,t),needeps)
    end do

    if(dt > 0) then
        t = dt
        if(jt < p) then
            call smooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
            tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,jt + 1,rt0,etahat(:,t),epshat(:,t),needeps)
        end if
        rt1 = 0.0_dp
        call dsmooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
        tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,jt,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)
        do t = (dt - 1), 1, -1
            call dsmooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
            tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,p,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)
        end do
    end if

end subroutine smoothsimfast
