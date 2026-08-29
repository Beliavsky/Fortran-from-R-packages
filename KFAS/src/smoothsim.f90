! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/smoothsim.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! disturbance smoothing algorithm for simulation
subroutine smoothsim(yt, ymiss, timevar, zt, ht,tt, rtv,qt,rqr, a1, p1, p1inf, &
d, j, p, m, n, r,tol,rankp,ft,finf,kt,kinf,epshat,etahat,rt0,rt1,needeps)
    use kfas_kinds, only: dp

    implicit none

    logical, intent(in) :: needeps
    integer, intent(in) :: p, m, n,r
    integer, intent(inout) :: d, j,rankp
    integer :: t, i,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p,p,(n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rtv
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), dimension(m,m,(n - 1) * max(timevar(4),timevar(5)) + 1) :: rqr
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), intent(inout), dimension(p,n) :: ft,finf
    real(dp), intent(inout), dimension(m,p,n) :: kt,kinf
    real(dp), dimension(p,n) :: vt
    real(dp), dimension(m) :: at
    real(dp), dimension(m,m) :: pt,pinf
    real(dp), intent(in) :: tol
    real(dp), dimension(m,m) :: im
    real(dp), intent(inout), dimension(r,n) :: etahat
    real(dp), intent(inout), dimension(p,n) :: epshat
    real(dp), intent(inout), dimension(m) :: rt0,rt1
    real(dp) :: lik
    real(dp), external :: ddot

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2, dger

    tv = max(timevar(4),timevar(5))
    lik = 0.0_dp

    j = 0
    d = 0
    pinf = p1inf
    pt = p1
    at = a1
    if(rankp > 0) then
        diffuse: do while(d < n .and. rankp > 0)
            d = d + 1
            call dfilter1step(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),pinf,finf(:,d),kinf(:,:,d),rankp,lik,tol,0.0_dp,p,m,j)
        end do diffuse
        if(rankp == 0 .and. j < p) then
            call filter1step(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),lik,tol,0.0_dp,p,m,j)
        else
            j = p
        end if
    end if

    !Non-diffuse filtering continues from t=d+1, i=1

    do t = d + 1, n
        call filter1step(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),ht(:,:,(t - 1) * timevar(2) + 1),&
        tt(:,:,(t - 1) * timevar(3) + 1),rqr(:,:,(t - 1) * tv + 1),&
        at,pt,vt(:,t),ft(:,t),kt(:,:,t),lik,tol,0.0_dp,p,m,0)
    end do

    !smoothing begins

    im = 0.0_dp
    do i = 1, m
        im(i,i) = 1.0_dp
    end do

    rt0 = 0.0_dp

    do t = n, d + 1, -1
        call smooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
        tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,1,rt0,etahat(:,t),epshat(:,t),needeps)
    end do

    if(d > 0) then
        t = d
        if(j < p) then
            call smooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
            tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,j + 1,rt0,etahat(:,t),epshat(:,t),needeps)
        end if
        rt1 = 0.0_dp
        call dsmooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
        tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,j,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)
        do t = (d - 1), 1, -1
            call dsmooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
            tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,p,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat(:,t),needeps)
        end do
    end if

end subroutine smoothsim
