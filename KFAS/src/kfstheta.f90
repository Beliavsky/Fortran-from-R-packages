! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/kfstheta.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! signal smoothing algorithm for gaussian approximation algorithm
subroutine kfstheta(yt, ymiss, timevar, zt, ht,tt, rtv,qt,rqr, tv, a1, p1, p1inf, &
p, n, m, r,tol,rankp,thetahat,lik)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, n,r,tv
    integer, intent(inout) :: rankp
    integer :: t, i,d, j
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p,p,(n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rtv
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), dimension(m,m,(n - 1) * tv + 1) :: rqr
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), dimension(p,n) :: ft,finf
    real(dp), dimension(m,p,n) :: kt,kinf
    real(dp), dimension(p,n) :: vt
    real(dp), dimension(m) :: at,mhelp
    real(dp), dimension(m,m) :: pinf,pt
    real(dp), intent(in) :: tol
    real(dp), dimension(m) :: rt0,rt1
    real(dp), dimension(m,m) :: im
    real(dp), dimension(r,n) :: etahat
    real(dp) :: c
    real(dp), external :: ddot
    real(dp), intent(inout), dimension(n,p) :: thetahat
    real(dp), intent(inout) :: lik
    real(dp), dimension(p) :: epshat

    external dgemv, dsymv

    epshat = 0.0_dp
    lik = 0.0_dp
    c = 0.5_dp * log(8.0_dp * atan(1.0_dp))

    j = 0
    d = 0
    pt = p1
    pinf = p1inf
    at = a1
    if(rankp > 0) then
    !diffuse filtering
        diffuse: do while(d < n .and. rankp > 0)
            d = d + 1
            call dfilter1step(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),pinf,finf(:,d),kinf(:,:,d),rankp,lik,tol,c,p,m,j)
        end do diffuse

        if(rankp == 0 .and. j < p) then
            call filter1step(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at,pt,vt(:,d),ft(:,d),kt(:,:,d),lik,tol,c,p,m,j)

        else
            j = p
        end if

    end if

    !Non-diffuse filtering continues from t=d+1, i=1
    do t = d + 1, n
        call filter1step(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),ht(:,:,(t - 1) * timevar(2) + 1),&
        tt(:,:,(t - 1) * timevar(3) + 1),rqr(:,:,(t - 1) * tv + 1),&
        at,pt,vt(:,t),ft(:,t),kt(:,:,t),lik,tol,c,p,m,0)

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
        ft(:,t),kt(:,:,t), im,p,m,r,1,rt0,etahat(:,t),epshat,.FALSE.)
    end do

    if(d > 0) then
        t = d
        if(j < p) then
            call smooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
            tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,j + 1,rt0,etahat(:,t),epshat,.FALSE.)
        end if
        rt1 = 0.0_dp
        call dsmooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
        tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
        ft(:,t),kt(:,:,t), im,p,m,r,j,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat,.FALSE.)
        do t = (d - 1), 1, -1
            call dsmooth1step(ymiss(t,:), transpose(zt(:,:,(t - 1) * timevar(1) + 1)), ht(:,:,(t - 1) * timevar(2) + 1), &
            tt(:,:,(t - 1) * timevar(3) + 1), rtv(:,:,(t - 1) * timevar(4) + 1), qt(:,:,(t - 1) * timevar(5) + 1), vt(:,t), &
            ft(:,t),kt(:,:,t), im,p,m,r,p,rt0,rt1,finf(:,t),kinf(:,:,t),etahat(:,t),epshat,.FALSE.)
        end do
    end if

    at = a1

    call dsymv('l',m,1.0_dp,p1,m,rt0,1,1.0_dp,at,1)
    if(d > 0) then
        call dsymv('l',m,1.0_dp,p1inf,m,rt1,1,1.0_dp,at,1)
    end if
    call dgemv('n',p,m,1.0_dp,zt(:,:,1),p,at,1,0.0_dp,thetahat(1,:),1)

    do t = 2, n
        call dgemv('n',m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,at,1,0.0_dp,mhelp,1)
        at = mhelp
        call dgemv('n',m,r,1.0_dp,rtv(:,:,(t - 2) * timevar(4) + 1),m,etahat(:,t - 1),1,1.0_dp,at,1)
        call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,at,1,0.0_dp,thetahat(t,:),1)

    end do


end subroutine kfstheta
