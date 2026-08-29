! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/kfilter2.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

  ! Subroutine for Kalman filtering of linear gaussian state space model

subroutine kfilter2(yt, ymiss, timevar, zt, ht,tt, rt, qt, a1, p1, p1inf, p,n,m,r,d,j,&
at, pt, vt, ft,kt, pinf, finf, kinf, lik, tol,rankp,theta,thetavar,filtersignal, att, ptt)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, r, n,filtersignal
    integer, intent(inout) :: d, j, rankp
    integer :: t,tv
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p,p,(n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rt
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), intent(in) :: tol
    real(dp), intent(inout), dimension(m,n + 1) :: at
    real(dp), intent(inout), dimension(m,m,n + 1) :: pt,pinf
    real(dp), intent(inout), dimension(p,n) :: vt,ft,finf
    real(dp), intent(inout), dimension(m,p,n) :: kt,kinf
    real(dp), intent(inout) :: lik
    real(dp), intent(inout), dimension(p,p,n) :: thetavar
    real(dp), intent(inout), dimension(n,p) :: theta
    real(dp), intent(inout), dimension(m,n) :: att
    real(dp), intent(inout), dimension(m,m,n) :: ptt
    real(dp), dimension(m,r) :: mr
    real(dp), dimension(p,m) :: pm
    real(dp) :: c
    real(dp), external :: ddot
    real(dp), dimension(m,m,(n - 1) * max(timevar(4),timevar(5)) + 1) :: rqr
    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2

    c = 0.5_dp * log(8.0_dp * atan(1.0_dp))

    lik = 0.0_dp
    tv = max(timevar(4),timevar(5))
    do t = 1, (n - 1) * tv + 1
        call dsymm('r','l',m,r,1.0_dp,qt(:,:,(t - 1) * timevar(5) + 1),r,rt(:,:,(t - 1) * timevar(4) + 1),m,0.0_dp,mr,m)
        call dgemm('n','t',m,m,r,1.0_dp,mr,m,rt(:,:,(t - 1) * timevar(4) + 1),m,0.0_dp,rqr(:,:,t),m)
    end do

    j = 0
    d = 0
    pinf(:,:,1) = p1inf
    pt(:,:,1) = p1
    at(:,1) = a1
    ! diffuse initialization
    if(rankp > 0) then
        diffuse: do while(d < n .and. rankp > 0)
            d = d + 1
            at(:,d + 1) = at(:,d)
            pt(:,:,d + 1) = pt(:,:,d)
            pinf(:,:,d + 1) = pinf(:,:,d)
            call dfilter1step2(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at(:,d + 1),pt(:,:,d + 1),vt(:,d),ft(:,d),kt(:,:,d),pinf(:,:,d + 1),finf(:,d),kinf(:,:,d),rankp,lik,tol,c,p,m,j, &
            att(:,d),ptt(:,:,d))
        end do diffuse


        if(rankp == 0 .and. j < p) then
           !non-diffuse filtering begins
            call filter1step2(ymiss(d,:),yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at(:,d + 1),pt(:,:,d + 1),vt(:,d),ft(:,d),kt(:,:,d),lik,tol,c,p,m,j, att(:, d), ptt(:,:, d))
        else
            j = p
        end if
    end if

    !Non-diffuse filtering continues from t=d+1, i=1

    do t = d + 1, n
        at(:,t + 1) = at(:,t)
        pt(:,:,t + 1) = pt(:,:,t)
        call filter1step2(ymiss(t,:),yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),ht(:,:,(t - 1) * timevar(2) + 1),&
        tt(:,:,(t - 1) * timevar(3) + 1),rqr(:,:,(t - 1) * tv + 1),&
        at(:,t + 1),pt(:,:,t + 1),vt(:,t),ft(:,t),kt(:,:,t),lik,tol,c,p,m,0, att(:,t), ptt(:,:,t))
    end do

    if(filtersignal == 1) then
        do t = 1, n
            call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,at(:,t),1,0.0_dp,theta(t,:),1)
            call dsymm('r','u',p,m,1.0_dp,pt(:,:,t),m,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,pm,p)
            call dgemm('n','t',p,p,m,1.0_dp,pm,p,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,thetavar(:,:,t),p)
        end do
    end if
end subroutine kfilter2
