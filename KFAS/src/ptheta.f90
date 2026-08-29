! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/ptheta.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! functions for computing p(theta)

! Only differences with gloglik is what is stored/returned (ft,finf,kt,kinf)
! Also no missing values as yt is actually the signals, and no ht
subroutine pthetafirst(yt, timevar, zt, tt, rqr, a1, p1, p1inf,&
p, m, n, lik, tol,rankp2,kt,kinf,ft,finf,d,j)
    use kfas_kinds, only: dp


    implicit none

    integer, intent(in) :: p, m, n
    integer, intent(inout) :: rankp2,d,j
    integer :: t, tv,rankp
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), intent(in) :: tol
    real(dp), intent(inout) :: lik
    real(dp), dimension(m) :: at
    real(dp), dimension(p) :: vt
    real(dp), intent(inout), dimension(p,n) :: ft,finf
    real(dp), intent(inout), dimension(m,p,n) :: kt,kinf
    real(dp), dimension(m,m) :: pt,pinf
    real(dp), intent(inout), dimension(m,m,(n - 1) * max(timevar(4),timevar(5)) + 1) :: rqr
    integer, dimension(p) :: ymiss
    real(dp), dimension(p,p) :: ht

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2

    tv = max(timevar(4),timevar(5))
    ymiss = 0
    ht = 0.0_dp

    rankp = rankp2
    j = 0
    d = 0
    at = a1
    pt = p1
    pinf = p1inf

    ! Diffuse initialization
    if(rankp > 0) then
        diffuse: do while(d < n .and. rankp > 0)
            d = d + 1
            call dfilter1step(ymiss,yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht,&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at,pt,vt,ft(:,d),kt(:,:,d),pinf,finf(:,d),kinf(:,:,d),rankp,lik,tol,0.0_dp,p,m,j)
        end do diffuse
        if(rankp == 0 .and. j < p) then
            call filter1step(ymiss,yt(d,:),transpose(zt(:,:,(d - 1) * timevar(1) + 1)),ht,&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),&
            at,pt,vt,ft(:,d),kt(:,:,d),lik,tol,0.0_dp,p,m,j)
        else
            j = p
        end if
    end if

    !Non-diffuse filtering continues from t=d+1, i=1

    do t = d + 1, n
        call filter1step(ymiss,yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),ht,&
        tt(:,:,(t - 1) * timevar(3) + 1),rqr(:,:,(t - 1) * tv + 1),&
        at,pt,vt,ft(:,t),kt(:,:,t),lik,tol,0.0_dp,p,m,0)
    end do




end subroutine pthetafirst


! use output of pthetafirst, kt and ft do not change
subroutine pthetarest(yt, timevar, zt, tt, a1,&
p, m, n, lik, kt,kinf,ft,finf,dt,jt)
    use kfas_kinds, only: dp


    implicit none

    integer, intent(in) :: p, m, n,dt,jt
    integer :: t
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(inout) :: lik
    real(dp), dimension(m) :: at
    real(dp), dimension(p) :: vt
    real(dp), intent(in), dimension(p,n) :: ft,finf
    real(dp), intent(in), dimension(m,p,n) :: kt,kinf
    integer, dimension(p) :: ymiss
    real(dp), external :: ddot

    external dgemv
    ymiss = 0
    at = a1
    if(dt > 0) then
        !diffuse filtering begins
        do t = 1, dt - 1
            call dfilter1stepnv(ymiss,yt(t,:),&
            transpose(zt(:,:,(t - 1) * timevar(1) + 1)),tt(:,:,(t - 1) * timevar(3) + 1),&
            at,vt,ft(:,t),kt(:,:,t), finf(:,t),kinf(:,:,t),p,m,p,lik)
        end do

        t = dt
        call dfilter1stepnv(ymiss,yt(t,:),&
        transpose(zt(:,:,(t - 1) * timevar(1) + 1)),tt(:,:,(t - 1) * timevar(3) + 1),&
        at,vt,ft(:,t),kt(:,:,t),finf(:,t),kinf(:,:,t),p,m,jt,lik)
        !non-diffuse filtering begins
        if(jt < p) then
            call filter1stepnv(ymiss,yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
            tt(:,:,(t - 1) * timevar(3) + 1),at,vt,ft(:,t),kt(:,:,t),p,m,jt,lik)
        end if
    end if
    !Non-diffuse filtering continues from t=d+1, i=1
    do t = dt + 1, n
        call filter1stepnv(ymiss,yt(t,:),transpose(zt(:,:,(t - 1) * timevar(1) + 1)),&
        tt(:,:,(t - 1) * timevar(3) + 1),at,vt,ft(:,t),kt(:,:,t),p,m,0,lik)
    end do

end subroutine pthetarest

