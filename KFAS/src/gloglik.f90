! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/gloglik.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Subroutine for computing the log-Likelihood of general linear gaussian state space model

subroutine gloglik(yt, ymiss, timevar, zt, ht, tt, rt, qt, a1, p1, p1inf,&
p, m, r, n, lik, tol,rankp)
    use kfas_kinds, only: dp


    implicit none

    integer, intent(in) :: p, m, r, n
    integer, intent(inout) :: rankp
    integer :: t,d,j,tv
    integer, intent(in), dimension(p,n) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(p,n) :: yt
    real(dp), intent(in), dimension(m,p,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p,p,(n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rt
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), intent(in) :: tol
    real(dp), intent(inout) :: lik
    real(dp), dimension(m) :: at
    real(dp), dimension(p) :: vt,ft,finf
    real(dp), dimension(m,p) :: kt,kinf
    real(dp), dimension(m,m) :: pt,pinf
    real(dp), dimension(m,r) :: mr
    real(dp) :: c
    real(dp), external :: ddot
    real(dp), dimension(m,m,(n - 1) * max(timevar(4),timevar(5)) + 1) :: rqr

    external dgemm, dsymm, dgemv, dsymv, dsyr, dsyr2, marginalxx


    ! compute RQR'
    tv = max(timevar(4),timevar(5))
    do t = 1, (n - 1) * tv + 1
        call dsymm('r','l',m,r,1.0_dp,qt(:,:,(t - 1) * timevar(5) + 1),r,rt(:,:,(t - 1) * timevar(4) + 1),m,0.0_dp,mr,m)
        call dgemm('n','t',m,m,r,1.0_dp,mr,m,rt(:,:,(t - 1) * timevar(4) + 1),m,0.0_dp,rqr(:,:,t),m)
    end do


    ! constant term for log-likelihood
    c = 0.5_dp * log(8.0_dp * atan(1.0_dp))
    lik = 0.0_dp

    j = 0
    d = 0
    pt = p1
    at = a1
    pinf = p1inf
    ! Diffuse initialization
    if(rankp > 0) then
        diffuse: do while(d < n .and. rankp > 0)
            d = d + 1

            call dfilter1step(ymiss(:,d),yt(:,d),zt(:,:,(d - 1) * timevar(1) + 1),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1), at,pt,vt,ft,kt,pinf,finf,kinf,rankp,lik,tol,c,p,m,j)

        end do diffuse

        if(rankp == 0 .and. j < p) then
            !non-diffuse filtering begins
            call filter1step(ymiss(:,d),yt(:,d),zt(:,:,(d - 1) * timevar(1) + 1),ht(:,:,(d - 1) * timevar(2) + 1),&
            tt(:,:,(d - 1) * timevar(3) + 1),rqr(:,:,(d - 1) * tv + 1),at,pt,vt,ft,kt,lik,tol,c,p,m,j)

        else
            j = p

        end if
    end if



    !Non-diffuse filtering continues from t=d+1, i=1

    do t = d + 1, n
        call filter1step(ymiss(:,t),yt(:,t),zt(:,:,(t - 1) * timevar(1) + 1),ht(:,:,(t - 1) * timevar(2) + 1),&
        tt(:,:,(t - 1) * timevar(3) + 1),rqr(:,:,(t - 1) * tv + 1),at,pt,vt,ft,kt,lik,tol,c,p,m,0)

    end do

end subroutine gloglik
