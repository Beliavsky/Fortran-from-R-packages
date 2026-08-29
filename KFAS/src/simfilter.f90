! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/simfilter.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! simulation filter
subroutine simfilter(ymiss,timevar, yt, zt, ht, tt, rtv, qt, a1, p1, &
p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info,rankp,&
tol,sim,c,simwhat,simdim,antithetics)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, r, n, nsim,nnd,simdim,simwhat,antithetics
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: info,rankp
    integer :: t, i, d, j,k
    real(dp), intent(in) :: tol
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p,p,(n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rtv
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), intent(in), dimension(nsim) :: c
    real(dp), intent(inout), dimension(simdim,n,3 * nsim * antithetics + nsim) :: sim
    real(dp), intent(inout), dimension(p,n,nsim) :: epsplus
    real(dp), intent(inout), dimension(r,n,nsim) :: etaplus
    real(dp), intent(inout), dimension(m,nsim) :: aplus1

    real(dp), dimension(n,p) :: yplus
    real(dp), dimension(m,n + 1) :: aplus
    real(dp), dimension(p,n) :: ft,finf
    real(dp), dimension(m,p,n) :: kt,kinf
    real(dp), dimension(r,r,(n - 1) * timevar(5) + 1) :: cholqt
    real(dp), dimension(m,m) :: cholp1
    real(dp), dimension(r,r) :: rcholhelp
    real(dp), dimension(m,n,4) :: alphatmp

    real(dp), dimension(m,n + 1) :: at
    real(dp), dimension(m,n + 1) :: atplus
    real(dp), dimension(m,m,n + 1) :: pt,pinf
    real(dp), dimension(p,n) :: vt
    real(dp) :: lik
    real(dp), dimension(1,p) :: theta
    real(dp), dimension(p,p,1) :: thetavar

    real(dp), external :: ddot

    external kfilter, filtersimfast, ldl, dtrmv, dgemv


    at = 0.0_dp
    call kfilter(yt, ymiss, timevar, zt, ht,tt, rtv, qt, a1, p1, p1inf, &
    p,n,m,r,d,j, at, pt, vt, ft,kt, pinf, finf, kinf, lik, tol,rankp,theta,thetavar,0)

    do t = 1, (n - 1) * timevar(5) + 1
        if(r == 1) then
            cholqt(1,1,t) = sqrt(qt(1,1,t))
        else
            rcholhelp = qt(:,:,t)
            call ldl(rcholhelp,r,tol,info)
            if(info /= 0) then
                info = -2
                return
            end if
            do i = 1,r
                cholqt(i,i,t) = sqrt(rcholhelp(i,i))
            end do
            do i = 1,r - 1
                cholqt((i + 1):r,i,t) = rcholhelp((i + 1):r,i) * cholqt(i,i,t)
            end do
        end if
    end do


    if(nnd > 0) then
        if(m == 1) then
            cholp1(1,1) = sqrt(p1(1,1))
        else
            cholp1 = p1
            call ldl(cholp1,m,tol,info)
            if(info /= 0) then
                info = -3
                return
            end if
            do i = 1,m
                cholp1(i,i) = sqrt(cholp1(i,i))
            end do
            do i = 1,m - 1
                cholp1((i + 1):m,i) = cholp1((i + 1):m,i) * cholp1(i,i)
            end do
        end if
    end if

    do i = 1, nsim
        aplus = 0.0_dp
         aplus(:,1) = a1
        if(nnd > 0) then
            call dtrmv('l','n','n',m,cholp1,m,aplus1(:,i),1)
            aplus(:,1) = aplus(:,1) + aplus1(:,i)
        end if

        do t = 1, n
            do k = 1, p
                if(ymiss(t,k) == 0) then
                    yplus(t,k) = epsplus(k,t,i) * sqrt(ht(k,k,(t - 1) * timevar(2) + 1)) + &
                    ddot(m,zt(k,:,(t - 1) * timevar(1) + 1),1,aplus(:,t),1)
                end if
            end do
            call dtrmv('l','n','n',r,cholqt(:,:,(t - 1) * timevar(5) + 1),r,etaplus(:,t,i),1)
            call dgemv('n',m,m,1.0_dp,tt(:,:,(t - 1) * timevar(3) + 1),m,aplus(:,t),1,0.0_dp,aplus(:,t + 1),1)
            call dgemv('n',m,r,1.0_dp,rtv(:,:,(t - 1) * timevar(4) + 1),m,etaplus(:,t,i),1,1.0_dp,aplus(:,t + 1),1)
        end do


        atplus = 0.0_dp
        call filtersimfast(yplus, ymiss, timevar, zt,tt, a1, ft,kt,&
        finf, kinf, d, j, p, m, n,atplus)
        if(simwhat == 4) then
            do t = 1, n
                sim(:,t,i) = at(:,t) - atplus(:,t) + aplus(:,t)
                if(antithetics == 1) then
                    sim(:,t,i + nsim) = at(:,t) + atplus(:,t) - aplus(:,t)
                    sim(:,t,i + 2 * nsim) = at(:,t) + c(i) * (sim(:,t,i) - at(:,t))
                    sim(:,t,i + 3 * nsim) = at(:,t) + c(i) * (sim(:,t,i + nsim) - at(:,t))
                end if
            end do
        else
            do t = 1, n
                alphatmp(:,t,1) = at(:,t) - atplus(:,t) + aplus(:,t)
                call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,&
                alphatmp(:,t,1),1,0.0_dp,sim(:,t,i),1)
                if(antithetics == 1) then
                    alphatmp(:,t,2) = at(:,t) + atplus(:,t) - aplus(:,t)
                    alphatmp(:,t,3) = at(:,t) + c(i) * (alphatmp(:,t,1) - at(:,t))
                    alphatmp(:,t,4) = at(:,t) + c(i) * (alphatmp(:,t,2) - at(:,t))
                    call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,&
                    alphatmp(:,t,2),1,0.0_dp,sim(:,t,i + nsim),1)
                    call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,&
                    alphatmp(:,t,3),1,0.0_dp,sim(:,t,i + 2 * nsim),1)
                    call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,&
                    alphatmp(:,t,4),1,0.0_dp,sim(:,t,i + 3 * nsim),1)
                end if
            end do
        end if
    end do
end subroutine simfilter

