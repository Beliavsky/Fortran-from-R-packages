! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/isamplefilter.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Importance sampling filtering of non-gaussian model

subroutine isamplefilter(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, dist, &
p, n, m, r, theta, maxiter,rankp,convtol, nnd,nsim,epsplus,etaplus,&
aplus1,c,tol,info,antithetics,w,sim,simwhat,simdim, expected, htol)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p,m, r, n,nnd,antithetics,nsim&
    ,simwhat,simdim,rankp, expected
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: info, maxiter
    integer :: t, j,i,k,maxiter2,maxitermax,info2
    real(dp), intent(in) :: convtol,tol
    real(dp), intent(inout) :: htol
    real(dp), intent(in), dimension(n,p) :: u
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rtv
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), intent(in),dimension(nsim) :: c
    real(dp), intent(inout), dimension(p,n,nsim) :: epsplus
    real(dp), intent(inout), dimension(r,n,nsim) :: etaplus
    real(dp), intent(inout), dimension(m,nsim) :: aplus1
    real(dp), intent(inout), dimension(n,p) :: theta
    real(dp), dimension(p,p,n) :: ht
    real(dp), intent(inout), dimension(simdim,n,3 * nsim * antithetics + nsim) :: sim
    real(dp), intent(inout), dimension(n,3 * nsim * antithetics + nsim) :: w
    real(dp), dimension(n,p) :: ytilde
    real(dp), dimension(n) :: tmp
    real(dp), external :: ddot

    integer, dimension(n,p) :: ymiss2
    real(dp), dimension(p,n,nsim) :: epsplus2
    real(dp), dimension(r,n,nsim) :: etaplus2
    real(dp), dimension(m,nsim) :: aplus12
    real(dp) :: diff
    real(dp) :: lik
    real(dp), dimension(:,:), allocatable :: tsim
    real(dp), dimension(:,:,:), allocatable :: sim2
    external approx, simgaussian
    allocate(sim2(simdim,n,3 * nsim * antithetics + nsim))

    ht = 0.0_dp
    ytilde = 0.0_dp
    w = 1.0_dp

    epsplus2 = epsplus
    etaplus2 = etaplus
    aplus12 = aplus1
    ymiss2 = ymiss
    ymiss2(1,:) = 1

    call simgaussian(ymiss2(1,:),timevar, ytilde(1,:), zt(:,:,1), &
    ht(:,:,1), tt(:,:,1), rtv(:,:,1), &
    qt(:,:,1), a1, p1, p1inf, nnd,nsim, epsplus2(:,1,:), etaplus2(:,1,:), aplus12(:,:), &
    p, 1, m, r, info,rankp,tol,sim(:,1,:),c,simwhat,simdim,antithetics)


    if(info /= 0) then
        return
    end if
    maxitermax = 0

    do i = 1, (n - 1) ! increase time

        ht = 0.0_dp
        ytilde = 0.0_dp

        maxiter2 = maxiter
        info2 = 0
        ! approximate
        call approx(yt(1:i,:), ymiss(1:i,:), timevar, zt(:,:,1:((i - 1) * timevar(1) + 1)), &
        tt(:,:,1:((i - 1) * timevar(3) + 1)), rtv(:,:,1:((i - 1) * timevar(4) + 1)), ht(:,:,1:i),&
        qt(:,:,1:((i - 1) * timevar(5) + 1)), a1, p1,p1inf, p,i,m,r,&
        theta(1:i,:), u(1:i,:), ytilde(1:i,:), dist,maxiter2,tol,rankp,convtol,diff,lik,&
        info2, expected, htol)

        if(info2 /= 0 .and. info2 /= 3) then !check for errors in approximating algorithm
            info = info2
            return
        end if

        if(maxiter2 > maxitermax) then
            maxitermax = maxiter2
        end if
        epsplus2 = epsplus
        etaplus2 = etaplus
        aplus12 = aplus1
        ymiss2 = ymiss
        ymiss2(i + 1,:) = 1
        sim2 = 0.0_dp
        info2 = 0
        ! simulate signals
        call simgaussian(ymiss2(1:(i + 1),:),timevar, ytilde(1:(i + 1),:), zt(:,:,1:(i * timevar(1) + 1)), &
        ht(:,:,1:(i + 1)), tt(:,:,1:(i * timevar(3) + 1)), rtv(:,:,1:(i * timevar(4) + 1)), &
        qt(:,:,1:(i * timevar(5) + 1)), a1, p1, p1inf, nnd,nsim, epsplus2(:,1:(i + 1),:), &
        etaplus2(:,1:(i + 1),:), aplus12(:,:),p, i + 1, m, r, info2,rankp,tol,&
        sim2(:,1:(i + 1),:),c,simwhat,simdim,antithetics)

        if(info2 /= 0) then
            info = info2
            return
        end if

        ! compute importance weights

        if(simwhat == 5) then
            do j = 1,p
                select case(dist(j))
                    case(2) !poisson
                        tmp(1:i) = exp(theta(1:i,j))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then
                                w(i + 1,:) = w(i + 1,:) * exp(yt(t,j) * (sim2(j,t,:) - theta(t,j))-&
                                u(t,j) * (exp(sim2(j,t,:)) - tmp(t))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim2(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                    case(3) !binomial
                        tmp(1:i) = log(1.0_dp + exp(theta(1:i,j)))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then

                                w(i + 1,:) = w(i + 1,:) * exp(yt(t,j) * (sim2(j,t,:) - theta(t,j))-&
                                u(t,j) * (log(1.0_dp + exp(sim2(j,t,:))) - tmp(t))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim2(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                            end if
                        end do
                    case(4) ! gamma
                        tmp(1:i) = exp(-theta(1:i,j))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then
                                w(i + 1,:) = w(i + 1,:) * exp(u(t,j) * (yt(t,j) * (tmp(t) - exp(-sim2(j,t,:)))&
                                +theta(t,j) - sim2(j,t,:))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim2(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                    case(5) !negbin
                        tmp(1:i) = exp(theta(1:i,j))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then
                                w(i + 1,:) = w(i + 1,:) * exp(yt(t,j) * (sim2(j,t,:) - theta(t,j)) +&
                                (yt(t,j) + u(t,j)) * log((u(t,j) + tmp(t)) / (u(t,j) + exp(sim2(j,t,:))))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim2(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                end select
            end do

        else
            allocate(tsim(p,(3 * nsim * antithetics + nsim) * (5 - simwhat)))
            do j = 1,p
                select case(dist(j))
                    case(2) !poisson
                        tmp(1:i) = exp(theta(1:i,j))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then
                                do k = 1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim2(:,t,k),1)
                                end do
                                w(i + 1,:) = w(i + 1,:) * exp(yt(t,j) * (tsim(j,:) - theta(t,j))-&
                                u(t,j) * (exp(tsim(j,:)) - tmp(t))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                            end if
                        end do
                    case(3) !binomial
                        tmp(1:i) = log(1.0_dp + exp(theta(1:i,j)))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then
                                do k = 1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim2(:,t,k),1)
                                end do
                                w(i + 1,:) = w(i + 1,:) * exp(yt(t,j) * (tsim(j,:) - theta(t,j))-&
                                u(t,j) * (log(1.0_dp + exp(tsim(j,:))) - tmp(t))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                            end if
                        end do
                    case(4) ! gamma
                        tmp(1:i) = exp(-theta(1:i,j))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then
                                do k = 1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim2(:,t,k),1)
                                end do
                                w(i + 1,:) = w(i + 1,:) * exp(u(t,j) * (yt(t,j) * (tmp(t) - exp(-tsim(j, &
                                    &:))) + theta(t,j) - tsim(j,:))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                    case(5) !negbin
                        tmp(1:i) = exp(theta(1:i,j))
                        do t = 1,i
                            if(ymiss2(t,j) == 0) then
                                do k = 1,3 * nsim * antithetics + nsim
                                    tsim(j,k) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim2(:,t,k),1)
                                end do
                                w(i + 1,:) = w(i + 1,:) * exp(yt(t,j) * (tsim(j,:) - theta(t,j)) +&
                                (yt(t,j) + u(t,j)) * log((u(t,j) + tmp(t)) / (u(t,j) + exp(tsim(j,:))))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                end select
            end do
            deallocate(tsim)
        end if
        sim(:,i + 1,:) = sim2(:,i + 1,:)
    end do
    maxiter = maxitermax
    deallocate(sim2)
end subroutine isamplefilter
