! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/isample.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Importance sampling of non-gaussian model

subroutine isample(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, u, dist, &
p, n, m, r, theta, maxiter,rankp,convtol, nnd,nsim,epsplus,etaplus,&
aplus1,c,tol,info,antithetics,w,sim,simwhat,simdim, expected, htol)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p,m, r, n,nnd,antithetics,nsim,simwhat,simdim&
    , rankp, expected
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,info
    integer :: t, j,i,info2
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
    real(dp), dimension(n,p) :: ytilde
    real(dp), dimension(n) :: tmp
    real(dp), intent(inout), dimension(3 * nsim * antithetics + nsim) :: w
    real(dp) :: diff
    real(dp), external :: ddot
    real(dp) :: lik
    real(dp), dimension(:,:), allocatable :: tsim

    external approx, simgaussian

    ht = 0.0_dp

    ! approximate
    call approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p,n,m,r,&
    theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff,lik,info, expected, htol)

    if(info /= 0 .and. info /= 3) then
        return
    end if

    info2 = 0
    ! simulate signals
    call simgaussian(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
    p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info2,rankp,&
    tol,sim,c,simwhat,simdim,antithetics)

    if(info2 /= 0) then
        info = info2
        return
    end if

    ! compute importance weights

    w = 1.0_dp

    if(simwhat == 5) then
        do j = 1,p
            select case(dist(j))
                case(2) !poisson
                    tmp = exp(theta(:,j))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then

                            w = w * exp(yt(t,j) * (sim(j,t,:) - theta(t,j))-&
                            u(t,j) * (exp(sim(j,t,:)) - tmp(t))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                        end if
                    end do
                case(3) !binomial
                    tmp = log(1.0_dp + exp(theta(:,j)))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then

                            w = w * exp(yt(t,j) * (sim(j,t,:) - theta(t,j))-&
                            u(t,j) * (log(1.0_dp + exp(sim(j,t,:))) - tmp(t))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                        end if
                    end do
                case(4) ! gamma
                    tmp = exp(-theta(:,j))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then
                            w = w * exp(u(t,j) * (yt(t,j) * (tmp(t) - exp(-sim(j,t,:))) + theta(t,j) - sim(j,t,:))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                        end if
                    end do
                case(5) !negbin
                    tmp = exp(theta(:,j))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then
                            w = w * exp(yt(t,j) * (sim(j,t,:) - theta(t,j)) +&
                            (yt(t,j) + u(t,j)) * log((u(t,j) + tmp(t)) / (u(t,j) + exp(sim(j,t,:))))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - sim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                        end if
                    end do
            end select
        end do

    else

        allocate(tsim(p,(3 * nsim * antithetics + nsim) * (5 - simwhat)))
        do j = 1,p
            select case(dist(j))
                case(2) !poisson
                    tmp = exp(theta(:,j))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then
                            do i = 1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim(:,t,i),1)
                            end do
                            w = w * exp(yt(t,j) * (tsim(j,:) - theta(t,j))-&
                            u(t,j) * (exp(tsim(j,:)) - tmp(t))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                        end if
                    end do
                case(3) !binomial
                    tmp = log(1.0_dp + exp(theta(:,j)))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then
                            do i = 1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim(:,t,i),1)
                            end do
                            w = w * exp(yt(t,j) * (tsim(j,:) - theta(t,j))-&
                            u(t,j) * (log(1.0_dp + exp(tsim(j,:))) - tmp(t))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))

                        end if
                    end do
                case(4) ! gamma
                    tmp = exp(-theta(:,j))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then
                            do i = 1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim(:,t,i),1)
                            end do
                            w = w * exp(u(t,j) * (yt(t,j) * (tmp(t) - exp(-tsim(j,:))) + theta(t,j) - tsim(j,:))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                        end if
                    end do
                case(5) !negbin
                    tmp = exp(theta(:,j))
                    do t = 1,n
                        if(ymiss(t,j) == 0) then
                            do i = 1,3 * nsim * antithetics + nsim
                                tsim(j,i) = ddot(m,zt(j,:,(t - 1) * timevar(1) + 1),1,sim(:,t,i),1)
                            end do
                            w = w * exp(yt(t,j) * (tsim(j,:) - theta(t,j)) +&
                            (yt(t,j) + u(t,j)) * log((u(t,j) + tmp(t)) / (u(t,j) + exp(tsim(j,:))))) / &
                            exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                        end if
                    end do
            end select
        end do
        deallocate(tsim)
    end if

end subroutine isample
