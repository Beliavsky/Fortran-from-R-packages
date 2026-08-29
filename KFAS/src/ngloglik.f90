! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/ngloglik.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Non-Gaussian log-likelihood computation
subroutine ngloglik(yt, ymiss, timevar, zt, tt, rtv, qt, a1, p1,p1inf, p,m, &
  r, n, lik, theta, u, dist,maxiter,rankp,convtol, &
  nnd,nsim,epsplus,etaplus,aplus1,c,tol,info,antit,sim,nsim2,diff,marginal, &
  expected, htol)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, r, n,nnd,antit,nsim,sim,nsim2,&
    rankp, expected
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout) :: maxiter,marginal,info
    integer :: j,t,info2
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
    real(dp), intent(inout), dimension(m,nsim) :: aplus1
    real(dp), intent(in),dimension(nsim) :: c
    real(dp), dimension(p,p,n) :: ht
    real(dp), dimension(n,p) :: ytilde
    real(dp), intent(inout), dimension(n,p) :: theta
    real(dp), intent(inout), dimension(p,n,nsim) :: epsplus
    real(dp), intent(inout), dimension(r,n,nsim) :: etaplus
    real(dp), intent(inout) :: lik
    real(dp), dimension(n) :: tmp
    real(dp), intent(inout) :: diff
    real(dp), dimension(:), allocatable :: w
    real(dp), dimension(:,:,:), allocatable :: tsim
    real(dp), external :: ddot

    external approx, marginalxx, dpoisf, dnormf, dbinomf, dgammaf, dnbinomf, simgaussian
    ht = 0.0_dp
    !approximate
    call approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p, n, m, r,&
    theta, u, ytilde, dist,maxiter,tol,rankp,convtol,diff,lik, info, expected, htol)

    if(info /= 0 .and. info /= 3) then
        return
    end if

    if(marginal == 1) then
        j = int(sum(p1inf))
        if(j > 0) then
            call marginalxx(p1inf,zt,tt,m,p,n,j,timevar,lik,marginal)
        end if
        if(marginal == -1) then
            info = 5
            return
        end if
    end if


    do j = 1,p
        select case(dist(j))
            case(2)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dpoisf(yt(t,j), u(t,j) * exp(theta(t,j)), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(3)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dbinomf(yt(t,j), u(t,j), exp(theta(t,j)) / (1.0_dp + exp(theta(t,j))), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(4)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dgammaf(yt(t,j), u(t,j), exp(theta(t,j)) / u(t,j), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
            case(5)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dnbinomf(yt(t,j), u(t,j), exp(theta(t,j)), lik)
                        call dnormf(ytilde(t,j), theta(t,j),sqrt(ht(j,j,t)), lik)
                    end if
                end do
        end select
    end do

    if(sim == 1) then
        allocate(w(nsim))
        allocate(tsim(p, n, nsim2))
        w = 1.0_dp
        info2 = 0

        ! simulate signals
        call simgaussian(ymiss,timevar, ytilde, zt, ht, tt, rtv, qt, a1, p1, &
        p1inf, nnd,nsim, epsplus, etaplus, aplus1, p, n, m, r, info2,rankp,&
        tol,tsim,c,5,p,antit)

        if(info2 == 0) then
            ! Compute weights
            do j = 1,p
                select case(dist(j))
                    case(2) !poisson
                        tmp = exp(theta(:,j))
                        do t = 1,n
                            if(ymiss(t,j) == 0) then
                                !  do i=1,nsim2
                                w = w * exp(yt(t,j) * (tsim(j,t,:) - theta(t,j))-&
                                u(t,j) * (exp(tsim(j,t,:)) - tmp(t))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                              !  end do
                            end if
                        end do
                    case(3) !binomial
                        tmp = log(1.0_dp + exp(theta(:,j)))
                        do t = 1,n
                            if(ymiss(t,j) == 0) then
                                w = w * exp(yt(t,j) * (tsim(j,t,:) - theta(t,j))-&
                                u(t,j) * (log(1.0_dp + exp(tsim(j,t,:))) - tmp(t))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                    case(4) ! gamma
                        tmp = exp(-theta(:,j))
                        do t = 1,n
                            if(ymiss(t,j) == 0) then
                                w = w * exp(u(t,j) * (yt(t,j) * (tmp(t) - exp(-tsim(j,t,:))) + theta(t,j) - tsim(j,t,:))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                    case(5)
                        tmp = exp(theta(:,j))
                        do t = 1,n
                            if(ymiss(t,j) == 0) then
                                w = w * exp(yt(t,j) * (tsim(j,t,:) - theta(t,j)) +&
                                (yt(t,j) + u(t,j)) * log((u(t,j) + tmp(t)) / (u(t,j) + exp(tsim(j,t,:))))) / &
                                exp(-0.5_dp / ht(j,j,t) * ((ytilde(t,j) - tsim(j,t,:))**2 - (ytilde(t,j) - theta(t,j))**2))
                            end if
                        end do
                end select
            end do

            lik = lik + log(sum(w) / dble(nsim2))
        else
            info = info2
            return
        end if
        deallocate(tsim)
        deallocate(w)
    end if

end subroutine ngloglik
