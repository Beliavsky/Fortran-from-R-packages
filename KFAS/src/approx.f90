! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/approx.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Subroutine for computation of the approximating Gaussian model for non-Gaussian models

subroutine approx(yt, ymiss, timevar, zt, tt, rtv, ht, qt, a1, p1,p1inf, p, n, m, r,&
theta, u, ytilde, dist, maxiter, tol, rankp, convtol, diff, lik, info, expected, htol)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p,m, r, n,rankp, expected
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    integer, intent(in), dimension(p) :: dist
    integer, intent(inout) :: maxiter,info
    integer :: i, k,tvrqr,kk,jt,dt
    real(dp), intent(in) :: tol,convtol
    real(dp), intent(inout) :: htol
    real(dp), intent(in), dimension(n,p) :: u
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rtv
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m) :: a1
    real(dp), intent(in), dimension(m,m) :: p1,p1inf
    real(dp), intent(inout), dimension(n,p) :: theta
    real(dp), intent(inout), dimension(n,p) :: ytilde
    real(dp), intent(inout), dimension(p,p,n) :: ht
    real(dp), intent(inout) :: diff
    real(dp), dimension(m,r) :: mr
    real(dp), dimension(m,m,(n - 1) * max(timevar(4),timevar(5)) + 1) :: rqr
    real(dp), intent(inout) :: lik
    real(dp), external :: ddot
    integer, external :: finitex
    real(dp) dev, devold
    real(dp), dimension(n,p) :: thetanew, thetaold
    real(dp), dimension(p,n) :: ft,finf
    real(dp), dimension(m,p,n) :: kt,kinf

    external dgemm, pytheta, pthetafirst, approxloop, pthetarest

    !compute rqr
    tvrqr = max(timevar(4),timevar(5))
    do i = 1, (n - 1) * tvrqr + 1
        call dgemm('n','n',m,r,r,1.0_dp,rtv(:,:,(i - 1) * timevar(4) + 1),m,qt(:,:,(i - 1) * timevar(5) + 1),r,0.0_dp,mr,m)
        call dgemm('n','t',m,m,r,1.0_dp,mr,m,rtv(:,:,(i - 1) * timevar(4) + 1),m,0.0_dp,rqr(:,:,i),m)
    end do

    ! compute logp(theta) for the first time, no need to compute kt/kinf and ft/finf in successive calls
    ! in case of totally diffuse initialization term p(theta)=0 for all theta

    if(rankp /= m) then
        call pthetafirst(theta, timevar, zt, tt, rqr, a1, p1, p1inf, p, m, n, devold, tol,rankp,kt,kinf,ft,finf,dt,jt)
    end if
    thetaold = theta
    devold = -huge(devold)
    k = 0
    do while(k < maxiter)

        k = k + 1
        ! compute new guess thetanew
        call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
        theta, thetanew, u, ytilde, dist,tol,rankp,lik, expected)
        ! and log(p(theta|y))
        call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)
        if(rankp /= m) then
            call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, kt,kinf,ft,finf,dt,jt)
        end if
        !non-finite value in linear predictor or muhat
        if(finitex(sum(thetanew, MASK = ymiss == 0)) == 0 .or. &
           finitex(maxval(exp(thetanew), MASK = ymiss == 0)) == 0) then
            if(k > 1) then
                kk = 0
                do while(finitex(sum(thetanew, MASK = ymiss == 0)) == 0 .or. &
                         finitex(maxval(exp(thetanew), MASK = ymiss == 0)) == 0)
                    kk = kk + 1
                    if(kk > maxiter) then
                        info = 1
                        maxiter = k
                        return
                    end if
                    !backtrack
                    theta = 0.5_dp * (thetaold + theta)
                    call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                    theta, thetanew, u, ytilde, dist,tol,rankp,lik, expected)

                    call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)
                    if(rankp /= m) then
                        call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, kt,kinf,ft,finf,dt,jt)
                    end if
                end do
            else !cannot correct step size as we have just began
                info = 1
                maxiter = k
                return
            end if
        end if

        if(finitex(dev) == 0) then !non-finite value of objective function
            if(k > 1) then
                kk = 0
                do while(finitex(dev) == 0)
                    kk = kk + 1
                    if(kk > maxiter) then !did not find valid likelihood
                        info = 2
                        maxiter = k
                        return
                    end if

                    theta = 0.5_dp * (thetaold + theta)
                    call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                    theta, thetanew, u, ytilde, dist,tol,rankp,lik, expected)

                    call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)

                    if(rankp /= m) then
                        call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, kt,kinf,ft,finf,dt,jt)
                    end if
                end do
            else !cannot correct step size as we have just began
                info = 2
                maxiter = k
                return
            end if
        end if


        ! decreasing deviance
        if((dev - devold) / (0.1_dp + abs(dev)) < -convtol .and. k > 1) then
            kk = 0
            do while((dev - devold) / (0.1_dp + abs(dev)) < convtol .and. kk < maxiter)
                kk = kk + 1
                ! previous theta produced too 'big' thetanew
                ! new guess by halving the last try
                theta = 0.5_dp * (thetaold + theta)
                call approxloop(yt, ymiss, timevar, zt, tt, rtv, ht, qt, rqr, tvrqr, a1, p1,p1inf, p,n,m,r, &
                theta, thetanew, u, ytilde, dist,tol,rankp,lik, expected)

                call pytheta(thetanew, dist, u, yt, ymiss, dev, p, n)
                if(rankp /= m) then
                    call pthetarest(thetanew, timevar, zt, tt, a1, p, m, n, dev, kt,kinf,ft,finf,dt,jt)
                end if

            end do
        end if

        diff = abs(dev - devold) / (0.1_dp + abs(dev))
        if(diff < convtol) then !convergence
            theta = thetanew
            info = 0
            exit
        else
            thetaold = theta
            theta = thetanew
            devold = dev
        end if
    end do
    if(maxiter == k) then
        info = 3
    end if
    maxiter = k

    ! check if the the approximation lead to potentially zero signal-to-noise ratio
    if(maxval(ht) > htol) then
      info = -5
      htol = maxval(ht)
      return
    end if

end subroutine approx
