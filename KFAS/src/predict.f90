! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/predict.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! auxiliary functions for computing Zalpha and ZVZ

subroutine signaltheta(tvz, zt, ahat, vt, p, n, m, theta, thetavar,d,states,m2)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, n,tvz,d,m2 !,tvh
    integer, intent(in), dimension(m2) :: states
    integer :: t
    real(dp), intent(in), dimension(p,m,(n - 1) * tvz + 1) :: zt
    real(dp), intent(in), dimension(m,n) :: ahat
    real(dp), intent(in), dimension(m,m,n) :: vt
    real(dp), intent(inout), dimension(n,p) :: theta
    real(dp), intent(inout), dimension(p,p,n) :: thetavar
    real(dp), dimension(p,m2) :: pm

    external dgemv, dsymm, dgemm

    do t = (d + 1), n
        call dgemv('n',p,m2,1.0_dp,zt(:,states,(t - 1) * tvz + 1),p,ahat(states,t),1,0.0_dp,theta(t,:),1)
        call dsymm('r','u',p,m2,1.0_dp,vt(states,states,t),m2,zt(:,states,(t - 1) * tvz + 1),p,0.0_dp,pm,p)
        call dgemm('n','t',p,p,m2,1.0_dp,pm,p,zt(:,states,(t - 1) * tvz + 1),p,0.0_dp,thetavar(:,:,t),p)
    end do

end subroutine signaltheta

subroutine zalpha(timevar, zt, alpha,theta,p,m,n,nsim,m2,states)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, n, nsim,timevar,m2
    integer, intent(in), dimension(m2) :: states
    integer :: t, i
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar + 1) :: zt
    real(dp), intent(in), dimension(n,m,nsim) :: alpha
    real(dp), intent(inout), dimension(n,p,nsim) :: theta

    external dgemv
    do i = 1,nsim
        do t = 1,n
            call dgemv('n',p,m2,1.0_dp,zt(:,states,(t - 1) * timevar + 1),p,&
            alpha(t,states,i),1,0.0_dp,theta(t,:,i),1)
        end do
    end do

end subroutine zalpha


