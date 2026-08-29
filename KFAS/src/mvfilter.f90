! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/mvfilter.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! function for computing the multivariate one step ahead prediction errors and their covariances
subroutine mvfilter(tvz, zt, p, m, n, d, at, pt, pinf, vt, ft,finf)
    use kfas_kinds, only: dp


    implicit none

    integer, intent(in) :: p, m, n,d,tvz
    integer :: t
    real(dp), intent(in), dimension(p,m,(n - 1) * tvz + 1) :: zt
    real(dp), intent(in), dimension(n,m) :: at
    real(dp), intent(inout), dimension(n,p) :: vt
    real(dp), intent(inout), dimension(p,p,n) :: ft
    real(dp), intent(inout), dimension(p,p,d) :: finf
    real(dp), intent(in), dimension(m,m,n) :: pt
    real(dp), intent(in), dimension(m,m,d) :: pinf
    real(dp), dimension(p,m) :: pm
    external dgemv, dgemm, dsymm

    do t = 1,d
        call dgemv('n',p,m,-1.0_dp,zt(:,:,(t - 1) * tvz + 1),p,at(t,:),1,1.0_dp,vt(t,:),1)
        call dsymm('r','u',p,m,1.0_dp,pt(:,:,t),m,zt(:,:,(t - 1) * tvz + 1),p,0.0_dp,pm,p)
        call dgemm('n','t',p,p,m,1.0_dp,pm,p,zt(:,:,(t - 1) * tvz + 1),p,1.0_dp,ft(:,:,t),p)

        call dsymm('r','u',p,m,1.0_dp,pinf(:,:,t),m,zt(:,:,(t - 1) * tvz + 1),p,0.0_dp,pm,p)
        call dgemm('n','t',p,p,m,1.0_dp,pm,p,zt(:,:,(t - 1) * tvz + 1),p,0.0_dp,finf(:,:,t),p)
    end do

    do t = d + 1,n
        call dgemv('n',p,m,-1.0_dp,zt(:,:,(t - 1) * tvz + 1),p,at(t,:),1,1.0_dp,vt(t,:),1)
        call dsymm('r','u',p,m,1.0_dp,pt(:,:,t),m,zt(:,:,(t - 1) * tvz + 1),p,0.0_dp,pm,p)
        call dgemm('n','t',p,p,m,1.0_dp,pm,p,zt(:,:,(t - 1) * tvz + 1),p,1.0_dp,ft(:,:,t),p)
    end do


end subroutine mvfilter
